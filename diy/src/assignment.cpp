//---------------------------------------------------------------------------
//
// assignment class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2011 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//----------------------------------------------------------------------------

#include "assignment.hpp"

//----------------------------------------------------------------------------
//
// RoundRobinAssignment derived class
//
//----------------------------------------------------------------------------
//
// constructor
//
// start_b: starting block global id (number of blocks in prior domains)
// tot_b: total number of blocks in this domain
// nb: number of blocks in my process (output)
// max_b: maximum number of blocks in any process (output)
// comm: MPI communicator
//
RoundRobinAssignment::RoundRobinAssignment(int start_b, int tot_b, int &nb, 
					   int &max_b, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->start_b = start_b;
  this->tot_b = tot_b;
  this->comm = comm;

  nb = tot_b / groupsize + (rank < tot_b % groupsize ? 1 : 0);
  this->nb = nb;
  max_b = (int)(ceil(tot_b / (float)groupsize));

}
//----------------------------------------------------------------------------
//
// ProcOrderAssignment derived class
//
//----------------------------------------------------------------------------
//
// constructor
//
// start_b: starting block global id (number of blocks in prior domains)
// tot_b: total number of blocks in this domain
// nb: number of blocks in my process (output)
// max_b: maximum number of blocks in any process (output)
// comm: MPI communicator
//
ProcOrderAssignment::ProcOrderAssignment(int start_b, int tot_b, int &nb, 
					 int &max_b, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->start_b = start_b;
  this->tot_b = tot_b;
  this->comm = comm;

  if (rank < groupsize - 1)
    nb = tot_b / groupsize;
  else // last rank gets the remaining blocks
    nb = tot_b - rank * tot_b / groupsize;
  this->nb = nb;
  max_b = (int)(ceil(tot_b / (float)groupsize));

}
//----------------------------------------------------------------------------
//
// ExistingAssignment derived class
//
//----------------------------------------------------------------------------
//
// constructor
//
// start_b: starting block global id (number of blocks in prior domains)
// nb: number of blocks in my process
// max_b: maximum number of blocks in any process (output)
// tot_b: total number of blocks in this domain (output)
// comm: MPI communicator
//
ExistingAssignment::ExistingAssignment(int start_b, int nb, int &max_b, 
				       int &tot_b, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->start_b = start_b;
  this->nb = nb;
  this->comm = comm;

  // communicate to find max_b and tot_b
  // would be faster with one allreduce and a custom reduction op
  MPI_Allreduce(&nb, &max_b, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&nb, &tot_b, 1, MPI_INT, MPI_SUM, comm);

}
//----------------------------------------------------------------------------
//
// reassignment functions based on Zoltan library
// not part of a class, so that various functions can be used as callbacks
// todo: see if some of these can be included in the assignment class
//
//---------------------------------------------------------------------------

// todo: replace Lattice with diy

#ifdef ZOLTAN

static int *NumSeeds; // number of seeds this process has in each block
static int NumBlks; // number of blocks this process has
static int time_group; // current time group
static Zoltan_Struct *zz;
static int Wgts[MAX_NUM_BLOCKS]; // weight of each of my blocks

//---------------------------------------------------------------------------
//
// inits zoltan repartitioning for use by Lattice4D
//
// lat: pointer to lattice
// comm: MPI communicator
//
void InitRepartition4D(void *lat, MPI_Comm comm) {

  float ver; // zoltan version number

  int retval = Zoltan_Initialize(1, NULL, &ver);
  assert(retval == ZOLTAN_OK);

  // zoltan structure
  zz = Zoltan_Create(comm);

  // general params
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");

  // callback functions
  Zoltan_Set_Num_Obj_Fn(zz, GetNumberofAssignedObjects4D, lat);
  Zoltan_Set_Obj_List_Fn(zz, GetObjectList4D, lat);
  Zoltan_Set_Num_Geom_Fn(zz, GetObjectSize4D, lat);
  Zoltan_Set_Geom_Multi_Fn(zz, GetObjects4D, lat);

}
//-----------------------------------------------------------------------
//
// inits zoltan repartitioning for use by LatticeAMR
//
// lat: pointer to lattice
// comm: MPI communicator
//
void InitRepartitionAMR(void *lat, MPI_Comm comm) {

  float ver; // zoltan version number

  int retval = Zoltan_Initialize(1, NULL, &ver);
  assert(retval == ZOLTAN_OK);

  // zoltan structure
  zz = Zoltan_Create(comm);

  // general params
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");

  // callback functions
  Zoltan_Set_Num_Obj_Fn(zz, GetNumberofAssignedObjectsAMR, lat);
  Zoltan_Set_Obj_List_Fn(zz, GetObjectListAMR, lat);
  Zoltan_Set_Num_Geom_Fn(zz, GetObjectSizeAMR, lat);
  Zoltan_Set_Geom_Multi_Fn(zz, GetObjectsAMR, lat);

}
//-----------------------------------------------------------------------
//
// Computes a new partition and implements it
//
// grp: current time group
// nb: number of local blocks
// block_ranks: global ranks of local blocks
// neighbor_ranks: global ranks of neighbors of local blocks
// neighbor_procs: process ids of neighbors of local blocks
// part: partition data structure
// seeds: locations to store received points, indexed by local block number
// size_seeds: sizes of seed arrays (will be grown automatically if necessary)
// num_seeds: number of seeds stored for each block
// avg_neigh: average number of neighbors per each of my blocks
// alloc_blocks: number of local blocks allocated
// alloc_neighbors: number of neighbors allocated for each local block
// comm: MPI communicator
// add_neighbor: pointer to a function that adds a neighbor
// wgts: weights of blocks
//
// pass NULL for seeds, size_seeds, and num_seeds
// if no seeds have been assigned yet
//
void ChangePartition(int grp, int *nb, int **block_ranks, int ***neighbor_ranks,
		     int ***neighbor_procs, Partition *part, VECTOR4 ***seeds, 
		     int **size_seeds, int **num_seeds, int *avg_neigh,
		     int *alloc_blocks, int **alloc_neighbors, MPI_Comm comm,
		     OSUFlow ***osuflow,
		     void (*add_neighbor)(int, int, int, int *, int *, int ***,
					  int ***, Partition *),
		     int *wgts) {

  int changes; // whether zoltan changed to partition
  int ngid_items; // number of items per global id
  int nlid_items; // number of items per local id
  int nimport; // number of items to be imported
  int nexport; // number of items to be exported
  ZOLTAN_ID_PTR import_gids; // global ids being imported
  ZOLTAN_ID_PTR import_lids; // local ids being imported
  ZOLTAN_ID_PTR export_gids; // global ids being exported
  ZOLTAN_ID_PTR export_lids; // local ids being exported
  int *import_procs; // source processors of objects being imported
  int *import_to_part; // parts where imported objects should go
  int *export_procs; // source processors of objects being exported
  int *export_to_part; // parts where exported objects should go
  static int alloc_nb_seeds = 0; // blocks alloc'd to hold number of seeds
  int *n_msg; // neighbors message buffer
  float *s_msg; // seeds message buffer
  MPI_Status status;
  int myproc;
  int i, j, k, m, n;

  time_group = grp;

  MPI_Comm_rank(comm, &myproc);

  // save weights
  assert(*nb <= MAX_NUM_BLOCKS);
  if (wgts != NULL) {
    for (i = 0; i < *nb; i++)
      Wgts[i] = wgts[i];
  }
  else {
    for (i = 0; i < *nb; i++)
      Wgts[i] = 0;
  }

  // todo: determine the max size of the message in a better way
  n_msg = (int *)malloc((2 * MAX_NUM_NEIGHBORS + 1) * sizeof(int));
  s_msg = (float *)malloc((1 + 4 * MAX_NUM_SEEDS) * sizeof(float));

  // save the number of seeds for weighting the partition
  if (alloc_nb_seeds == 0) {
    NumSeeds = (int *)malloc(sizeof(int));
    alloc_nb_seeds = 1;
  }
  while (alloc_nb_seeds < *nb) {
    NumSeeds = (int *)realloc(NumSeeds, alloc_nb_seeds * sizeof(int) * 2);
    alloc_nb_seeds *= 2;
  }
  for (i = 0; i < *nb; i++) {
    if (seeds != NULL)
      NumSeeds[i] = (*num_seeds)[i];
    else
      NumSeeds[i] = 0;
  }
  NumBlks = *nb;

  // compute the partition
  int retval = Zoltan_LB_Partition (zz, &changes, &ngid_items, &nlid_items, 
				    &nimport, &import_gids, &import_lids, 
				    &import_procs, &import_to_part, &nexport, 
				    &export_gids, &export_lids, &export_procs, 
				    &export_to_part);
  assert(retval == 0);

#ifdef DEBUG
  for (i = 0; i < nexport; i++)
    fprintf(stderr,"rank %d to proc %d; ",export_gids[i], export_procs[i]);
  fprintf(stderr,"\n");
#endif

  // exchange process assignments
  ExchangeExports(nexport, export_gids, export_lids, export_procs, *nb, 
		  *block_ranks, *neighbor_ranks, *neighbor_procs, part, comm);

  // send outgoing neighbor lists
  for (i = 0; i < nexport; i++) {

    // find local block number of the global partition id
    for (k = 0; k < *nb; k++) {
      if ((*block_ranks)[k] == (int)(export_gids[i]))
	break;
    }
    assert(k < *nb); // sanity check

    // pack and send neighbors info
    n_msg[0] = part->parts[export_gids[i]].NumNeighbors;
    n = 1;
    assert(n_msg[0] < MAX_NUM_NEIGHBORS); // msg has space?
    for (j = 0; j < n_msg[0]; j++) {
      n_msg[n++] = (*neighbor_ranks)[k][j];
      n_msg[n++] = (*neighbor_procs)[k][j];
    }
    MPI_Send(n_msg, 2 * n_msg[0] + 1, MPI_INT, export_procs[i], 0, comm);

    // pack and send seed info
    n = 0;
    if (seeds != NULL) {
      assert((*num_seeds)[k] < MAX_NUM_SEEDS); // msg has space?
      s_msg[n++] = (float)(*num_seeds)[k];
      for (j = 0; j < (*num_seeds)[k]; j++) {
	for (m = 0; m < 4; m++)
	  s_msg[n++] = (*seeds)[k][j][m];
      }
    }
    else
      s_msg[0] = 0.0f;
    MPI_Send(s_msg, n, MPI_FLOAT, export_procs[i], 0, comm);

    // remove the block
    RemoveBlock(export_gids[i], export_procs[i], nb, block_ranks, 
		alloc_neighbors, neighbor_ranks, neighbor_procs, part,
		seeds, size_seeds, num_seeds);

  }

  // receive incoming neighbor lists
  for (i = 0; i < nimport; i++) {

    // neighbors message
    MPI_Recv(n_msg, 2 * MAX_NUM_NEIGHBORS + 1, MPI_INT, import_procs[i], 
	     0, comm, &status);

    // add the block
    AddBlock(import_gids[i], n_msg[0], &(n_msg[1]), alloc_blocks, nb, 
	     block_ranks, alloc_neighbors, neighbor_ranks, neighbor_procs,
	     part, myproc, seeds, size_seeds, num_seeds, osuflow, add_neighbor);

    // seeds message
    MPI_Recv(s_msg, 1 + 4 * MAX_NUM_SEEDS, MPI_FLOAT, import_procs[i], 0, comm,
	     &status);

    // find local block number of the global partition id
    for (k = 0; k < *nb; k++) {
      if ((*block_ranks)[k] == (int)(import_gids[i]))
	break;
    }
    assert(k < *nb); // sanity check

    // unpack and copy seeds
    if (seeds != NULL) {
      n = 0;
      (*num_seeds)[k] = (int)s_msg[n++];
      // grow the size of seeds if necessary
      while ((*size_seeds)[k] < (*num_seeds)[k] * (int)(sizeof(VECTOR4))) {
	(*seeds)[k] = (VECTOR4 *)realloc((*seeds)[k], (*size_seeds)[k] * 2);
	(*size_seeds)[k] *= 2;
      }
      for (j = 0; j < (*num_seeds)[k]; j++) {
	(*seeds)[k][j].Set(s_msg[n], s_msg[n + 1], s_msg[n + 2], s_msg[n + 3]);
	n += 4;
      }

    }
  }

  // update perf stats
  int tot_neighbors = 0;
  for (i = 0; i < *nb; i++)
    tot_neighbors += part->parts[(*block_ranks)[i]].NumNeighbors;
  *avg_neigh = (*nb > 0 ? tot_neighbors / *nb : 0);

  // cleanup
//   Zoltan_LB_Free_Part(&import_gids, &import_lids, 
//                       &import_procs, &import_to_part);
//   Zoltan_LB_Free_Part(&export_gids, &export_lids, 
//                       &export_procs, &export_to_part);

}
//---------------------------------------------------------------------------
//
// exchanges exported assignments with all neighbors
//
// nexport: number of imported blocks
// export_gids: global ranks of exported blocks
// export_procs: new process assignments of exported blocks
// nb: number of local blocks
// block_ranks: ranks of local blocks
// neighbor_ranks: ranks of neighbors of local blocks
// neighbor_procs: process assignments of local blocks
// part: partition data structure
// commm: MPI communicator
//
void ExchangeExports(int nexport, ZOLTAN_ID_PTR export_gids, 
		     ZOLTAN_ID_PTR export_lids, int *export_procs,
		     int nb, int *block_ranks, int **neighbor_ranks, 
		     int **neighbor_procs, Partition *part, MPI_Comm comm) {

  int *Send, *Recv; // message buffers
  int *SendSizes, *RecvSizes; // sizes index into message buffers
  int *SendDispls, *RecvDispls; // displacements index into count information
  int p; // process number
  int nproc; // total number of processes
  int ns, nr; // number of send, receive elements
  int i, j, k, n;

  MPI_Comm_size(comm, &nproc);

  SendSizes = (int *)malloc(nproc * sizeof(int));
  SendDispls = (int *)malloc(nproc * sizeof(int));
  RecvSizes = (int *)malloc(nproc * sizeof(int));
  RecvDispls = (int *)malloc(nproc * sizeof(int));

  // debug
  int myproc;
  MPI_Comm_rank(comm, &myproc);

  // exchange number of exports going to each process
  for (p = 0; p < nproc; p++) { // all processes
    SendSizes[p] = 0;
    for (k = 0; k < nexport; k++) { // my exports
      i = export_lids[k];
      // neighbors of exported block
      for (j = 0; j < part->parts[block_ranks[i]].NumNeighbors; j++) {
	if (neighbor_procs[i][j] == p) {
	  SendSizes[p]++;
	  break;
	}
      }				  
    }
  }

  MPI_Alltoall(SendSizes, 1, MPI_INT, RecvSizes, 1, MPI_INT, comm);

  // index and message vectors for alltoallv
  for (i = 0; i < nproc; i++) {
    SendSizes[i] *= 2;
    RecvSizes[i] *= 2;
  }
  SendDispls[0] = RecvDispls[0] = 0;
  for (i = 1; i < nproc; i++) {
    SendDispls[i] = SendDispls[i - 1] + SendSizes[i - 1];
    RecvDispls[i] = RecvDispls[i - 1] + RecvSizes[i - 1];
  }
  ns = nr = 0;
  for (i = 0; i < nproc; i++) {
    ns += SendSizes[i];
    nr += RecvSizes[i];
  }
  if (ns > 0)
    Send = (int *)malloc(ns * sizeof(int));
  else // malloc at least 1 int to quite valgrind unitialized error
    Send = (int *)malloc(sizeof(int));
  if (nr > 0)
    Recv = (int *)malloc(nr * sizeof(int));
  else // malloc at least 1 int to quite valgrind unitialized error
    Recv = (int *)malloc(sizeof(int));

  // exchange exports going to each process
  n = 0;
  for (p = 0; p < nproc; p++) { // all processes
    for (k = 0; k < nexport; k++) { // my exports
      i = export_lids[k];
      // neighbors of exported block
      for (j = 0; j < part->parts[block_ranks[i]].NumNeighbors; j++) {
	if (neighbor_procs[i][j] == p) {
	  Send[n++] = export_gids[k];
	  Send[n++] = export_procs[k];
	  break;
	}
      }				  
    }
  }

  MPI_Alltoallv(Send, SendSizes, SendDispls, MPI_INT, 
		Recv, RecvSizes, RecvDispls, MPI_INT, comm);

  // process the received info: revise the process assignments of my neighbors
  n = 0;
  for (p = 0; p < nproc; p++) { // all processes
    for(k = RecvDispls[p]; k < RecvDispls[p] + RecvSizes[p]; k += 2) {
      for (i = 0; i < nb; i++) { // my blocks
	for (j = 0; j < part->parts[block_ranks[i]].NumNeighbors; j++) {
	  if (neighbor_ranks[i][j] == Recv[k])
	    neighbor_procs[i][j] = Recv[k + 1];
	}				  
      }
    }
  }

  // cleanup
  if (ns > 0)
    free(Send);
  if (nr > 0)
    free(Recv);
  free(SendSizes);
  free(RecvSizes);
  free(SendDispls);
  free(RecvDispls);

}
//------------------------------------------------------------------------
//
// removes a block from my process and assigns it to another process
//
// myrank: global partition rank of the block
// newproc: new process id where the block is heading
// nb: number of local blocks
// block_ranks: global ranks of local blocks
// alloc_neighbors: number of neighbors allocated for each local block
// neighbor_ranks: global ranks of neighbors of local blocks
// neighbor_procs: process ids of neighbors of local blocks
// part: partition data structure
// seeds: seed points in each block
// size_seeds: allocated size (bytes) of seeds for each block
// num_seeds: number of seeds in each block
// add_neighbor: pointer to a function that adds a neighbor
//
// pass NULL for seeds, size_seeds, and num_seeds
// if no seeds have been assigned yet
//
void RemoveBlock(int myrank, int newproc, int *nb, int **block_ranks,
		 int **alloc_neighbors, int ***neighbor_ranks,
		 int ***neighbor_procs, Partition *part,
		 VECTOR4 ***seeds, int **size_seeds, int **num_seeds) {

  int i, j, k;

  // j will be the local block number that is being removed
  for (j = 0; j < *nb; j++) {
    if ((*block_ranks)[j] == myrank)
      break;
  }
  assert(j < *nb); // sanity check; should always be true

  // block_ranks: move subsequent blocks up by one
  for (k = j; k < *nb - 1; k++)
    (*block_ranks)[k] = (*block_ranks)[k + 1];

  // neighbor_ranks, neighbor_procs: move subsequent blocks up by one
  for (k = j; k < *nb - 1; k++) {
    for (i = 0; i < part->parts[(*block_ranks)[k]].NumNeighbors; i++) {

      // the preceding entry may not be allocated large enough
      if ((*alloc_neighbors)[k] < (*alloc_neighbors)[k + 1]) {
	if ((*alloc_neighbors)[k] == 0) {
	  (*neighbor_ranks)[k] = (int *)malloc((*alloc_neighbors)[k + 1] * 
				sizeof(int));
	  (*neighbor_procs)[k] = (int *)malloc((*alloc_neighbors)[k + 1] * 
				sizeof(int));
	}
	else {
	  (*neighbor_ranks)[k] = 
		  (int *)realloc((*neighbor_ranks)[k],
				 (*alloc_neighbors)[k + 1] * 
				 sizeof(int));
	  (*neighbor_procs)[k] = 
		  (int *)realloc((*neighbor_procs)[k],
				 (*alloc_neighbors)[k + 1] * 
				 sizeof(int));
	}
	(*alloc_neighbors)[k] = (*alloc_neighbors)[k + 1];
      }

      (*neighbor_ranks)[k][i] = (*neighbor_ranks)[k + 1][i];
      (*neighbor_procs)[k][i] = (*neighbor_procs)[k + 1][i];

    }
  }

  // seeds: remove seeds and move subsequent blocks up by one
  if (num_seeds != NULL) {
    (*num_seeds)[j] = 0;
    for (k = j; k < *nb - 1; k++) {
      // grow the size of seeds if necessary
      while ((*size_seeds)[k] < (*num_seeds)[k + 1] * (int)(sizeof(VECTOR4))) {
	(*seeds)[k] = (VECTOR4 *)realloc((*seeds)[k], (*size_seeds)[k] * 2));
	(*size_seeds)[k] *= 2;
      }
      (*num_seeds)[k] = (*num_seeds)[k + 1];
      for (i = 0; i < (*num_seeds)[k]; i++)
	(*seeds)[k][i] = (*seeds)[k + 1][i];
    }
  }

  // update the procs of other blocks in neighbor_procs
  for (k = 0; k < *nb - 1; k++) {
    for (i = 0; i < part->parts[(*block_ranks)[k]].NumNeighbors; i++) {
      if ((*neighbor_ranks)[k][i] == myrank)
	(*neighbor_procs)[k][i] = newproc;
    }
  }

  part->RemoveBlock(myrank); // update partition object

  (*nb)--; // my number of blocks

}
//---------------------------------------------------------------------------
//
// adds a block to my process
//
// myrank: global partition rank of the block being added
// num_neighbors: number of neighbors the block will have
// neighbors: list of (neighbor_rank, neighbor_proc pairs)
// alloc_blocks: number of local blocks allocated
// nb: number of local blocks
// block_ranks: global ranks of local blocks
// alloc_neighbors: number of neighbors allocated for each local block
// neighbor_ranks: global ranks of neighbors of local blocks
// neighbor_procs: process ids of neighbors of local blocks
// part: partition data structure
// myproc: my process id
// seeds: seed points in each block
// size_seeds: allocated size (bytes) of seeds for each block
// num_seeds: number of seeds in each block
// add_neighbor: pointer to a function that adds a neighbor
//
// pass NULL for seeds, size_seeds, and num_seeds
// if no seeds have been assigned yet
//
void AddBlock(int myrank, int num_neighbors, int *neighbors, int *alloc_blocks,
	      int *nb, int **block_ranks, int **alloc_neighbors, 
	      int ***neighbor_ranks, int ***neighbor_procs, Partition *part,
	      int myproc, VECTOR4 ***seeds, int **size_seeds, int **num_seeds,
	      OSUFlow ***osuflow,
	      void (*add_neighbor)(int, int, int, int *, int *, 
				   int ***, int ***, Partition *)) {

  int old_n = *alloc_blocks;

  int i, j, k, n;

  if (old_n < *nb + 1) {

    if (old_n == 0) {
      n = 1;
      *block_ranks = (int *)malloc(n * sizeof(int));
      *alloc_neighbors = (int *)malloc(n * sizeof(int));
      *neighbor_ranks = (int **)malloc(n * sizeof(int *));
      *neighbor_procs = (int **)malloc(n * sizeof(int *));
      *size_seeds = (int *)malloc(n * sizeof(int));
      *num_seeds = (int *)malloc(n * sizeof(int));
      *seeds = (VECTOR4 **)malloc(n * sizeof(VECTOR4 *));
      *osuflow = (OSUFlow **)malloc(n * sizeof(OSUFlow));
    }
    else {
      n = old_n * 2;
      *block_ranks = (int *)realloc(*block_ranks, n * sizeof(int));
      *alloc_neighbors = (int *)realloc(*alloc_neighbors, n * sizeof(int));
      *neighbor_ranks = (int **)realloc(*neighbor_ranks, n * sizeof(int *));
      *neighbor_procs = (int **)realloc(*neighbor_procs, n * sizeof(int *));
      *size_seeds = (int *)realloc(*size_seeds, n * sizeof(int));
      *num_seeds = (int *)realloc(*num_seeds, n * sizeof(int));
      *seeds = (VECTOR4 **)realloc(*seeds, n  * sizeof(VECTOR4 *));
      *osuflow = (OSUFlow **)realloc(*osuflow, n * sizeof(OSUFlow));

    }

    *alloc_blocks = n;
    for(i = old_n; i < n; i++) {
      (*alloc_neighbors)[i] = 0;
      (*seeds)[i] = (VECTOR4 *)malloc(sizeof(VECTOR4));
      (*num_seeds)[i] = 0;
      (*size_seeds)[i] = sizeof(VECTOR4);
      (*osuflow)[i] = new OSUFlow;
    }

  }

  part->AddBlock(myrank); // update partition object

  (*block_ranks)[*nb] = myrank;

  // update neighbors
  n = 0;
  for (j = 0; j < num_neighbors; j++) {
    add_neighbor(*nb, neighbors[n], neighbors[n + 1], *block_ranks, 
		 *alloc_neighbors, neighbor_ranks, neighbor_procs, part);
    n += 2;
  }

  (*nb)++;

  // update the procs of other blocks in neighbor_procs
  for (k = 0; k < *nb; k++) {
    for (i = 0; i < part->parts[(*block_ranks)[k]].NumNeighbors; i++) {
      // see if I own this neighbor block
      for (j = 0; j < *nb; j++) {
	if ((*block_ranks)[j] == (*neighbor_ranks)[k][i])
	  break;
      }
      // if so, sync neighbor_procs
      if (j < *nb)
	(*neighbor_procs)[k][i] = myproc;
    }
  }

}
//---------------------------------------------------------------------------
//
// Callback functions for zoltan to use in conjunction w/ Lattice4D
//
//-----------------------------------------------------------------------
//
// returns number of objects assigned to my process
//
int GetNumberofAssignedObjects4D(void *lat, int *err) {

  int n = 0;
  int i;
  *err = 0;
  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroup4D(lat, time_group, i))
      n++;
  }
  return n;

}
//-----------------------------------------------------------------------
//
// sets the global and local identifiers of the objects assigned to my process
// also sets the weight of each object (optionally)
//
void GetObjectList4D(void *lat, int ngids, int nlids, 
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
		   int wgt_dim, float *obj_wgts, int *err) {

  int i, n = 0;

  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroup4D(lat, time_group, i)) {
      gids[n] = ((Lattice4D *)lat)->GetRank(i);
      lids[n] = i;
      obj_wgts[n] = Wgts[i];
      n++;
    }
  }

  *err = 0;

}
//-----------------------------------------------------------------------
//
// returns the dimensionality of an object
//
int GetObjectSize4D(void *lat, int *err) {

  *err = 0;
  return 3;

}
//-----------------------------------------------------------------------
//
// sets the geometry of the objects assigned to my process
// using the minimum corner of a block as the geometry (for now)
//
void GetObjects4D(void *lat, int ngids, int nlids, int nobjs, 
		ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int ndim, 
		double *pts, int *err) {

  float from[3], to[3]; // block spatial extent
  int min_t, max_t; // block temporal extent
  int i, n = 0;

  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroup4D(lat, time_group, i)) {
      ((Lattice4D *)lat)->GetGlobalVB(((Lattice4D *)lat)->GetRank(i), 
				      from, to, &min_t, &max_t);
      pts[n++] = (double)((from[0] + to[0]) * 0.5);
      pts[n++] = (double)((from[1] + to[1]) * 0.5);
      pts[n++] = (double)((from[2] + to[2]) * 0.5);
    }
  }

  *err = 0;

}
//-----------------------------------------------------------------------
//
// Callback functions for zoltan to use in conjunction w/ LatticeAMR
//
//-----------------------------------------------------------------------
//
// returns number of objects assigned to my process
//
int GetNumberofAssignedObjectsAMR(void *lat, int *err) {

  int n = 0;
  int i;
  *err = 0;
  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroupAMR(lat, time_group, i))
      n++;
  }
  return n;

}
//-----------------------------------------------------------------------
//
// sets the global and local identifiers of the objects assigned to my process
// also sets the weight of each object (optionally)
//
void GetObjectListAMR(void *lat, int ngids, int nlids, 
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
		   int wgt_dim, float *obj_wgts, int *err) {

  int i, n = 0;

  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroupAMR(lat, time_group, i)) {
      gids[n] = ((LatticeAMR *)lat)->GetRank(i);
      lids[n] = i;
      obj_wgts[n] = Wgts[i];
      n++;
    }
  }

  *err = 0;

}
//-----------------------------------------------------------------------
//
// returns the dimensionality of an object
//
int GetObjectSizeAMR(void *lat, int *err) {

  *err = 0;
  return 3;

}
//-----------------------------------------------------------------------
//
// sets the geometry of the objects assigned to my process
// using the minimum corner of a block as the geometry (for now)
//
void GetObjectsAMR(void *lat, int ngids, int nlids, int nobjs, 
		ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int ndim, 
		double *pts, int *err) {

  float from[3], to[3]; // block spatial extent
  int min_t, max_t; // block temporal extent
  int i, n = 0;

  for (i = 0; i < NumBlks; i++) {
    if (IsBlockInTimeGroup4D(lat, time_group, i)) {
      ((LatticeAMR *)lat)->GetGlobalVB(((LatticeAMR *)lat)->GetRank(i), 
				       from, to, &min_t, &max_t);
      pts[n++] = (double)((from[0] + to[0]) * 0.5);
      pts[n++] = (double)((from[1] + to[1]) * 0.5);
      pts[n++] = (double)((from[2] + to[2]) * 0.5);
    }
  }

  *err = 0;

}
//-----------------------------------------------------------------------
//
// end of zoltan callbacks
//
//-----------------------------------------------------------------------

#endif
