
#include "Partition.h"

//--------------------------------------------------------------------------
//
// Partition
//
// constructs and initializes the partitions
//
// npart: total (global) number of partitions (time and and space)
// nproc: total number of processes
//
Partition::Partition(int npart, int nproc) {

  int nn = ceil(npart / nproc); // max number of parts per proc
  int i, j;

  this->npart = npart;
  this->nproc = nproc;

  // allocate partitions list
  assert((parts = new Partition4D[npart]) != NULL);

  // init the partitions list
  for (j = 0; j < npart; j++)  {

    // allocate sending and receiving lists for one neighbor
    assert((parts[j].NumSendPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].SizeSendPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].SendPoints = (float **)malloc(sizeof(float *))) != NULL);
    assert((parts[j].NumRecvPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].SizeRecvPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].RecvPoints = (float **)malloc(sizeof(float *))) != NULL);

    // allocate sending and receiving points for the one neighbor
    assert((parts[j].SendPoints[0] = (float *)malloc(4 * sizeof(float)))
	   != NULL);
    assert((parts[j].RecvPoints[0] = (float *)malloc(4 * sizeof(float)))
	   != NULL);

    // number of points used and allocated
    parts[j].NumSendPoints[0] = parts[j].NumRecvPoints[0] = 0;
    parts[j].SizeSendPoints[0] = parts[j].SizeRecvPoints[0] = 
      4 * sizeof(float);

#ifdef MPI
    // request lists, with one neighbor initially
    assert((parts[j].SendReqs = (MPI_Request *)
	    malloc(sizeof(MPI_Request))) != NULL);
    assert((parts[j].RecvReqs = (MPI_Request *)
	    malloc(sizeof(MPI_Request))) != NULL);
#endif

    // default procs
    parts[j].Proc = j;

    // number of neighbors used and allocated
    parts[j].NumNeighbors = 0;
    parts[j].AllocNeighbors = 1;

    // other status info
#ifdef MPI
    parts[j].NumSendReqs = 0;
    parts[j].NumRecvReqs = 0;
#endif
    ClearReq(j);
    ClearLoad(j);
    ClearComp(j);

  }

  // allocate process indices
  assert((proc_parts = (int **)malloc(nproc * sizeof(int *))) != NULL);
  for (i = 0; i < nproc; i++)
    assert((proc_parts[i] = (int *)malloc(nn * sizeof(int))) != NULL);
  assert((proc_nparts = (int *)malloc(nproc * sizeof(int))) != NULL);

  assert((proc_neighbors = (int ***)malloc(nproc * sizeof(int **))) != NULL);
  for (i = 0; i < nproc; i++) {
    assert((proc_neighbors[i] = (int **)malloc(nn * sizeof(int *))) != NULL);
    for (j = 0; j < nn; j++)
      assert((proc_neighbors[i][j] = (int *)malloc(sizeof(int))) != NULL);
  }

  assert((proc_nneighbors = (int **)malloc(nproc * sizeof(int *))) != NULL);
  for (i = 0; i < nproc; i++)
    assert((proc_nneighbors[i] = (int *)malloc(nn * sizeof(int))) != NULL);
  for (i = 0; i < nproc; i++) {
    for (j = 0; j < nn; j++)
      proc_nneighbors[i][j] = 0;
  }

  assert((proc_aneighbors = (int **)malloc(nproc * sizeof(int *))) != NULL);
  for (i = 0; i < nproc; i++)
    assert((proc_aneighbors[i] = (int *)malloc(nn * sizeof(int))) != NULL);
  for (i = 0; i < nproc; i++) {
    for (j = 0; j < nn; j++)
      proc_aneighbors[i][j] = 1;
  }

}
//--------------------------------------------------------------------------
//
// ~Partition
// frees memory
//
Partition::~Partition()
{

  int i, j;

  for (i = 0; i < nproc; i++)
    free(proc_parts[i]);
  free(proc_parts);

  for (j = 0; j < npart; j++) {

    for (i = 0; i < parts[j].AllocNeighbors; i++) {
      free(parts[j].SendPoints[i]);
      free(parts[j].RecvPoints[i]);
    }

    free(parts[j].NumSendPoints);
    free(parts[j].SizeSendPoints);
    free(parts[j].SendPoints);
    free(parts[j].NumRecvPoints);
    free(parts[j].SizeRecvPoints);
    free(parts[j].RecvPoints);

  }

  if (parts != NULL)
    delete [] parts; 

}
//---------------------------------------------------------------------------
//
// GetProc
//
// gets the process number for a given rank
//
// myrank: my global partition number
//
int Partition::GetProc(int myrank) {

  return(parts[myrank].Proc); 

}
//--------------------------------------------------------------------------
//
// PostPoint
//
// posts a point for sending to a neighbor
// myrank: my global partition number
// p: 4D point
// neighbor: number of the neighbor
//
void Partition::PostPoint(int myrank, VECTOR4 p, int neighbor) {

  while (parts[myrank].SizeSendPoints[neighbor] < 
      (parts[myrank].NumSendPoints[neighbor] + 1) * 4 * sizeof(float)) {

    assert((parts[myrank].SendPoints[neighbor] = (float *)realloc(
	parts[myrank].SendPoints[neighbor],
	parts[myrank].SizeSendPoints[neighbor] * 2)) != NULL);

    parts[myrank].SizeSendPoints[neighbor] *= 2;

  }

  parts[myrank].SendPoints[neighbor][4 * 
       parts[myrank].NumSendPoints[neighbor] + 0] = p[0];
  parts[myrank].SendPoints[neighbor][4 * 
       parts[myrank].NumSendPoints[neighbor] + 1] = p[1];
  parts[myrank].SendPoints[neighbor][4 * 
       parts[myrank].NumSendPoints[neighbor] + 2] = p[2];
  parts[myrank].SendPoints[neighbor][4 * 
       parts[myrank].NumSendPoints[neighbor] + 3] = p[3];

  parts[myrank].NumSendPoints[neighbor]++;

}
//------------------------------------------------------------------------
//
// PrintPost
//
// prints the posted points
// myrank: global partition number
// 
void Partition::PrintPost(int myrank) {

  int i, j;

  fprintf(stderr, "\nPosted points list for rank %d\n", myrank);

  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    if (parts[myrank].NumSendPoints[i])
      fprintf(stderr, "rank %d posted %d points to neighbor %d\n", 
	      myrank, parts[myrank].NumSendPoints[i], i);

    if (parts[myrank].NumSendPoints[i]) {
      for (j = 0; j < parts[myrank].NumSendPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\t%.3f\n",
            parts[myrank].SendPoints[i][4 * j + 0],
	    parts[myrank].SendPoints[i][4 * j + 1], 
	    parts[myrank].SendPoints[i][4 * j + 2],
            parts[myrank].SendPoints[i][4 * j + 3]);
    }

  }

  fprintf(stderr,"\n");

}
//--------------------------------------------------------------------------
//
// PrintRecv
//
// prints the received points
// myrank: global partition number
//
void Partition::PrintRecv(int myrank) {

  int i, j;

  fprintf(stderr, "\nReceived points list for rank %d\n", myrank);

  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    if (parts[myrank].NumRecvPoints[i]) {
      fprintf(stderr, 
           "rank %d received %d points from neighbor %d:\n", 
            myrank, parts[myrank].NumRecvPoints[i], i);
      for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\t%.3f\n",
		parts[myrank].RecvPoints[i][4 * j + 0],
		parts[myrank].RecvPoints[i][4 * j + 1],
		parts[myrank].RecvPoints[i][4 * j + 2],
		parts[myrank].RecvPoints[i][4 * j + 3]);
    }

  }

  fprintf(stderr,"\n");

}
//---------------------------------------------------------------------------
//
// GetRecvPts
//
// copies received points from partition to user supplied list
// myrank: global partition number
//
// caller must ensure that list has enough room
//
//
void Partition::GetRecvPts(int myrank, VECTOR4 *ls) {

  int num = 0;
  int i, j;

  // copy points
  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++)
      (ls[num++]).Set(
		      parts[myrank].RecvPoints[i][4 * j + 0], 
		      parts[myrank].RecvPoints[i][4 * j + 1], 
		      parts[myrank].RecvPoints[i][4 * j + 2],
		      parts[myrank].RecvPoints[i][4 * j + 3]);

  }

  // clear receive lists
  for (i = 0; i < parts[myrank].NumNeighbors; i++)
    parts[myrank].NumRecvPoints[i] = 0;

}
//---------------------------------------------------------------------------
//
// adds a neighbor to global partition structure and
// process index into neighbors
// grows size of data structures if necessary
//
// myrank: my global partition number
// myblock: my local block number
// neighrank: neighbor global partition number
//
void Partition::AddNeighbor(int myrank, int myblock, int neighrank) {

  int p; // process
  int n; // number of neighbors

  // grow the parts data structure
  while (parts[myrank].AllocNeighbors < parts[myrank].NumNeighbors + 1)
    GrowNeighbors(myrank);

  // grow the proc_neighbors index
  p = parts[neighrank].Proc;
  n = proc_nneighbors[p][myblock]++;
  while (proc_aneighbors[p][myblock] < proc_nneighbors[p][myblock] + 1) {
    assert((proc_neighbors[p][myblock] = 
	    (int *)realloc(proc_neighbors[p][myblock], 
			   proc_aneighbors[p][myblock] * 2 * 
			   sizeof(int))) != NULL);
    proc_aneighbors[p][myblock] *= 2;
  }

  // update proc_neighbors and parts
  proc_neighbors[p][myblock][n] = parts[myrank].NumNeighbors++;

}
//---------------------------------------------------------------------------
//
// GrowNeighbors
//
// doubles the size of the neighbor-dependent structures
//
// myrank: my global partition number
//
void Partition::GrowNeighbors(int myrank) {

  int nn = parts[myrank].AllocNeighbors;  // number of neighbors allocated
  int i;

  // send points list for all neighbors
  assert((parts[myrank].NumSendPoints =
	  (int *)realloc(parts[myrank].NumSendPoints,
			 nn * sizeof(int) * 2)) != NULL);
  assert((parts[myrank].SizeSendPoints =
	  (int *)realloc(parts[myrank].SizeSendPoints, 
			 nn * sizeof(int) * 2)) != NULL);
  assert((parts[myrank].SendPoints =
	  (float **)realloc(parts[myrank].SendPoints, 
			    nn * sizeof(float *) * 2)) != NULL);

  // receive points list for all neighbors
  assert((parts[myrank].NumRecvPoints =
	  (int *)realloc(parts[myrank].NumRecvPoints,
			 nn * sizeof(int) * 2)) != NULL);
  assert((parts[myrank].SizeRecvPoints =
	  (int *)realloc(parts[myrank].SizeRecvPoints,
			 nn * sizeof(int) * 2)) != NULL);
  assert((parts[myrank].RecvPoints =
	  (float **)realloc(parts[myrank].RecvPoints,
			    nn * sizeof(float *) * 2)) != NULL);

  // sending and receiving points for individual neighbors
  for (i = nn; i < 2 * nn; i++) {
    assert((parts[myrank].SendPoints[i] = (float *)malloc(4 * sizeof(float)))
	   != NULL);
    assert((parts[myrank].RecvPoints[i] = (float *)malloc(4 * sizeof(float)))
	   != NULL);
    parts[myrank].NumSendPoints[i] = parts[myrank].NumRecvPoints[i] = 0;
    parts[myrank].SizeSendPoints[i] = parts[myrank].SizeRecvPoints[i] = 
      4 * sizeof(float);
  }

  // requests
#ifdef MPI
  assert((parts[myrank].SendReqs =
	  (MPI_Request *)realloc(parts[myrank].SendReqs,
				 nn * sizeof(MPI_Request) * 2)) != NULL);
  assert((parts[myrank].RecvReqs =
	  (MPI_Request *)realloc(parts[myrank].RecvReqs,
				 nn * sizeof(MPI_Request) * 2)) != NULL);
#endif

  parts[myrank].AllocNeighbors *= 2;

}
//---------------------------------------------------------------------------

// MPI versions of communication
// todo: make shared memory versions

#ifdef MPI

//---------------------------------------------------------------------------
//
// ExchangeNeighbors
//
// exhanges points with all neighbors
//
// neighbor_ranks: ranks (global partition numbers) of all my neighbors
// seeds: locations to store received points
// size_seeds: sizes of seed arrays (will be grown automatically if necessary)
// commm: MPI communicator
//
void Partition::ExchangeNeighbors(int **neighbor_ranks, VECTOR4 **seeds, int *size_seeds, MPI_Comm comm) {

  int nn = 0; // total number of neighbors
  int groupsize; // size of comm
  int *SendCounts, *RecvCounts; // count information
  int *CountSizes; // sizes index into count information
  int *CountDispls; // displacements index into count information
  int *SendPoints, *RecvPoints; // points
  int *SendPointDispls, *RecvPointDispls; // point displacements
  int *SendPointSizes, *RecvPointSizes; // point sizes
  int p; // process number
  int b; // block number
  int nps, npr; // number of sending and receiving points
  int rank; // global partition rank
  int myproc; // my MPI process
  int neigh_block; // (my) block number of my neighbor
  int i, j, k, n;

  MPI_Comm_size(comm, &groupsize);
  MPI_Comm_rank(comm, &myproc);

  // allocate count information arrays
  for (i = 0; i < proc_nparts[myproc]; i++)
    nn += parts[proc_parts[myproc][i]].NumNeighbors;
  assert((SendCounts = (int *)malloc(nn * 2 * sizeof(int))) != NULL);
  assert((RecvCounts = (int *)malloc(nn * 2 * sizeof(int))) != NULL);
  assert((CountSizes = (int *)malloc(groupsize * sizeof(int))) != NULL);
  assert ((CountDispls = (int *)malloc(groupsize * sizeof(int))) != NULL);

  // aggregate my counts into one message, get sizes and displacements
  n = 0;
  for (p = 0; p < groupsize; p++) { // all processes
    CountSizes[p] = 0;
    CountDispls[p] = n;
    for (i = 0; i < proc_nparts[myproc]; i++) { // my blocks
      rank = proc_parts[myproc][i];
      // blocks belonging to process p that are neighbors of block i
      for (j = 0; j < proc_nneighbors[p][i]; j++) {
	neigh_block = proc_neighbors[p][i][j];
	SendCounts[n++] = neighbor_ranks[i][neigh_block];
	SendCounts[n++] = parts[rank].NumSendPoints[neigh_block];
	CountSizes[p] += 2;
      }				  
    }
  }

  // exchange the count information
  // receive sizes and displacements same as send --
  // neighbor relation is symmetric
  MPI_Alltoallv(SendCounts, CountSizes, CountDispls, MPI_INT, 
		RecvCounts, CountSizes, CountDispls, MPI_INT, comm);

  // allocate point arrays
  nps = 0;
  npr = 0;
  for (i = 0; i < nn; i++) {
    nps += SendCounts[i * 2 + 1];
    npr += RecvCounts[i * 2 + 1];
  }
  assert((SendPoints = (int *)malloc(nps * 4 * sizeof(float))) != NULL);
  assert((RecvPoints = (int *)malloc(npr * 4 * sizeof(float))) != NULL);
  assert((SendPointSizes = (int *)malloc(groupsize * sizeof(int))) != NULL);
  assert((RecvPointSizes = (int *)malloc(groupsize * sizeof(int))) != NULL);
  assert ((SendPointDispls = (int *)malloc(groupsize * sizeof(int))) != NULL);
  assert ((RecvPointDispls = (int *)malloc(groupsize * sizeof(int))) != NULL);

  // get sizes and displacements for points exchange
  SendPointDispls[0] = 0;
  RecvPointDispls[0] = 0;
  for (p = 0; p < groupsize; p++) {

    // displacements
    if (p > 0) {
      SendPointDispls[p] = SendPointDispls[p - 1] + SendPointSizes[p - 1];
      RecvPointDispls[p] = RecvPointDispls[p - 1] + RecvPointSizes[p - 1];
    }

    // sizes (total number of points to send to a process)
    SendPointSizes[p] = 0;
    RecvPointSizes[p] = 0;
    for (i = 0; i < CountSizes[p] / 2; i++) {
      SendPointSizes[p] += 4 * SendCounts[CountDispls[p] + i * 2 + 1];
      RecvPointSizes[p] += 4 * RecvCounts[CountDispls[p] + i * 2 + 1];
    }

  }

  // aggregate my points into one message
  n = 0;
  for (p = 0; p < groupsize; p++) { // all processes
    for (i = 0; i < proc_nparts[myproc]; i++) { // my blocks
      rank = proc_parts[myproc][i];
      // blocks belonging to process p that are neighbors of block i
      for (j = 0; j < proc_nneighbors[p][i]; j++) {
	neigh_block = proc_neighbors[p][i][j];
	// points going to this neighbor
	for (k = 0; k < parts[rank].NumSendPoints[neigh_block]; k++) {
	  SendPoints[n++] =
	    parts[rank].SendPoints[neigh_block][4 * k + 0];
	  SendPoints[n++] =
	    parts[rank].SendPoints[neigh_block][4 * k + 1];
	  SendPoints[n++] =
	    parts[rank].SendPoints[neigh_block][4 * k + 2];
	  SendPoints[n++] =
	    parts[rank].SendPoints[neigh_block][4 * k + 3];
	}
      }				  
    }
  }

  // exchange the points
  MPI_Alltoallv(SendPoints, SendPointSizes, SendPointDispls, MPI_FLOAT, 
		RecvPoints, RecvPointSizes, RecvPointDispls, MPI_FLOAT, comm);

  // unpack the received points
  j = 0;
  for (n = 0; n < nn; n++) {
    // find my block number from the partition rank
    for (b = 0; b < proc_nparts[myproc]; b++) {
      if (proc_parts[myproc][b] == RecvCounts[n * 2 + 0])
	break;
    }
    assert(b < proc_nparts[myproc]);
    // grow size of seeds
    if (!size_seeds[b]) {
      assert((seeds[b] = (VECTOR4 *)
	      malloc(RecvCounts[n * 2 + 1] * sizeof(VECTOR4))) != NULL);
      size_seeds[b] = RecvCounts[n * 2 + 1] * sizeof(VECTOR4);
    }
    while (size_seeds[b] < RecvCounts[n * 2 + 1] * sizeof(VECTOR4)) {
      assert((seeds[b] = (VECTOR4 *)realloc(seeds[b], size_seeds[b] * 2))
	     != NULL);
      size_seeds[b] *= 2;
    }
    // copy points to seeds
    for (i = 0; i < RecvCounts[n * 2 + 1]; i++) {
      seeds[b][i].Set(RecvPoints[4 * j + 0], RecvPoints[4 * j + 1],
		      RecvPoints[4 * j + 2], RecvPoints[4 * j + 3]);
      j++;
    }
  }

}
//------------------------------------------------------------------------

#endif
