//
// THIS ENTIRE FILE IS DEPRECATED

//
//---------------------------------------------------------------------------
//
// Partition.C
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// All rights reserved. May not be used, modified, or copied
// without permission
//
//--------------------------------------------------------------------------

#include "Partition.h"

//--------------------------------------------------------------------------
//
// constructs and initializes the partitions
//
// npart: total (global) number of partitions (time and and space)
// nproc: total number of processes (pass 1 for serial code)
// myproc: my process id (pass 0 for serial code)
// track_ids: keep track of ids, defaults to false
//
Partition::Partition(int npart, int nproc, int myproc, bool track_ids) {

  int j;

  this->npart = npart;
  this->nproc = nproc;
  this->myproc = myproc;
  this->track_ids = track_ids;

  // allocate partitions list
  assert((parts = new Partition4D[npart]) != NULL);

  // init the partitions list
  nb = 0;
  for (j = 0; j < npart; j++)  {

    // sending lists
    assert((parts[j].NumSendPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].SizeSendPoints = (int *)malloc(sizeof(int))) != NULL);
    assert((parts[j].SendPoints = (float **)malloc(sizeof(float *))) != NULL);
    if (this->track_ids)
      assert((parts[j].SendPointIds = 
	      (int64_t **)malloc(sizeof(int64_t *))) != NULL);
    
    parts[j].NumSendPoints[0] = 0;
    parts[j].SizeSendPoints[0] = 0;

    // number of neighbors used and allocated
    parts[j].NumNeighbors = 0;
    parts[j].AllocNeighbors = 1;

    // other status info
    ClearLoad(j);

  }

  tag = 0;

}
//--------------------------------------------------------------------------
//
// destructor
//
Partition::~Partition()
{

  int i, j;

  for (j = 0; j < npart; j++) {

    for (i = 0; i < parts[j].NumNeighbors; i++) {
      if (parts[j].SizeSendPoints[i]) {
	free(parts[j].SendPoints[i]);
	if (this->track_ids)
		free(parts[j].SendPointIds[i]);
			}
    }

    free(parts[j].NumSendPoints);
    free(parts[j].SizeSendPoints);
    free(parts[j].SendPoints);
		if (this->track_ids)
			free(parts[j].SendPointIds);

  }

  if (parts != NULL)
    delete [] parts; 

}
//---------------------------------------------------------------------------
//
// posts a point for sending to a neighbor
//
// gid: my global partition number
// p: 4D point
// neighbor: number of the neighbor
//
void Partition::PostPoint(int gid, float *p, int neighbor, 
			  int64_t seed_id) {

  // seed id must be given if id tracking is turned on
  assert((this->track_ids == true && seed_id != -1) 
	 || (this->track_ids == false && seed_id == -1));

  unsigned int size = parts[gid].SizeSendPoints[neighbor];

  while (size < 
	 (parts[gid].NumSendPoints[neighbor] + 1) * 4 * sizeof(float)) {

    if (size == 0) {
      size = 4 * sizeof(float);
      assert((parts[gid].SendPoints[neighbor] = 
	      (float *)malloc(size)) != NULL);
      // assume int64_t = point / 2
      if (this->track_ids)
	assert((parts[gid].SendPointIds[neighbor] = 
	      	(int64_t *)malloc(size / 2)) != NULL);
    }
    else {
      size *= 2;
      assert((parts[gid].SendPoints[neighbor] = 
	      (float *)realloc(parts[gid].SendPoints[neighbor], 
			       size)) != NULL);
      // assume int64_t = point / 2
      if (this->track_ids)
	assert((parts[gid].SendPointIds[neighbor] = 
		(int64_t *)realloc(parts[gid].SendPointIds[neighbor], 
				   size / 2)) != NULL);
    }

  }

  parts[gid].SizeSendPoints[neighbor] = size;

  parts[gid].SendPoints[neighbor]
    [4 * parts[gid].NumSendPoints[neighbor] + 0] = p[0];
  parts[gid].SendPoints[neighbor]
    [4 * parts[gid].NumSendPoints[neighbor] + 1] = p[1];
  parts[gid].SendPoints[neighbor]
    [4 * parts[gid].NumSendPoints[neighbor] + 2] = p[2];
  parts[gid].SendPoints[neighbor]
    [4 * parts[gid].NumSendPoints[neighbor] + 3] = p[3];

  if (this->track_ids)
    parts[gid].SendPointIds[neighbor]
      [parts[gid].NumSendPoints[neighbor]] = seed_id;

  parts[gid].NumSendPoints[neighbor]++;

}
//------------------------------------------------------------------------
// //
// // prints the posted points
// // gid: global partition number
// // 
// void Partition::PrintPost(int gid) {

//   int i, j;


//   for (i = 0; i < parts[gid].NumNeighbors; i++) {

//     if (parts[gid].NumSendPoints[i])
//       fprintf(stderr, "gid %d posted %d points to neighbor %d\n", 
// 	      gid, parts[gid].NumSendPoints[i], i);

//     if (parts[gid].NumSendPoints[i]) {
//       for (j = 0; j < parts[gid].NumSendPoints[i]; j++)
// 	fprintf(stderr, "%.3f\t%.3f\t%.3f\t%.3f\n",
//             parts[gid].SendPoints[i][4 * j + 0],
// 	    parts[gid].SendPoints[i][4 * j + 1], 
// 	    parts[gid].SendPoints[i][4 * j + 2],
//             parts[gid].SendPoints[i][4 * j + 3]);
//     }

//   }

//   fprintf(stderr,"\n");

// }
//--------------------------------------------------------------------------
// //
// // removes a block from the partition structure and proc_neighbors
// //
// // gid: my global partition number
// //
// void Partition::RemoveBlock(int gid) {

//   int i;

//   parts[gid].Proc = -1;
//   parts[gid].NumNeighbors = 0;
//   for (i = 0; i < parts[gid].AllocNeighbors; i++)
//     parts[gid].NumSendPoints[i] = 0;
//   nb--;

// }
// //---------------------------------------------------------------------------
// //
// // adds a block to the partition structure
// //
// // gid: my global partition number
// //
// void Partition::AddBlock(int gid) {

//   int i;

//   parts[gid].Proc = myproc;
//   parts[gid].NumNeighbors = 0;
//   for (i = 0; i < parts[gid].AllocNeighbors; i++)
//     parts[gid].NumSendPoints[i] = 0;
//   nb++;

// }
//---------------------------------------------------------------------------
//
// adding a neighbor involves growing size of data structures if necessary
//
// gid: my global partition number
//
void Partition::AddNeighbor(int gid) {

  int n; // new number of neighbors allocated
  int i;

  // grow the parts data structure
  while (parts[gid].AllocNeighbors < parts[gid].NumNeighbors + 1) {

    if (parts[gid].AllocNeighbors == 0) {

      n = 1;

      // send list for all neighbors
      assert((parts[gid].NumSendPoints =
	      (int *)malloc(n * sizeof(int))) != NULL);
      assert((parts[gid].SizeSendPoints =
	      (int *)malloc(n * sizeof(int))) != NULL);
      assert((parts[gid].SendPoints =
	      (float **)malloc(n * sizeof(float *))) != NULL);
			if (this->track_ids)
				assert((parts[gid].SendPointIds =
	 	     (int64_t **)malloc(n * sizeof(int64_t *))) != NULL);

    }
    else {

      n = parts[gid].AllocNeighbors * 2;

      // send list for all neighbors
      assert((parts[gid].NumSendPoints =
	      (int *)realloc(parts[gid].NumSendPoints, 
			     n * sizeof(int))) != NULL);
      assert((parts[gid].SizeSendPoints =
	      (int *)realloc(parts[gid].SizeSendPoints, 
			     n * sizeof(int))) != NULL);
      assert((parts[gid].SendPoints =
	      (float **)realloc(parts[gid].SendPoints, 
				n * sizeof(float *))) != NULL);
			if (this->track_ids)
				assert((parts[gid].SendPointIds =
		      (int64_t **)realloc(parts[gid].SendPointIds, 
					n * sizeof(int64_t *))) != NULL);

    }

    // sending and receiving points for individual neighbors
    for (i = parts[gid].AllocNeighbors; i < n; i++) {
      parts[gid].NumSendPoints[i] = 0;
      parts[gid].SizeSendPoints[i] = 0;
    }

    parts[gid].AllocNeighbors = n;

  }

  parts[gid].NumNeighbors++;

}
//---------------------------------------------------------------------------

// // MPI version of communication

// #ifdef _MPI

// //---------------------------------------------------------------------------
// //
// // exhanges points with all neighbors
// //
// // all-all-version
// //
// // DEPRECATED; remove when no longer needed
// //
// // block_gids: global ids of those blocks
// // neighbor_gids: gids (global partition numbers) of all my neighbors
// // seeds: locations to store received points, indexed by local block number
// // alloc_seeds: allocated number of seeds for each block 
// // (grown automatically if necessary)
// // num_seeds: actual number of seeds stored for each block
// // commm: MPI communicator
// //
// // returns: total number of points received by this process
// //
// int Partition::OldExchangeNeighbors(int *block_gids, int **neighbor_gids, 
// 				 int **neighbor_procs, float ***seeds, 
// 				 int *alloc_seeds, int *num_seeds, 
// 				 MPI_Comm comm, int64_t **seed_ids) {
	
//   int nn = 0; // total number of neighbors
//   int *SendCounts = NULL; // send count information
//   int *RecvCounts = NULL; // receive count information
//   int *CountSizes; // sizes index into count information
//   int *CountDispls; // displacements index into count information
//   float *SendPoints = NULL; // send points
//   float *RecvPoints = NULL; // receive points
//   int64_t *SendPointIds = NULL; // send point ids
//   int64_t *RecvPointIds = NULL; // receive point ids
//   int *SendPointDispls, *RecvPointDispls; // point displacements
//   int *SendPointSizes, *RecvPointSizes; // point sizes
//   int p; // process number
//   int b; // block number
//   int nps, npr; // number of sending and receiving points
//   int r; // global partition gid
//   int np; // number of points from a neighbor
//   int ofst; // offset into received points
//   int i, j, k, n, m;

//   // allocate count information arrays
//   for (i = 0; i < nb; i++)
//     nn += parts[block_gids[i]].NumNeighbors;
//   if (nn > 0) {
//     assert((SendCounts = (int *)malloc(nn * 2 * sizeof(int))) != NULL);
//     assert((RecvCounts = (int *)malloc(nn * 2 * sizeof(int))) != NULL);
//   }
//   assert((CountSizes = (int *)malloc(nproc * sizeof(int))) != NULL);
//   assert((CountDispls = (int *)malloc(nproc * sizeof(int))) != NULL);

//   // aggregate my counts into one message, get sizes and displacements
//   n = 0;
//   for (p = 0; p < nproc; p++) { // all processes
//     CountSizes[p] = 0;
//     CountDispls[p] = n;
//     for (i = 0; i < nb; i++) { // my blocks
//       r = block_gids[i]; // gid of my (sending) block
//       for (j = 0; j < parts[r].NumNeighbors; j++) { // neighbors of block i
// 	if (neighbor_procs[i][j] == p) {
// 	  SendCounts[n++] = neighbor_gids[i][j];
// 	  SendCounts[n++] = parts[r].NumSendPoints[j];
// 	  CountSizes[p] += 2;
// 	}
//       }				  
//     }
//   }

//   // exchange the count information
//   // receive sizes and displacements same as send --
//   // neighbor relation is symmetric
//   MPI_Alltoallv(SendCounts, CountSizes, CountDispls, MPI_INT, 
// 		RecvCounts, CountSizes, CountDispls, MPI_INT, comm);

//   // allocate point arrays
//   nps = 0;
//   npr = 0;
//   for (p = 0; p < nproc; p++) {
//     for (i = 0; i < CountSizes[p] / 2; i++) {
//       nps += SendCounts[CountDispls[p] + i * 2 + 1];
//       npr += RecvCounts[CountDispls[p] + i * 2 + 1];
//     }
//   }
//   if (nps > 0) {
//     assert((SendPoints = (float *)malloc(nps * 4 * sizeof(float))) != NULL);
//     if (this->track_ids)
//       assert((SendPointIds = (int64_t *)malloc(nps * sizeof(int64_t))) != NULL);
//   }
//   if (npr > 0) {
//     assert((RecvPoints = (float *)malloc(npr * 4 * sizeof(float))) != NULL);
//     if (this->track_ids)
//       assert((RecvPointIds = (int64_t *)malloc(npr * sizeof(int64_t))) != NULL);
//   }
//   assert((SendPointSizes  = (int *)malloc(nproc * sizeof(int))) != NULL);
//   assert((RecvPointSizes  = (int *)malloc(nproc * sizeof(int))) != NULL);
//   assert((SendPointDispls = (int *)malloc(nproc * sizeof(int))) != NULL);
//   assert((RecvPointDispls = (int *)malloc(nproc * sizeof(int))) != NULL);

//   // pack my points for sending
//   SendPointDispls[0] = RecvPointDispls[0] = 0;
//   for (p = 0; p < nproc; p++) { // all processes

//     // Displacements
//     if (p > 0) {
//       SendPointDispls[p] = SendPointDispls[p - 1] + SendPointSizes[p - 1];
//       RecvPointDispls[p] = RecvPointDispls[p - 1] + RecvPointSizes[p - 1];
//     }

//     // sizes
//     SendPointSizes[p] = 0;
//     RecvPointSizes[p] = 0;
//     for (i = 0; i < CountSizes[p] / 2; i++) {
//       SendPointSizes[p] += 4 * SendCounts[CountDispls[p] + i * 2 + 1];
//       RecvPointSizes[p] += 4 * RecvCounts[CountDispls[p] + i * 2 + 1];
//     }

//     np = 0; // total number of points I have packed
//     for (i = 0; i < nb; i++) { // my blocks

//       r = block_gids[i]; // gid of my (sending) block
//       for (j = 0; j < parts[r].NumNeighbors; j++) { // neighbors of block i
// 	if (neighbor_procs[i][j] == p) {
// 	  ofst = SendPointDispls[p] + 4 * np;
// 	  for (k = 0; k < parts[r].NumSendPoints[j]; k++) {
// 	    SendPoints[ofst + 4 * k + 0] = parts[r].SendPoints[j][4 * k + 0];
// 	    SendPoints[ofst + 4 * k + 1] = parts[r].SendPoints[j][4 * k + 1];
// 	    SendPoints[ofst + 4 * k + 2] = parts[r].SendPoints[j][4 * k + 2];
// 	    SendPoints[ofst + 4 * k + 3] = parts[r].SendPoints[j][4 * k + 3];
// 	    if (this->track_ids)
// 	      // ofst / 4 because ofst = pts = 4 floats
// 	      SendPointIds[(ofst / 4) + k] = parts[r].SendPointIds[j][k];
// 	  }
// 	  np += parts[r].NumSendPoints[j];
// 	  parts[r].NumSendPoints[j] = 0;
// 	}

//       } // neighbors of block i				  

//     } // myblocks

//   } // all procs

//   // exchange the points
//   MPI_Alltoallv(SendPoints, SendPointSizes, SendPointDispls, MPI_FLOAT, 
// 		RecvPoints, RecvPointSizes, RecvPointDispls, MPI_FLOAT, comm);

//   if (this->track_ids) {
//     // the send point sizes and recv point sizes are in terms of points
//     // divide them by to get them in terms of seed ids
//     int *SendPointIdSizes = (int *)malloc(sizeof(int) * nproc);	
//     int *SendPointIdDispls = (int *)malloc(sizeof(int) * nproc);	
//     int *RecvPointIdSizes = (int *)malloc(sizeof(int) * nproc);	
//     int *RecvPointIdDispls = (int *)malloc(sizeof(int) * nproc);	
//     assert(SendPointIdSizes != NULL && SendPointIdDispls != NULL && 
// 	   RecvPointIdSizes != NULL && RecvPointIdDispls != NULL);
//     for (i = 0; i < nproc; i++) {
//       SendPointIdSizes[i] = SendPointSizes[i] / 4;
//       SendPointIdDispls[i] = SendPointDispls[i] / 4;
//       RecvPointIdSizes[i] = RecvPointSizes[i] / 4;
//       RecvPointIdDispls[i] = RecvPointDispls[i] / 4;
//     }
//     // exchange seed ids
//     MPI_Alltoallv(
// 		  SendPointIds, SendPointIdSizes, SendPointIdDispls, MPI_LONG_LONG, 
// 		  RecvPointIds, RecvPointIdSizes, RecvPointIdDispls, MPI_LONG_LONG, comm);
//     free(SendPointIdSizes);
//     free(SendPointIdDispls);
//     free(RecvPointIdSizes);
//     free(RecvPointIdDispls);
//   }	

//   // unpack the received points
//   for (p = 0; p < nproc; p++) { // all procs

//     np = 0; // total number of points I have unpacked
//     for (n = 0; n < CountSizes[p] / 2; n++) { // neighbors in each proc

//       // find my block number from the partition gid
//       for (b = 0; b < nb; b++) {
// 	if (block_gids[b] == RecvCounts[CountDispls[p] + n * 2])
// 	  break;
//       }
//       assert(b < nb); // sanity
//       m = RecvCounts[CountDispls[p] + n * 2 + 1]; // number of points recv'd
//                                                   // from this neighbor
//       // grow size of seeds
//       if (!alloc_seeds[b] && num_seeds[b] + m > 0) {
// 	assert((seeds[b] = (float **)
// 		malloc((num_seeds[b] + m) * sizeof(float *))) != NULL);
// 	for (i = 0; i < num_seeds[b] + m; i++)
// 	  assert((seeds[b][i] = (float *)
// 		  malloc(4 * sizeof(float))) != NULL);
// 	if (this->track_ids)
// 	  assert((seed_ids[b] = (int64_t *)
// 		  malloc((num_seeds[b] + m) * sizeof(int64_t))) != NULL);
// 	alloc_seeds[b] = num_seeds[b] + m;
//       }

//       while (alloc_seeds[b] < num_seeds[b] + m) {
// 	assert((seeds[b] = 
// 		(float **)realloc(seeds[b], 
// 				  alloc_seeds[b] * 2 * 
// 				  sizeof(float *))) != NULL);
// 	for (i = 0; i < alloc_seeds[b]; i++)
// 	  assert((seeds[b][alloc_seeds[b] + i] = (float *)
// 		  malloc(4 * sizeof(float))) != NULL);
// 	if (this->track_ids)
// 	  assert((seed_ids[b] = 
// 		  (int64_t *)realloc(seed_ids[b],
// 				     alloc_seeds[b] * 2 * 
// 				     sizeof(int64_t *))) != NULL);
// 	alloc_seeds[b] *= 2;
//       }

//       // copy points to seeds
//       for (i = 0; i < m; i++) {
// 	ofst = RecvPointDispls[p] + 4 * np;
// 	seeds[b][num_seeds[b]][0] = RecvPoints[ofst + 4 * i];
// 	seeds[b][num_seeds[b]][1] = RecvPoints[ofst + 4 * i + 1];
// 	seeds[b][num_seeds[b]][2] = RecvPoints[ofst + 4 * i + 2];
// 	seeds[b][num_seeds[b]][3] = RecvPoints[ofst + 4 * i + 3];
// 	if (this->track_ids)
// 	  seed_ids[b][num_seeds[b]] = RecvPointIds[ofst / 4 + i];
// 	num_seeds[b]++;
//       }

//       np += m;

//     }

//   }

//   // cleanup
//   if (nn > 0) {
//     free(SendCounts);
//     free(RecvCounts);
//   }
//   free(CountSizes);
//   free(CountDispls);
//   if(nps > 0) {
//     free(SendPoints);
//     if (this->track_ids)
//       free(SendPointIds);
//   }
//   if(npr > 0) {
//     free(RecvPoints);
//     if (this->track_ids)
//       free(RecvPointIds);
//   }
//   free(SendPointSizes);
//   free(RecvPointSizes);
//   free(SendPointDispls);
//   free(RecvPointDispls);

//   return npr;

// }
// //------------------------------------------------------------------------
// //
// // exchanges points with all neighbors
// //
// // point-to-point synchronous version (current)
// //
// // block_gids: global ids of my blocks
// // neighbor_gids: global ids of all my neighbor blocks
// // neighbor_procs: process ids of all my neighbor blocks
// // nproc: number of actual distinct processes in my neighborhood (output)
// // pts: points received in each block, from each process (output)
// // cts: counts of points received in each block, from each process (ouput)
// // pids: point ids received in each block, from each process (output)
// // commm: MPI communicator
// //
// // returns: total number of points received by this process
// //
// int Partition::SyncExchangeNeighbors(int *block_gids, int **neighbor_gids, 
// 					int **neighbor_procs, int *nproc, 
// 					float ***pts, int ***cts, 
// 					int64_t ***pids, MPI_Comm comm) {

//   int *nps, *npr; // numbers of sending and receiving points
//   int tot_npr = 0; // total number of received points
//   int *proc; // process id of all my neighbors
//   int nn; // my total number of neighbor blocks
//   int i, j, k, p, r;

//   // total number of neighbor blocks
//   nn = 0;
//   for (i = 0; i < nb; i++)
//     nn += parts[block_gids[i]].NumNeighbors;

//   if (this->nproc < nn)
//     assert((proc = (int *)malloc(this->nproc * sizeof(int))) != NULL);
//   else
//     assert((proc = (int *)malloc(nn * sizeof(int))) != NULL);

//   // count the number of procs with whom I need to talk
//   *nproc = 0;
//   for (i = 0; i < nb; i++) { // all my blocks
//     r = block_gids[i]; // gid of my (sending) block
//     for (j = 0; j < parts[r].NumNeighbors; j++) { // neighbors of block i
//       p = neighbor_procs[i][j]; // process of this neighbor
//       for (k = 0; k < *nproc; k++) { 
// 	if (proc[k] == p)
// 	  break;
//       }
//       if (k == *nproc) {
// 	proc[k] = p;
// 	(*nproc)++;
//       }
//     }
//   }

//   // allocate sending and receiving count messages
//   assert((nps = (int *)malloc(*nproc * sizeof(int))) != NULL);
//   assert((npr = (int *)malloc(*nproc * sizeof(int))) != NULL);

//   // exchange counts
//   SyncExchangeNeighborCounts(block_gids, neighbor_gids, neighbor_procs, nn, 
// 			 nps, npr, cts, proc, nproc, comm);

//   // exchange points
//   SyncExchangeNeighborPoints(block_gids, neighbor_gids, neighbor_procs, proc,
// 			 nproc, nps, npr, pts, pids, comm);

//   // total number of points recvd
//   for (i = 0; i < *nproc; i++)
//     tot_npr += npr[i];

// //   // debug
// //   fprintf(stderr, "pts recvd * %d from proc 0 %.3f %.3f %.3f * %d from proc 1 %.3f %.3f %.3f\n", (*cts)[0][0], (*pts)[0][0], (*pts)[0][1], (*pts)[0][2], (*cts)[1][0], (*pts)[1][0], (*pts)[1][1], (*pts)[1][2]);

//   free(proc);
//   free(nps);
//   free(npr);

//   return tot_npr;

// }
// //---------------------------------------------------------------------------
// //
// // exchanges counts with all neighbors
// //
// // block_gids: global ids of those blocks
// // neighbor_gids: gids (global partition numbers) of all my neighbors
// // neighbor_procs: process ids of all my neighbors
// // nn: the total number of my neighboring blocks
// // nps: total number of points sent to each process (ouput)
// // npr: total number of points received from each process (ouput)
// // recv_cts: number of points received in each block, from each process (ouput)
// // proc: the process ids with whom I communicate
// // nproc: the number of processes with whom I communicate
// // commm: MPI communicator
// //
// void Partition::SyncExchangeNeighborCounts(int *block_gids, 
// 					   int **neighbor_gids, 
// 					   int **neighbor_procs, int nn, 
// 					   int *nps, int *npr, int ***recv_cts, 
// 					   int *proc, int *nproc, 
// 					   MPI_Comm comm) {

//   int **send_cts; // sending messages for counts
//   int *mlen; // length of sending messages
//   MPI_Request *reqs; // array of requests, one per process communicating
//   MPI_Status *stats; // array of statuses, oner per process communicating
//   static int first = 1; // first time flag
//   static int old_nproc; // previous value of nproc
//   static int old_nn; // previous value of nn
//   int r; // block gid
//   int p; // process number
//   int i, j, k;

//   // allocate sending and receiving count messages
//   assert((mlen = new int[*nproc]) != NULL);
//   assert((send_cts = new int*[*nproc]) != NULL);
//   assert((reqs = new MPI_Request[2 * *nproc]) != NULL); // both send and recv
//   assert((stats = new MPI_Status[2 * *nproc]) != NULL);
//   for (i = 0; i < *nproc; i++)
//     assert((send_cts[i] = new int[nn * 2 + 1]) != NULL);
//   if (first) {
//     assert((*recv_cts = (int**)malloc(*nproc * sizeof(int *))) != NULL);
//     for (i = 0; i < *nproc; i++)
//       assert(((*recv_cts)[i] = (int *)malloc((nn * 2 + 1) * sizeof(int))) 
// 	     != NULL);
//     old_nproc = *nproc;
//     old_nn = nn;
//     first = 0;
//   }
//   if (*nproc > old_nproc) {
//     assert((*recv_cts = (int **)realloc(*recv_cts, *nproc * sizeof(int *))) 
// 	   != NULL);
//     for (i = old_nproc; i < *nproc; i++)
//       assert(((*recv_cts)[i] = (int *)malloc((nn * 2 + 1) * 
// 					   sizeof(int))) != NULL);
//     old_nproc = *nproc;
//   }
//   if (nn > old_nn) {
//     for (i = 0; i < old_nproc; i++) 
//       assert(((*recv_cts)[i] = (int *)realloc((*recv_cts)[i], (nn * 2 + 1) * 
// 					    sizeof(int))) != NULL);
//     old_nn = nn;
//   }

//   // init
//   for (i = 0; i < *nproc; i++) {
//     memset(send_cts[i], 0, (nn * 2 + 1) * sizeof(int));
//     mlen[i] = 1;
//     nps[i] = 0;
//     npr[i] = 0;
//   }

//   // pack count messages
//   for (i = 0; i < nb; i++) { // all my blocks
//     r = block_gids[i]; // gid of my (sending) block
//     for (j = 0; j < parts[r].NumNeighbors; j++) { // neighbors of block i
//       p = neighbor_procs[i][j]; // process of this neighbor
//       for (k = 0; k < *nproc; k++) {      
// 	if (proc[k] == p)
// 	  break;
//       }

//       // add to the message
//       send_cts[k][mlen[k]++] = neighbor_gids[i][j]; // part. gid of dest block
//       send_cts[k][mlen[k]++] = parts[r].NumSendPoints[j]; // number of pts
//       send_cts[k][0]++;
//     }
//   }


//   // exchange counts
//   for (p = 0; p < *nproc; p++)
//     MPI_Isend(send_cts[p], mlen[p], MPI_INT, proc[p], 0, comm, 
// 	      &reqs[2 * p]);
//   for (p = 0; p < *nproc; p++)
//     MPI_Irecv((*recv_cts)[p], nn * 2 + 1, MPI_INT, proc[p], 0, comm, 
// 	      &reqs[2 * p + 1]);
//   MPI_Waitall(*nproc * 2, reqs, stats);

//   // number of points to receive from each process
//   for (p = 0; p < *nproc; p++) {
//     for (i = 0; i < send_cts[p][0]; i++)
//       nps[p] += send_cts[p][i * 2 + 2];
//     for (i = 0; i < (*recv_cts)[p][0]; i++)
//       npr[p] += (*recv_cts)[p][i * 2 + 2];
//   }

//   for (i = 0; i < *nproc; i++)
//     delete[] send_cts[i];
//   delete[] send_cts;
//   delete[] mlen;
//   delete[] reqs;
//   delete[] stats;

// }
// //---------------------------------------------------------------------------
// //
// // exchanges points with all neighbors
// //
// // block_gids: global ids of those blocks
// // neighbor_gids: gids (global partition numbers) of all my neighbors
// // neighbor_procs: process ids of all my neighbors
// // proc: the process ids in my neighborhood
// // nproc: the number of processes in my neighborhood
// // nps: number of points sent in each block, to each process
// // npr: number of points received in each block, from each process
// // recv_pts: points received in each block, from each process (output)
// // recv_pids: point ids received in each block, from each process (output)
// // commm: MPI communicator
// //
// void Partition::SyncExchangeNeighborPoints(int *block_gids, 
// 					   int **neighbor_gids, 
// 					   int **neighbor_procs, int *proc, 
// 					   int *nproc, int *nps, int *npr, 
// 					   float ***recv_pts, 
// 					   int64_t ***recv_pids,
// 					   MPI_Comm comm) {

//   float **send_pts; // sending messages for points
//   int64_t **send_pids; // sending messages for point ids
//   int *mlen; // length of sending messages
//   MPI_Request *reqs; // array of requests, one per process communicating
//   MPI_Status *stats; // array of statuses, oner per process communicating
//   static int first = 1; // first time flag
//   static int *alloc_recv; // allocated size of send, recv messages
//   static int old_nproc; // previous value of nproc
//   int r; // block gid
//   int p; // process number
//   int i, j, k, n;

//   // allocate receiving point messages
//   if (first) {
//     assert((*recv_pts = (float **)malloc(*nproc * sizeof(float *)))
// 	   != NULL);
//     assert((alloc_recv = (int *)malloc(*nproc * sizeof(int)))
// 	   != NULL);
//     if (track_ids)
//       assert((*recv_pids = (int64_t **)malloc(*nproc * sizeof(int64_t *))) 
// 	     != NULL);
//     for (p = 0; p < *nproc; p++)
//       alloc_recv[p] = 0;
//     old_nproc = *nproc;
//     first = 0;
//   }
//   if (*nproc > old_nproc) {
//     assert((alloc_recv = (int *)realloc(alloc_recv, *nproc * sizeof(int)))
// 	   != NULL);
//     assert((*recv_pts = (float **)realloc(*recv_pts, *nproc * sizeof(float *)))
// 	   != NULL);
//     if (track_ids)
//       assert((*recv_pids = (int64_t **)realloc(*recv_pids, *nproc * 
// 					       sizeof(int64_t *))) != NULL);
//     for (p = old_nproc; p < *nproc; p++)
//       alloc_recv[p] = 0;
//     old_nproc = *nproc;
//   }

//   // allocate sending point messages
//   assert((mlen = new int[*nproc]) != NULL);
//   assert((send_pts = new float*[*nproc]) != NULL);
//   assert((send_pids = new int64_t*[*nproc]) != NULL);
//   assert((reqs = new MPI_Request[2 * *nproc]) != NULL); // both send and recv
//   assert((stats = new MPI_Status[2 * *nproc]) != NULL);
//   for (p = 0; p < *nproc; p++) {
//     assert((send_pts[p] = (float *)malloc(nps[p] * 
// 					  4 * sizeof(float))) != NULL);
//     if (track_ids)
//       assert((send_pids[p] = (int64_t *)malloc(nps[p] *
// 					       sizeof(int64_t))) != NULL);
//   }
//   for (i = 0; i < *nproc; i++)
//     mlen[i] = 0;

//   // grow size of received points
//   for (p = 0; p < *nproc; p++) {
//     if (!alloc_recv[p] && npr[p]) {
//       assert(((*recv_pts)[p] = (float *)malloc(npr[p] * 
// 					    4 * sizeof(float))) != NULL);
//       if (track_ids)
// 	assert(((*recv_pids)[p] = (int64_t *)malloc(npr[p] *
// 						 sizeof(int64_t))) != NULL);
//       alloc_recv[p] = npr[p];
//     }
//     while (alloc_recv[p] < npr[p]) {
//       assert(((*recv_pts)[p] = (float *)
// 	      realloc((*recv_pts)[p], alloc_recv[p] * 
// 		      2 * 4 * sizeof(float))) != NULL);
//       if (track_ids)
// 	assert(((*recv_pids)[p] = (int64_t *)
// 		realloc((*recv_pids)[p], alloc_recv[p] *
// 			2 * sizeof(int64_t))) != NULL);
//       alloc_recv[p] *= 2;
//     }

//   }

//   // pack messages containing points
//   for (i = 0; i < nb; i++) { // all my blocks
//     r = block_gids[i]; // gid of my (sending) block
//     for (j = 0; j < parts[r].NumNeighbors; j++) { // neighbors of block i
//       p = neighbor_procs[i][j]; // process this goes to
//       for (k = 0; k < *nproc; k++) {
// 	if (proc[k] == p)
// 	  break;
//       }
//       proc[k] = p;
//       // add to the message
//       for (n = 0; n < parts[r].NumSendPoints[j]; n++) { // the points
// 	send_pts[k][mlen[k]++] = parts[r].SendPoints[j][4 * n + 0];
// 	send_pts[k][mlen[k]++] = parts[r].SendPoints[j][4 * n + 1];
// 	send_pts[k][mlen[k]++] = parts[r].SendPoints[j][4 * n + 2];
// 	send_pts[k][mlen[k]++] = parts[r].SendPoints[j][4 * n + 3];
// 	if (track_ids)
// 	  send_pids[k][mlen[k]++] = parts[r].SendPointIds[j][n];
//       }
//       parts[r].NumSendPoints[j] = 0;
//     } // neighbors of block i
//   } // all my blocks

//   // exchange points
//   for (p = 0; p < *nproc; p++)
//     MPI_Isend(send_pts[p], nps[p] * 4, MPI_FLOAT, proc[p], 0, comm, 
// 	      &reqs[2 * p]);
//   for (p = 0; p < *nproc; p++)
//     MPI_Irecv((*recv_pts)[p], npr[p] * 4, MPI_FLOAT, proc[p], 0, comm, 
// 	      &reqs[2 * p + 1]);
//   MPI_Waitall(*nproc * 2, reqs, stats);

//   // exchange point ids
//   if (track_ids) {
//     for (p = 0; p < *nproc; p++)
//       MPI_Isend(send_pids[p], nps[p], MPI_LONG_LONG, proc[p], 0, comm, 
// 		&reqs[2 * p]);
//     for (p = 0; p < *nproc; p++)
//       MPI_Irecv((*recv_pids)[p], npr[p], MPI_LONG_LONG, proc[p], 0, comm, 
// 		&reqs[2 * p + 1]);
//     MPI_Waitall(*nproc * 2, reqs, stats);
//   }

//   // clean up
//   for (p = 0; p < *nproc; p++) {
//     free(send_pts[p]);
//     if (track_ids)
//       free(send_pids[p]);
//   }
//   delete[] mlen;
//   delete[] send_pts;
//   delete[] send_pids;
//   delete[] reqs;
//   delete[] stats;

// }
// //---------------------------------------------------------------------------

// #endif

