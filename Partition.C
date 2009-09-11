
#include "Partition.h"

//--------------------------------------------------------------------------
//
// Partition
//
// constructs and initializes the partitions
//
// nsp: total (global) number of spatial partitions
// ntp: total (global) number of temporal partitions
//
Partition::Partition(int nsp, int ntp) {

  int i, j;

  assert((npart  = nsp * ntp) <= MAX_PARTS);
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
    // 1-d request list, with one neighbor intially and 4 requests per neighbor
    assert((parts[j].Reqs = (MPI_Request *)
	    malloc(4 * sizeof(MPI_Request))) != NULL);
#endif

    // default procs
    parts[j].Proc = j;

    // number of neighbors used and allocated
    parts[j].NumNeighbors = 0;
    parts[j].AllocNeighbors = 1;

    // other status info
#ifdef MPI
    parts[j].NumReqs = 0;
#endif
    ClearReq(j);
    ClearLoad(j);
    ClearComp(j);

  }

}
//--------------------------------------------------------------------------
//
// ~Partition
// frees memory
//
Partition::~Partition()
{

  for (int j = 0; j < npart; j++) {

    for (int i = 0; i < parts[j].AllocNeighbors; i++) {
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
// GrowNeighbors
//
// doubles the size of the neighbor-dependent structures
//
// myrank: global partition number
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
  // 4 requests per neighbor
  assert((parts[myrank].Reqs =
	  (MPI_Request *)realloc(parts[myrank].Reqs,
				 nn * 4 * sizeof(MPI_Request) * 2)) != NULL);
#endif

  parts[myrank].AllocNeighbors *= 2;

}
//---------------------------------------------------------------------------

// MPI versions of SendNeighbors and ReceiveNeighbors 
// todo: make shared memory versions

#ifdef MPI


//---------------------------------------------------------------------------
//
// SendNeigbors()
//
// sends points to all neighbors
//
// myrank: global partition number
// ranks: ranks (global partition numbers) of all neighbors
// commm: MPI communicator
//
void Partition::SendNeighbors(int myrank, int *ranks, MPI_Comm comm) {

  int proc, myproc;
  int i, j;

  MPI_Comm_rank(comm, &myproc);

  int err;

//   // debug
//   if (myrank == 284 || myrank == 281) {
//     for (i = 0; i < parts[myrank].NumNeighbors; i++)
//       fprintf(stderr," rank %d has neighbor %d\n",myrank,ranks[i]);
//   }

  // for all neighbors
  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    // if neighbor exists (not beyond domain boundary)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // send only to remote locations
      if (proc != myproc) {

	j = parts[myrank].NumReqs++;

	// tag for header message offset by MAX_PARTS to differentiate
	// from points message
// 	MPI_Isend(&(parts[myrank].NumSendPoints[i]), 1, MPI_INT, proc, 
// 		 MAX_PARTS + ranks[i], comm, &(parts[myrank].Reqs[j]));

	err = MPI_Send(&(parts[myrank].NumSendPoints[i]), 1, MPI_INT, proc, 
		 MAX_PARTS + ranks[i], comm);

	// debug
	if (myrank == 284 && ranks[i] == 279)
	  fprintf(stderr,"rank 284 sent %d points to rank 279, proc %d tag %d err = %d\n",parts[myrank].NumSendPoints[i],proc,MAX_PARTS+ranks[i],err);
	if (myrank == 281 && ranks[i] == 272)
	  fprintf(stderr,"rank 281 sent %d points to rank 272, proc %d tag %d err = %d\n",parts[myrank].NumSendPoints[i],proc,MAX_PARTS+ranks[i],err);
	if (ranks[i] == 272)
	  fprintf(stderr,"-1:\n");
	if (ranks[i] == 279)
	  fprintf(stderr,"-2:\n");

	if (parts[myrank].NumSendPoints[i]) {
// 	  fprintf(stderr, "rank %d sending %d points to remote rank %d\n",
// 		  myrank, parts[myrank].NumSendPoints[i], ranks[i]);

	  j = parts[myrank].NumReqs++;
// 	  MPI_Isend(parts[myrank].SendPoints[i], 
// 		   parts[myrank].NumSendPoints[i] * 4, MPI_FLOAT, proc,
// 		    ranks[i], comm, &(parts[myrank].Reqs[j]));

	  MPI_Send(parts[myrank].SendPoints[i], 
		   parts[myrank].NumSendPoints[i] * 4, MPI_FLOAT, proc,
		    ranks[i], comm);

	  parts[myrank].NumSendPoints[i] = 0;
	}

      } // if location is remote

      // debug
//       else {
// 	if (parts[myrank].NumSendPoints[i])
// 	  fprintf(stderr, "rank %d sending %d points to local rank %d\n",
// 		  myrank, parts[myrank].NumSendPoints[i], ranks[i]);
// 	}

    } // if neighbor exists

  } // for all neighbors

}
//--------------------------------------------------------------------------
//
// ReceiveNeighbors
//
// receives points from all neighbors
//
// block_ranks: global partition numbers for blocks in this process
// neighbor_ranks: ranks of neighbors of all blocks for this process
// block: current block number
// nb: number of blocks in this process
// comm: MPI communicator
//
// returns total number of points received
//
int Partition::ReceiveNeighbors(int *block_ranks, int **neighbor_ranks,
				int block, int nb, MPI_Comm comm) {

  int nr; // rank of neighbor
  int neigh_block; // block number of neighbor
  int num = 0; // total number of points received
  int proc, myproc;  
  int i, j, k;
  int myrank = block_ranks[block];

  MPI_Comm_rank(comm, &myproc);

  // debug
  MPI_Status status;

  // for all neighbors
  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    nr = neighbor_ranks[block][i];

    // if neighbor exists (not outside of domain)
    if (nr >= 0) {

      proc = GetProc(nr);

      // remote: get number of received points and tag
      if (proc != myproc) {

	j = parts[myrank].NumReqs++;	
	// tag for header message offset by MAX_PARTS to differentiate 
	// from points message
// 	MPI_Irecv(&(parts[myrank].NumRecvPoints[i]), 1, MPI_INT, 
// 		  proc, MAX_PARTS + myrank, comm, &(parts[myrank].Reqs[j]));
// 	MPI_Waitall(parts[myrank].NumReqs, parts[myrank].Reqs, 
// 		    MPI_STATUSES_IGNORE);

	if (myrank == 279 && nr == 284)
	  fprintf(stderr,"myrank %d expecting num pts from nr %d proc %d tag %d\n",myrank,nr,proc,MAX_PARTS+myrank);

	if (myrank == 272 && nr == 281)
	  fprintf(stderr,"myrank %d expecting num pts from nr %d proc %d tag %d\n",myrank,nr,proc,MAX_PARTS+myrank);

	MPI_Recv(&(parts[myrank].NumRecvPoints[i]), 1, MPI_INT, 
		 proc, MAX_PARTS + myrank, comm, &status);

	if (myrank == 279 && nr == 284)
	fprintf(stderr,"myrank %d received num pts from nr %d\n",myrank,nr);

	if (myrank == 272 && nr == 281)
	fprintf(stderr,"myrank %d received num pts from nr %d\n",myrank,nr);


	if (myrank == 272)
	  fprintf(stderr,"1:\n");
	if (myrank == 279)
	  fprintf(stderr,"2:\n");

      }

      // local: get number of received points
      if (proc == myproc) {

	// block number of my neighbor
	// eventually use hash table or more efficient reverse lookup
	for (neigh_block = 0; neigh_block < nb; neigh_block++) {
	  if (block_ranks[neigh_block] == nr)
	    break;
	}
	assert(neigh_block < nb);

	// I am neighbor k of my neighbor
	for (k = 0; k < parts[nr].NumNeighbors; k++) {
	  if (neighbor_ranks[neigh_block][k] == myrank)
	    break;
	}

	// I must appear as a neighbor of my neighbor (symmetric)
	assert(k < parts[nr].NumNeighbors);

	// my number received = the number my neighbor sends
	parts[myrank].NumRecvPoints[i] = parts[nr].NumSendPoints[k];

      }

      // get the points

      // if something to receive
      if (parts[myrank].NumRecvPoints[i]) {

	while (parts[myrank].SizeRecvPoints[i] < 
            parts[myrank].NumRecvPoints[i] * 4 * sizeof(float)) {
	  parts[myrank].RecvPoints[i] = (float *)realloc(
            parts[myrank].RecvPoints[i], parts[myrank].SizeRecvPoints[i] * 2);
	  assert(parts[myrank].RecvPoints[i] != NULL);
	  parts[myrank].SizeRecvPoints[i] *= 2;

	}

	// remote
	if (proc != myproc) {

	  j = parts[myrank].NumReqs++;
// 	  MPI_Irecv(parts[myrank].RecvPoints[i],
// 		    parts[myrank].NumRecvPoints[i] * 4, MPI_FLOAT, proc, 
// 		    myrank, comm, &(parts[myrank].Reqs[j]));
// 	  MPI_Waitall(parts[myrank].NumReqs, parts[myrank].Reqs, 
// 		      MPI_STATUSES_IGNORE);

	  MPI_Recv(parts[myrank].RecvPoints[i],
		    parts[myrank].NumRecvPoints[i] * 4, MPI_FLOAT, proc, 
		   myrank, comm, &status);

	  // debug
// 	  fprintf(stderr, "rank %d received %d points from remote rank %d\n", myrank, parts[myrank].NumRecvPoints[i], nr);

	} // remote

	// local
	if (proc == myproc) {

	  for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++) {
	    parts[myrank].RecvPoints[i][4 * j + 0] = 
            parts[nr].SendPoints[k][4 * j + 0];
	    parts[myrank].RecvPoints[i][4 * j + 1] = 
            parts[nr].SendPoints[k][4 * j + 1];
	    parts[myrank].RecvPoints[i][4 * j + 2] = 
            parts[nr].SendPoints[k][4 * j + 2];
	    parts[myrank].RecvPoints[i][4 * j + 3] = 
            parts[nr].SendPoints[k][4 * j + 3];
	  }

	  parts[nr].NumSendPoints[k] = 0;

// 	  // debug
// 	  fprintf(stderr, "rank %d received %d points from local rank %d\n", myrank, parts[myrank].NumRecvPoints[i], nr);

	} // local


	num += parts[myrank].NumRecvPoints[i];

      } // something to receive

    } // if neighbor exists

  } // for all neighbors

  parts[myrank].NumReqs = 0;

  return num;

}
//------------------------------------------------------------------------

#endif
//---------------------------------------------------------------------------
//
// Check
//
// debug routine prints numbers of received points
//
// caller must ensure that list has enough room
//
//
void Partition::Check(int myrank) {

  int i;

  for (i = 0; i < parts[myrank].NumNeighbors; i++) {

    fprintf(stderr,"myrank = %d NumRecvPoints[%d] = %d\n",myrank,i,parts[myrank].NumRecvPoints[i]);

  }

}
//---------------------------------------------------------------------------
