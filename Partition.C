
#include "Partition.h"

//--------------------------------------------------------------------------
//
// Partition
//
// constructs and initializes the partitions
//
// nsp: total (global) number of spatial partitions
// ntp: total (global) number of temporal partitions
// d: actual dimensions used 3 (3D steady state) or 4(3D time varying)
//
Partition::Partition(int nsp, int ntp, int d) {

  int i, j;

  assert((npart  = nsp * ntp) <= MAX_PARTS);
  assert((parts = new Partition4D[npart]) != NULL);

  // neighborhood size
  if (d == 3)
    nbhd = 27;
  else
    nbhd = MAX_NEIGHBORS;

  // init the partitions list
  for (int j = 0; j < npart; j++)  {

    // default procs
    parts[j].Proc = j;

    // sending and receiving lists
    for (int i = 0; i < nbhd; i++){

      parts[j].NumSendPoints[i] = parts[j].NumRecvPoints[i] = 0;
      parts[j].SendPoints[i] = (float *)malloc(4 * sizeof(float));
      assert(parts[j].SendPoints[i] != NULL);
      parts[j].RecvPoints[i] = (float *)malloc(4 * sizeof(float));
      assert(parts[j].RecvPoints[i] != NULL);
      parts[j].SizeSendPoints[i] = parts[j].SizeRecvPoints[i] = 
	4 * sizeof(float);

    }

    // other status info
    parts[j].NumReqs = 0;
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

  if (nbhd == 27 || nbhd == MAX_NEIGHBORS) {

    for (int j = 0; j < npart; j++) {
      for (int i = 0; i < nbhd; i++) {
	free(parts[j].SendPoints[i]);
	free(parts[j].RecvPoints[i]);
      }
    }

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
// neighbor: number of the neighbor (0 - MAX_NEIGHBORS)
//
void Partition::PostPoint(int myrank, VECTOR4 p, int neighbor) {

  while (parts[myrank].SizeSendPoints[neighbor] < 
      (parts[myrank].NumSendPoints[neighbor] + 1) * 4 * sizeof(float)) {

    parts[myrank].SendPoints[neighbor] = (float *)realloc(
	parts[myrank].SendPoints[neighbor],
	parts[myrank].SizeSendPoints[neighbor] * 2);

    assert(parts[myrank].SendPoints[neighbor] != NULL);
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

  for (i = 0; i < MAX_NEIGHBORS; i++) {

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

  for (i = 0; i < MAX_NEIGHBORS; i++) {

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
  for (i = 0; i < MAX_NEIGHBORS; i++) {

    for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++)
      (ls[num++]).Set(
		      parts[myrank].RecvPoints[i][4 * j + 0], 
		      parts[myrank].RecvPoints[i][4 * j + 1], 
		      parts[myrank].RecvPoints[i][4 * j + 2],
		      parts[myrank].RecvPoints[i][4 * j + 3]);

  }

  // clear receive lists
  for (i = 0; i < MAX_NEIGHBORS; i++)
    parts[myrank].NumRecvPoints[i] = 0;

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
//
void Partition::SendNeighbors(int myrank, int *ranks, MPI_Comm comm) {

  int proc, myproc;
  int i, j;

  MPI_Comm_rank(comm, &myproc);

  // for all neighbors
  for (i = 0; i < MAX_NEIGHBORS; i++) {

    // if neighbor exists (not beyond domain boundary)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // send only to remote locations
      if (proc != myproc) {

	j = parts[myrank].NumReqs++;
	// tag for header message offset by MAX_PARTS to differentiate
	// from points message
	MPI_Isend(&(parts[myrank].NumSendPoints[i]), 1, MPI_INT, proc, 
		 MAX_PARTS + ranks[i], comm, &(parts[myrank].Reqs[j]));

	if (parts[myrank].NumSendPoints[i]) {
	  fprintf(stderr, "rank %d sending %d points to remote rank %d\n",
		  myrank, parts[myrank].NumSendPoints[i], ranks[i]);
	  j = parts[myrank].NumReqs++;
	  MPI_Isend(parts[myrank].SendPoints[i], 
		   parts[myrank].NumSendPoints[i] * 4, MPI_FLOAT, proc,
		    ranks[i], comm, &(parts[myrank].Reqs[j]));
	  parts[myrank].NumSendPoints[i] = 0;
	}

      } // if location is remote

      // debug
      else {
	if (parts[myrank].NumSendPoints[i])
	  fprintf(stderr, "rank %d sending %d points to local rank %d\n",
		  myrank, parts[myrank].NumSendPoints[i], ranks[i]);
	}

    } // if neighbor exists

  } // for all neighbors

}
//--------------------------------------------------------------------------
//
// ReceiveNeighbors
//
// receives points from all neighbors
//
// myrank: global partition  number
// ranks: ranks (global partition numbers) of all neighbors
//
// returns total number of points received
//
int Partition::ReceiveNeighbors(int myrank, int *ranks, MPI_Comm comm) {

  int num = 0;
  int proc, myproc;  
  int i, j, k;
  int n;

  MPI_Comm_rank(comm, &myproc);

  // for all neighbors
  for (i = 0; i < MAX_NEIGHBORS; i++) {

    // if neighbor exists (not outside of domain)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // remote: get number of received points and tag
      if (proc != myproc) {
	j = parts[myrank].NumReqs++;	
	// tag for header message offset by MAX_PARTS to differentiate
	// from points message
	MPI_Irecv(&(parts[myrank].NumRecvPoints[i]), 1, MPI_INT, 
		  proc, MAX_PARTS + myrank, comm, &(parts[myrank].Reqs[j]));

      }

      // local: get number of received points
      if (proc == myproc) {
	k = nbhd - 1 - i; // I am neighbor k of my neighbor
	if (!parts[ranks[i]].NumSendPoints[k])
	  k = -1;
	if (k >= 0)
	  parts[myrank].NumRecvPoints[i] = parts[ranks[i]].NumSendPoints[k];
	else
	  parts[myrank].NumRecvPoints[i] = 0;
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
	  MPI_Irecv(parts[myrank].RecvPoints[i],
		    parts[myrank].NumRecvPoints[i] * 4, MPI_FLOAT, proc, 
		    myrank, comm, &(parts[myrank].Reqs[j]));
	  // debug
// 	  fprintf(stderr, "rank %d received %d points from remote rank %d\n", myrank, parts[myrank].NumRecvPoints[i], ranks[i]);

	} // remote

	// local
	if (proc == myproc) {

	  for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++) {
	    parts[myrank].RecvPoints[i][4 * j + 0] = 
            parts[ranks[i]].SendPoints[k][4 * j + 0];
	    parts[myrank].RecvPoints[i][4 * j + 1] = 
            parts[ranks[i]].SendPoints[k][4 * j + 1];
	    parts[myrank].RecvPoints[i][4 * j + 2] = 
            parts[ranks[i]].SendPoints[k][4 * j + 2];
	    parts[myrank].RecvPoints[i][4 * j + 3] = 
            parts[ranks[i]].SendPoints[k][4 * j + 3];
	  }

	  parts[ranks[i]].NumSendPoints[k]  = 0;

// 	  // debug
// 	  fprintf(stderr, "rank %d received %d points from local rank %d\n", myrank, parts[myrank].NumRecvPoints[i], ranks[i]);

	} // local

	num += parts[myrank].NumRecvPoints[i];

      } // something to receive

    } // if neighbor exists

  } // for all neighbors

  // flush all pending messages
  MPI_Waitall(parts[myrank].NumReqs, parts[myrank].Reqs, 
	      MPI_STATUSES_IGNORE);
  parts[myrank].NumReqs = 0;

  return num;

}
//------------------------------------------------------------------------

#endif
