
#include "Lattice.h"

//---------------------------------------------------------------------------
//
// Lattice
//
// constructs and initializes the lattice class
//
// xlen, ylen, zlen: total size of the data
// ghost: size of ghost layer per side
// np: total number of partitions in the lattice
//
// Tom Peterka, 12/8/08
//
Lattice::Lattice(int xlen, int ylen, int zlen, int ghost, int np) {

  int i, j;

  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  vb_list = calc_subvolume(xlen, ylen, zlen, ghost, np, idim, jdim, kdim); 
  npart = np; 
  parts = new LatPartition[np]; 

  // init the partitions list
  for (int j = 0; j < np; j++)  {

    // assign default procs
    parts[j].Proc = j;

    // init sending and receiving lists
    for (int i = 0; i < 6; i++){

      parts[j].NumSendPoints[i] = parts[j].NumRecvPoints[i] = 0;
      if ((parts[j].SendPoints[i] = (float *)malloc(3 * sizeof(float))) == NULL)
	Error("Error: Lattice() cannot allocate memory for SendPoints\n");
      if ((parts[j].RecvPoints[i] = (float *)malloc(3 * sizeof(float))) == NULL)
	Error("Error: Lattice() cannot allocate memory for RecvPoints\n");
      parts[j].SizeSendPoints[i] = parts[j].SizeRecvPoints[i] = 3 * sizeof(float);

    }

  }

}
//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void Lattice::RoundRobin_proc(int nproc) {

  for (int i = 0; i < npart; i++) 
    parts[i].Proc = i % nproc; 

}
//---------------------------------------------------------------------------
//
// look up the lattice[i,j,k] element to get the processor for the subdomain number (rank)
//
int Lattice::GetProc(int i, int j, int k) {

  if (i < 0 || i >= idim)
    return(-1); 
  else if (j < 0 || j >= jdim)
    return(-1); 
  else if (k < 0 || k >= kdim)
    return(-1); 

  int idx = k * idim * jdim + j * idim + i; 

  return (parts[idx].Proc); 

}
//---------------------------------------------------------------------------
//
int Lattice::GetProc(int rank) {

  return(parts[rank].Proc); 

}
//----------------------------------------------------------------------------
//
// query the paritions that are assigned to processor 'proc' 
//
void Lattice::GetPartitions(int proc, int**p_list, int& num) {

  int cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      cnt++;
  }

  num = cnt; 
  (*p_list) = new int[cnt]; 
  cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      (*p_list)[cnt++] = i; 
  }

}
//----------------------------------------------------------------------------
//
// GetPartitions
//
// query the paritions that are assigned to processor 'proc' 
//
// p_list must be allocated large enough prior to calling
//
// Tom Peterka, 12/5/08
//
void Lattice::GetPartitions(int proc, int*p_list, int& num) {

  int cnt = 0; 

  for (int i = 0; i < npart; i++)
    if (parts[i].Proc == proc)
      cnt ++; 

  num = cnt; 
  cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      p_list[cnt++] = i;
  }

}
//---------------------------------------------------------------------------
//
// look up the lattice[i,j,k] element to get the subdomain number (rank)
//
int Lattice::GetRank(int i, int j, int k) {

  if (i < 0 || i >= idim) 
    return(-1); 
  else if (j < 0 || j >= jdim) 
    return(-1); 
  else if (k < 0 || k >= kdim) 
    return(-1); 

  int idx = k * idim * jdim + j * idim + i; 

  return (idx); 

}
//----------------------------------------------------------------------------
//
// GetNeighborRanks
//
// gets ranks of 6 neighbors
// the neighbor of an edge gets rank -1
//
// Tom Peterka, 11/26/08
//
void Lattice::GetNeighborRanks(int myrank, int *neighbor_ranks) {

  int i, j, k; // my lattice coords

  GetIndices(myrank, i, j, k);

  // -x, +x
  neighbor_ranks[0] = (i == 0        ? -1 : GetRank(i - 1, j, k));
  neighbor_ranks[1] = (i == xdim - 1 ? -1 : GetRank(i + 1, j, k));

  // -y, +y	    
  neighbor_ranks[2] = (j == 0        ? -1 : GetRank(i, j - 1, k));
  neighbor_ranks[3] = (j == ydim - 1 ? -1 : GetRank(i, j + 1, k));

  // -z, +z
  neighbor_ranks[4] = (k == 0        ? -1 : GetRank(i, j, k - 1));
  neighbor_ranks[5] = (k == zdim - 1 ? -1 : GetRank(i, j, k + 1));

}
//---------------------------------------------------------------------------
//
// Find the subdomain that contains the physical location (x,y,z) 
//
int Lattice::GetRank(float x, float y, float z) {
  
  for (int i = 0; i < npart; i++) {

    if (vb_list[i].xmin > x || vb_list[i].xmax < x) 
      continue; 
    else if (vb_list[i].ymin > y || vb_list[i].ymax < y) 
      continue; 
    else if (vb_list[i].zmin > z || vb_list[i].zmax < z) 
      continue; 
    return(i); 
 
 }

  return(-1); 

}
//---------------------------------------------------------------------------
//
// find the indices of the lattice element that has its subdomain number = rank 
//
int Lattice::GetIndices(int rank, int &iidx, int &jidx, int &kidx) {
  
  if (rank < 0 || rank >= npart) 
    return(-1); 

  kidx = rank / (idim * jdim) ; 
  jidx = (rank % (idim * jdim)) / idim; 
  iidx = rank % idim; 

  return(1); 

}
//--------------------------------------------------------------------------
//
// find the indices of the lattice element that contanis (x,y,z) 
//
int Lattice::GetIndices(float x, float y, float z, int &iidx, int &jidx, int &kidx) {

  int rank = GetRank(x, y, z); 

  if (rank != -1) {
    GetIndices(rank, iidx, jidx, kidx); 
    return(rank); 
  }
  else 
    return(-1); 

}
//--------------------------------------------------------------------------
//
// look up the volume bounds of lattice[i,j,k] 
//
int Lattice::GetBounds(int i, int j, int k, volume_bounds_type& vb)  {

  if (i < 0 || i >= idim) 
    return(-1); 
  else if (j < 0 || j >= jdim) 
    return(-1); 
  else if (k < 0 || k >= kdim) 
    return(-1); 
  
  int idx = k * idim * jdim + j * idim + i; 
  vb = vb_list[idx]; 
  return(1); 

}
//--------------------------------------------------------------------------
//
// look up the volume bounds of the subdomain 'rank' 
//
int Lattice::GetBounds(int rank, volume_bounds_type &vb) {

  if (rank < 0 || rank >= npart) 
    return(-1); 
  vb = vb_list[rank]; 
  return(1); 

}
//--------------------------------------------------------------------------
//
//check if the point (x,y,z) is in the lattice element [i,j,k]
//
bool Lattice::isIn(float x, float y, float z, int i, int j, int k) {

  if (i < 0 || i > idim - 1 || j < 0 || j > jdim - 1 || k < 0 || k > kdim - 1) 
    return(false); 

  int idx = k * idim * jdim + j * idim + i; 
  volume_bounds_type vb = vb_list[idx]; 

  if (vb.xmin > x || vb.xmax < x) 
    return (false); 
  else if (vb.ymin > y || vb.ymax < y) 
    return (false); 
  else if (vb.zmin > z || vb.zmax < z) 
    return (false); 

  return(true); 

}
//--------------------------------------------------------------------------
//
// returns neighbor number (0-5) of x,y,z point
//
int Lattice::CheckNeighbor(int myrank, float x, float y, float z) {

  int i,j,k; 

  GetIndices(myrank, i, j, k); 

  if (isIn(x, y, z, i-1, j, k) == true) return(0);  //-X
  if (isIn(x, y, z, i+1, j, k) == true) return(1);  //+X
  if (isIn(x, y, z, i, j-1, k) == true) return(2);  //-Y
  if (isIn(x, y, z, i, j+1, k) == true) return(3);  //+Y
  if (isIn(x, y, z, i, j, k-1) == true) return(4);  //-Z
  if (isIn(x, y, z, i, j, k+1) == true) return(5);  //+Z

  return(-1); 

}
//---------------------------------------------------------------------------
//
// GetNeighbor
//
// returns neighbor number (0-5) of x,y,z point and sets i,j,k of neighbor
// myrank: global partition number
//
// Tom Peterka, 11/26/08
//
int Lattice::GetNeighbor(int myrank, float x, float y, float z, int &ei, int &ej, int &ek) {

  int neighbor;
  int si, sj, sk;

  neighbor = CheckNeighbor(myrank, x, y, z);
  GetIndices(myrank, si, sj, sk);

  if (neighbor == 0) {
    ei = si - 1; ej = sj; ek = sk;
  }
  else if (neighbor == 1) {
    ei = si + 1; ej = sj; ek = sk;
  }
  else if (neighbor == 2) {
    ei = si; ej = sj - 1; ek = sk;
  }
  else if (neighbor == 3) {
    ei = si; ej = sj + 1; ek = sk;
  }
  else if (neighbor == 4) {
    ei = si; ej = sj; ek = sk - 1;
  }
  else if (neighbor == 5) {
    ei = si; ej = sj; ek = sk + 1;
  }

  return neighbor;

}
//----------------------------------------------------------------------------

void Lattice::InitSeedLists() {

  seedlists = new list<VECTOR3>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

void Lattice::ResetSeedLists() {

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

void Lattice::InsertSeed(int i, int j, int k, VECTOR3 p) {

  int rank = GetRank(i,j,k); 

  seedlists[rank].push_back(p); 

}

void Lattice::ClearSeedList(int rank) {
  seedlists[rank].clear(); 
}
//-------------------------------------------------------------------------
//
// PostPoint
//
// posts a point for sending to a neighbor
// myrank: global partition number
//
// Tom Peterka, 12/4/08
//
void Lattice::PostPoint(int myrank, VECTOR3 p, int neighbor) {

  while (parts[myrank].SizeSendPoints[neighbor] < 
      (parts[myrank].NumSendPoints[neighbor] + 1) * 3 * sizeof(float)) {

    parts[myrank].SendPoints[neighbor] = (float *)realloc(
	parts[myrank].SendPoints[neighbor],
	parts[myrank].SizeSendPoints[neighbor] * 2);

    if (parts[myrank].SendPoints[neighbor] == NULL)
      Error("Error: PostPoint() cannot reallocate memory\n");

    parts[myrank].SizeSendPoints[neighbor] *= 2;

  }

  parts[myrank].SendPoints[neighbor][3 * 
       parts[myrank].NumSendPoints[neighbor] + 0] = p[0];
  parts[myrank].SendPoints[neighbor][3 * 
       parts[myrank].NumSendPoints[neighbor] + 1] = p[1];
  parts[myrank].SendPoints[neighbor][3 * 
       parts[myrank].NumSendPoints[neighbor] + 2] = p[2];

  parts[myrank].NumSendPoints[neighbor]++;

}
//------------------------------------------------------------------------
//
// PrintPost
//
// prints the posted points
// myrank: global partition number
// 
// Tom Peterka, 12/4/08
//
void Lattice::PrintPost(int myrank) {

  int ranks[6];
  int i, j;

  GetNeighborRanks(myrank, ranks);

  fprintf(stderr, "\nPosted messages list (gp = global partition #)\n");

  for (i = 0; i < 6; i++) {

    if (parts[myrank].NumSendPoints[i])
      fprintf(stderr, "gp %d posted %d points to gp %d\n", 
	      myrank, parts[myrank].NumSendPoints[i], ranks[i]);

    if (parts[myrank].NumSendPoints[i]) {
      for (j = 0; j < parts[myrank].NumSendPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\n",
            parts[myrank].SendPoints[i][3 * j + 0],
	    parts[myrank].SendPoints[i][3 * j + 1], 
            parts[myrank].SendPoints[i][3 * j + 2]);
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
// Tom Peterka, 12/4/08
//
void Lattice::PrintRecv(int myrank) {

  int ranks[6];
  int i, j;

  GetNeighborRanks(myrank, ranks);

  fprintf(stderr, "\nReceived points list (gp = global partition #)\n");

  for (i = 0; i < 6; i++) {

    if (parts[myrank].NumRecvPoints[i]) {
      fprintf(stderr, 
           "gp %d received %d points from gp %d:\n", 
            myrank, parts[myrank].NumRecvPoints[i], ranks[i]);
      for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\n",
		parts[myrank].RecvPoints[i][3 * j + 0],
		parts[myrank].RecvPoints[i][3 * j + 1],
		parts[myrank].RecvPoints[i][3 * j + 2]);
    }

  }

  fprintf(stderr,"\n");

}
//---------------------------------------------------------------------------
//
// GetNumRecv
//
// returns number of received points ready for processing
// myrank: global partition number
//
// Tom Peterka, 12/4/08
//
int Lattice::GetNumRecv(int myrank) {

  int num = 0;
  int i;

  for (i = 0; i < 6; i++) {
    if (parts[myrank].NumRecvPoints[i] > 0)
      num += parts[myrank].NumRecvPoints[i];
  }

  return num;

}
//--------------------------------------------------------------------------
//
// GetRecvPts
//
// copies received points from partition to user supplied list
// myrank: global partition number
//
// caller must ensure that list has enough room
//
//
// Tom Peterka, 12/4/08
//
void Lattice::GetRecvPts(int myrank, VECTOR3 *ls) {

  int num = 0;
  int i, j;

  // copy points
  for (i = 0; i < 6; i++) {
    for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++)
      (ls[num++]).Set(parts[myrank].RecvPoints[i][3 * j + 0], 
	   parts[myrank].RecvPoints[i][3 * j + 1], 
           parts[myrank].RecvPoints[i][3 * j + 2]);
  }

  // clear receive lists
  for (i = 0; i < 6; i++)
    parts[myrank].NumRecvPoints[i] = 0;

}
//---------------------------------------------------------------------------

// MPI functions

#ifdef MPI

//---------------------------------------------------------------------------
//
// SendNeigbors()
//
// myrank: global partition number
// sends points to all neighbors
//
// Tom Peterka, 12/4/08
//
void Lattice::SendNeighbors(int myrank, MPI_Comm comm) {

  MPI_Request req;
  int ranks[6];
  int proc, myproc;
  int i;

  MPI_Comm_rank(comm, &myproc);
  GetNeighborRanks(myrank, ranks);

  // for all neighbors
  for (i = 0; i < 6; i++) {

    // if neighbor exists (not edge)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // send only to remote locations
      if (proc != myproc) {

	// if there are points to send
	if (parts[myrank].NumSendPoints[i]) {

	  MPI_Isend(parts[myrank].SendPoints[i], 
	    parts[myrank].NumSendPoints[i] * 3, MPI_FLOAT, proc, 0, comm, &req);

	  parts[myrank].NumSendPoints[i] = 0;

	} // if there are points to send

      } // if location is remote

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
// returns total number of points received
//
// Tom Peterka, 12/4/08
//
int Lattice::ReceiveNeighbors(int myrank, MPI_Comm comm) {

  MPI_Status status;
  MPI_Request req;
  int flag;
  int ready;
  int ranks[6]; // ranks of my neighbors
  int num = 0;
  int proc, myproc;  
  int i, j, k;
  int n;

  MPI_Comm_rank(comm, &myproc);
  GetNeighborRanks(myrank, ranks);

  // for all neighbors
  for (i = 0; i < 6; i++) {

    // if neighbor exists (not outside of domain)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // get number of received points

      // remote
      if (proc != myproc) {

	MPI_Iprobe(proc, 0, comm, &flag, &status);
	if (flag) {
	  MPI_Get_count(&status, MPI_FLOAT, &n);
	  parts[myrank].NumRecvPoints[i] = (float)n / 3.0f;
	}
	else
	  parts[myrank].NumRecvPoints[i] = 0;

      } // remote

      // local
      if (proc == myproc) {

	if (i % 2 == 0)
	  k = i + 1;
	else
	  k = i - 1;
	if (!parts[ranks[i]].NumSendPoints[k])
	  k = -1;
	if (k >= 0)
	  parts[myrank].NumRecvPoints[i] = parts[ranks[i]].NumSendPoints[k];
	else
	  parts[myrank].NumRecvPoints[i] = 0;

      } // local

      // get the points

      // if something to receive
      if (parts[myrank].NumRecvPoints[i]) {

	while (parts[myrank].SizeRecvPoints[i] < 
            parts[myrank].NumRecvPoints[i] * 3 * sizeof(float)) {
	  parts[myrank].RecvPoints[i] = (float *)realloc(
            parts[myrank].RecvPoints[i], parts[myrank].SizeRecvPoints[i] * 2);
	  if (parts[myrank].RecvPoints[i] == NULL)
	    Error("Error: ReceivePoints() cannot reallocate memory\n");
	  parts[myrank].SizeRecvPoints[i] *= 2;

	}

	// remote
	if (proc != myproc) {
	  MPI_Recv(parts[myrank].RecvPoints[i], 
             parts[myrank].NumRecvPoints[i] * 3, 
             MPI_FLOAT, proc, 0, comm, &status);
	} // remote

	// local
	if (proc == myproc) {

	  for (j = 0; j < parts[myrank].NumRecvPoints[i]; j++) {
	    parts[myrank].RecvPoints[i][3 * j + 0] = 
            parts[ranks[i]].SendPoints[k][3 * j + 0];
	    parts[myrank].RecvPoints[i][3 * j + 1] = 
            parts[ranks[i]].SendPoints[k][3 * j + 1];
	    parts[myrank].RecvPoints[i][3 * j + 2] = 
            parts[ranks[i]].SendPoints[k][3 * j + 2];
	  }

	  parts[ranks[i]].NumSendPoints[k]  = 0;

	} // local

	num += parts[myrank].NumRecvPoints[i];

      } // something to receive

    } // if neighbor exists

  } // for all neighbors

  return num;

}
//------------------------------------------------------------------------



#endif

//
// Error()
// error handler
//
void Lattice::Error(const char *fmt, ...){

}
