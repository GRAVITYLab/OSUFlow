
#include "Lattice4D.h"
//---------------------------------------------------------------------------
//
// Lattice 4D 
//
// constructs and initializes the time-varying lattice class
//
// xlen, ylen, zlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of parition in time 
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp) {

  int i, j;

  volume_bounds_type *vbs; 

  nbhd = 0;

  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  ldim = tlen;  // number of time steps
  // spatial domain partitioning first 
  vbs = calc_subvolume(xlen, ylen, zlen, ghost, nsp, idim, jdim, kdim); 
  tdim = ntp; 
  npart = nsp*ntp; 
  parts = new Partition4D[npart]; 
  flowMatrix = new int[npart*npart]; 
  memset(flowMatrix, '\0', npart*npart*sizeof(int)); 
  vb_list = new volume_bounds_type[npart]; 

  for (int t=0; t<tdim; t++)
    for(int k=0; k<kdim; k++)
      for (int j=0; j<jdim; j++)
	for (int i=0; i<idim; i++) {
	  int idx = t *idim*jdim*kdim+ k*idim*jdim + j*idim+i; 
	  int idx2 = k*idim*jdim + j*idim+i; 
	  vb_list[idx].xmin=vbs[idx2].xmin; 
	  vb_list[idx].xmax=vbs[idx2].xmax; 
	  vb_list[idx].ymin=vbs[idx2].ymin; 
	  vb_list[idx].ymax=vbs[idx2].ymax; 
	  vb_list[idx].zmin=vbs[idx2].zmin; 
	  vb_list[idx].zmax=vbs[idx2].zmax; 
	  vb_list[idx].tmin= tlen/ntp * t-1; 
	  if (vb_list[idx].tmin <0) vb_list[idx].tmin = 0; 
	  vb_list[idx].tmax= tlen/ntp * (t+1)+1; 
	  if (vb_list[idx].tmax > tlen-1) vb_list[idx].tmax = tlen-1; 
	}

  delete [] vbs; 

  // init the partitions list
  for (int j = 0; j < npart; j++)  {
    // assign default procs
    parts[j].Proc = j;
  }

  // init the block status
  for (i = 0; i < MAX_BLOCKS; i++) {
    ClearLoad(i);
    ClearComp(i);
  }

}
//--------------------------------------------------------------------------
//
// destructor
//
Lattice4D::~Lattice4D()
{

  // delete the MPI-related memory
#ifdef MPI
  if (nbhd == 27 || nbhd == MAX_NEIGHBORS) {

    for (int j = 0; j < npart; j++) {

      for (int i = 0; i < nbhd; i++) {
	free(parts[j].SendPoints[i]);
	free(parts[j].RecvPoints[i]);
      }

    }

  }
#endif

  if (flowMatrix!=NULL) delete [] flowMatrix; 
  if (parts!=NULL) delete [] parts; 
  if (vb_list!=NULL) delete [] vb_list; 
  if (seedlists!=NULL) delete [] seedlists; 

}
//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void Lattice4D::RoundRobin_proc(int nproc) {

  for (int i = 0; i < npart; i++) 
    parts[i].Proc = i % nproc; 

}
//---------------------------------------------------------------------------
//
// look up the lattice[i,j,k] element for the processor of the subdomain
//
int Lattice4D::GetProc(int i, int j, int k, int t) {

  if (i < 0 || i >= idim)
    return(-1); 
  else if (j < 0 || j >= jdim)
    return(-1); 
  else if (k < 0 || k >= kdim)
    return(-1); 
  else if (t < 0 || t >= tdim)
    return(-1); 

  int idx = t* idim*jdim*kdim + k * idim * jdim + j * idim + i; 

  return (parts[idx].Proc); 

}
//---------------------------------------------------------------------------
//
int Lattice4D::GetProc(int rank) {

  return(parts[rank].Proc); 

}
//----------------------------------------------------------------------------
//
// query the paritions that are assigned to processor 'proc' 
//
void Lattice4D::GetPartitions(int proc, int**p_list, int& num) {

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
// look up the lattice[i,j,k] element to get the subdomain number (rank)
// returns -1 if out of the domain
//
int Lattice4D::GetRank(int i, int j, int k, int t) {

  if (i < 0 || i >= idim || j < 0 || j >= jdim || k < 0 || k >= kdim || 
     t < 0 || t >= tdim) 
    return(-1); 

  return(t * idim * jdim * kdim + k * idim * jdim + j * idim + i);

}
//----------------------------------------------------------------------------
//
// Find the subdomain that contains the physical location (x,y,z) 
//
int Lattice4D::GetRank(float x, float y, float z, float t) {
  
  for (int i = 0; i < npart; i++) {

    if (vb_list[i].xmin > x || vb_list[i].xmax < x) 
      continue; 
    else if (vb_list[i].ymin > y || vb_list[i].ymax < y) 
      continue; 
    else if (vb_list[i].zmin > z || vb_list[i].zmax < z) 
      continue; 
    else if (vb_list[i].tmin > t || vb_list[i].tmax < t) 
      continue; 
    return(i); 
 
 }

  return(-1); 
}
//---------------------------------------------------------------------------
//
// find the indices of the lattice element that has its subdomain number = rank 
//
int Lattice4D::GetIndices(int rank, int &iidx, int &jidx, int &kidx, int& tidx) {
  
  if (rank < 0 || rank >= npart) 
    return(-1); 

  tidx = rank / (idim*jdim*kdim); 
  int r = rank % (idim*jdim*kdim); 
  kidx = r / (idim * jdim) ; 
  jidx = (r % (idim * jdim)) / idim; 
  iidx = r % idim; 

  return(1); 
}
//--------------------------------------------------------------------------
//
// find the indices of the lattice element that contanis (x,y,z) 
//
int Lattice4D::GetIndices(float x, float y, float z, float t, int &iidx, int &jidx, int &kidx, int&tidx) {

  int rank = GetRank(x, y, z, t); 

  if (rank != -1) {
    GetIndices(rank, iidx, jidx, kidx, tidx); 
    return(rank); 
  }
  else 
    return(-1); 

}
//--------------------------------------------------------------------------
//
// look up the volume bounds of lattice[i,j,k] 
//
int Lattice4D::GetBounds(int i, int j, int k, int t, volume_bounds_type& vb)  {

  if (i < 0 || i >= idim) 
    return(-1); 
  else if (j < 0 || j >= jdim) 
    return(-1); 
  else if (k < 0 || k >= kdim) 
    return(-1); 
  else if (t < 0 || t >=tdim) 
    return(-1); 
  
  int idx = t*idim*jdim*kdim + k * idim * jdim + j * idim + i; 
  vb = vb_list[idx]; 
  return(1); 

}
//--------------------------------------------------------------------------
//
// look up the volume bounds of the subdomain 'rank' 
//
int Lattice4D::GetBounds(int rank, volume_bounds_type &vb) {

  if (rank < 0 || rank >= npart) 
    return(-1); 
  vb = vb_list[rank]; 
  return(1); 

}
//--------------------------------------------------------------------------
//
// check if the point (x,y,z) is in the lattice element [i,j,k,l]
//
bool Lattice4D::isIn(float x, float y, float z, float t, int i, int j, int k, 
		   int l) {

  if (i < 0 || i > idim - 1 || j < 0 || j > jdim - 1 || k < 0 || k > kdim - 1 || l <0 || l> tdim-1) 
    return(false); 

  int idx = l*idim*jdim*kdim + k * idim * jdim + j * idim + i; 
  volume_bounds_type vb = vb_list[idx]; 

  if (vb.xmin > x || vb.xmax < x) 
    return (false); 
  else if (vb.ymin > y || vb.ymax < y) 
    return (false); 
  else if (vb.zmin > z || vb.zmax < z) 
    return (false); 
  else if (vb.tmin > t || vb.tmax < t) 
    return(false); 

  return(true); 

}
//--------------------------------------------------------------------------
//
// returns neighbor rank where x,y,z,t is in 
// this is a global partition number
// returns -1 if out of the domain
//
int Lattice4D::CheckNeighbor(int myrank, float x, float y, float z, float t) {

  int i,j,k,l; 

  GetIndices(myrank, i, j, k, l); 

  for (int nt=-1; nt<=1; nt++) 
    for (int nk=-1; nk<=1; nk++) 
      for (int nj=-1; nj<=1; nj++) 
	for (int ni=-1; ni<=1; ni++)  {
	  if (nt==0 && nk==0 && nj==0 && ni==0) continue; //self 
	  else if (nt+l<0 || nt+l>tdim-1) continue; 
	  else if (nk+k<0 || nk+k>kdim-1) continue;
	  else if (nj+j<0 || nj+j>jdim-1) continue;
	  else if (ni+i<0 || ni+i>idim-1) continue;
	  else if (isIn(x,y,z,t,i+ni, j+nj, k+nk, l+nt))
	    return(GetRank(ni+i, nj+j, nk+k, nt+l)); 
	}
  return(-1); 
}

//
// GetNeighbor
//
//
int Lattice4D::GetNeighbor(int myrank, float x, float y, float z, float t, 
                           int &ei, int &ej, int &ek, int &et) {

  int neighbor;
  int si, sj, sk;

  neighbor = CheckNeighbor(myrank, x, y, z, t);
  
  GetIndices(neighbor, ei, ej, ek, et); 
  
  return (neighbor); 

}


//---------------------------------------------------------------------------

void Lattice4D::InitSeedLists() {

  seedlists = new list<VECTOR4>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

void Lattice4D::ResetSeedLists() {

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}

void Lattice4D::ResetSeedLists(int i) {

  seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

bool Lattice4D::InsertSeed(int i, int j, int k, int t, VECTOR4 p) {

  int rank = GetRank(i,j,k, t); 

  if (rank ==-1) return(false); 
  else {
    seedlists[rank].push_back(p); 
    return(true); 
  }
}
//--------------------------------------------------------------------------

bool Lattice4D::InsertSeed(int from_i, int from_j, int from_k, int from_t, 
			   int i, int j, int k, int t, VECTOR4 p) {

  int from_rank = GetRank(from_i, from_j, from_k, from_t); 
  int to_rank = GetRank(i,j,k, t); 

  if (to_rank ==-1 || from_rank==-1) return(false); 
  else {
    seedlists[to_rank].push_back(p); 
    flowMatrix[from_rank*npart+to_rank]++; 
    return(true); 
  }
}
//--------------------------------------------------------------------------

bool Lattice4D::InsertSeed(int i, VECTOR4 p) {

  if (i>=npart) return(false); 
  else {
    seedlists[i].push_back(p); 
    return(true); 
  }
}
//--------------------------------------------------------------------------

bool Lattice4D::InsertSeed(int from_rank, int to_rank, VECTOR4 p) {

  if (from_rank >=npart || to_rank>=npart) return(false); 
  else {
    flowMatrix[from_rank*npart+to_rank]++; 
    seedlists[to_rank].push_back(p); 
  }    return(true); 
}

void Lattice4D::ResetFlowMatrix() 
{
  if (flowMatrix !=NULL) 
    memset(flowMatrix, '\0', npart*npart*sizeof(int)); 
}
//--------------------------------------------------------------------------
//
// MPI functions from this point on
// author: Tom Peterka
//

#ifdef MPI

//--------------------------------------------------------------------------
//
// Lattice 4D 
//
// constructs and initializes the time-varying lattice class
//
// xlen, ylen, zlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of parition in time 
// d: actual dimensions used 3 (3D steady state) or 4(3D time varying)
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp, int d) {

  int i, j;
  volume_bounds_type *vbs; 

  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  ldim = tlen;  // number of time steps

  // spatial domain partitioning first 
  vbs = calc_subvolume(xlen, ylen, zlen, ghost, nsp, idim, jdim, kdim); 
  tdim = ntp; 
  npart = nsp*ntp; 
  parts = new Partition4D[npart]; 
  flowMatrix = new int[npart*npart]; 
  memset(flowMatrix, '\0', npart*npart*sizeof(int)); 
  vb_list = new volume_bounds_type[npart]; 

  for (int t=0; t<tdim; t++)
    for(int k=0; k<kdim; k++)
      for (int j=0; j<jdim; j++)
	for (int i=0; i<idim; i++) {
	  int idx = t *idim*jdim*kdim+ k*idim*jdim + j*idim+i; 
	  int idx2 = k*idim*jdim + j*idim+i; 
	  vb_list[idx].xmin=vbs[idx2].xmin; 
	  vb_list[idx].xmax=vbs[idx2].xmax; 
	  vb_list[idx].ymin=vbs[idx2].ymin; 
	  vb_list[idx].ymax=vbs[idx2].ymax; 
	  vb_list[idx].zmin=vbs[idx2].zmin; 
	  vb_list[idx].zmax=vbs[idx2].zmax; 
	  vb_list[idx].tmin= tlen/ntp * t-1; 
	  if (vb_list[idx].tmin <0) vb_list[idx].tmin = 0; 
	  vb_list[idx].tmax= tlen/ntp * (t+1)+1; 
	  if (vb_list[idx].tmax > tlen-1) vb_list[idx].tmax = tlen-1; 
	}

  delete [] vbs; 

  // neighborhood size
  if (d == 3)
    nbhd = 27;
  else
    nbhd = MAX_NEIGHBORS;

  // init the partitions list
  for (int j = 0; j < npart; j++)  {

    // assign default procs
    parts[j].Proc = j;

    // init sending and receiving lists, request lists
    for (int i = 0; i < nbhd; i++){

      parts[j].NumSendPoints[i] = parts[j].NumRecvPoints[i] = 0;
      parts[j].SendPoints[i] = (float *)malloc(4 * sizeof(float));
      assert(parts[j].SendPoints[i] != NULL);
      parts[j].RecvPoints[i] = (float *)malloc(4 * sizeof(float));
      assert(parts[j].RecvPoints[i] != NULL);
      parts[j].SizeSendPoints[i] = parts[j].SizeRecvPoints[i] = 4 * sizeof(float);
      parts[j].HasData = 0;
      parts[j].NumReqs = 0;

    }

  }

  // init the block status
  for (i = 0; i < MAX_BLOCKS; i++) {
    ClearLoad(i);
    ClearComp(i);
  }

}
//--------------------------------------------------------------------------
//
// GetNumPartitions
//
// return the number of partitions that are assigned to processor 'proc' 
//
int Lattice4D::GetNumPartitions(int proc) {

  int n = 0; 
  int i;

  for (i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      n++;
  }

  return n;

}
//---------------------------------------------------------------------------
//
// GetPartitions
//
// query the partitions that are assigned to processor 'proc' 
// p_list must be allocated large enough prior to calling
//
void Lattice4D::GetPartitions(int proc, int*p_list) {

  int n = 0; 
  int i;

  for (i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      p_list[n++] = i;
  }

  return;

}
//---------------------------------------------------------------------------
//
// GetNeighbor
// returns neighbor number (0 - MAX_NEIGHBORS) of neighbor containing point
// returns -1 if point is not in one of the neighbors
//
int Lattice4D::GetNeighbor(int myrank, float x, float y, float z, float t) {

  int n; // neighbor number (0 - MAX_NEIGHBORS)
  int nr; // my neighbor's rank (global partition) number
  int r; // rank of point

  r = GetRank(x, y, z, t);

  // for all neighbors
  for (n = 0; n < nbhd; n++) {

    // offset my rank to get neighbor's rank
    nr = myrank + n - (nbhd - 1) / 2;

    if (r == nr)
      return n;

  }

  return -1;

}
//---------------------------------------------------------------------------
//
// GetNeighborRanks
//
// gets ranks of all neighbors
// ranks are the global partition numbers for all neighbors
// the neighbor of an edge gets rank -1
//
void Lattice4D::GetNeighborRanks(int myrank, int *neighbor_ranks) {

  int in, jn, kn, ln; // my neighbor's lattice coords
  int n; // my neighbor number (0-80)
  int nr; // my neighbor's rank (global partition) number

  // for all neighbors
  for (n = 0; n < nbhd; n++) {

    // offset my rank to get neighbor's rank
    nr = myrank + n - (nbhd - 1) / 2;

    // use GetIndices to filter out ranks outside of boundary
    // then convert the remaining indices back to rank
    if (GetIndices(nr, in, jn, kn, ln) == 1)
      neighbor_ranks[n] = GetRank(in, jn, kn, ln);
    else
      neighbor_ranks[n] = -1;

  }

}
//---------------------------------------------------------------------------
//
// PostPoint
//
// posts a point for sending to a neighbor
// myrank: my global partition number
// p: 4D point
// neighbor: number of the neighbor (0 - MAX_NEIGHBORS)
//
void Lattice4D::PostPoint(int myrank, VECTOR4 p, int neighbor) {

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
void Lattice4D::PrintPost(int myrank) {

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
void Lattice4D::PrintRecv(int myrank) {

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
void Lattice4D::GetRecvPts(int myrank, VECTOR4 *ls) {

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
//
// Error()
// error handler
//
void Lattice4D::Error(const char *fmt, ...){

  va_list argp;
  vfprintf(stderr, fmt, argp);
  sleep(5);
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD,0);
#else
  exit(0);
#endif

}
//---------------------------------------------------------------------------
//
// SendNeigbors()
//
// myrank: global partition number
// sends points to all neighbors
//
void Lattice4D::SendNeighbors(int myrank, MPI_Comm comm) {

  int ranks[MAX_NEIGHBORS];
  int proc, myproc;
  int i, j;

  MPI_Comm_rank(comm, &myproc);
  GetNeighborRanks(myrank, ranks);

  // for all neighbors
  for (i = 0; i < MAX_NEIGHBORS; i++) {

    // if neighbor exists (not beyond domain boundary)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // send only to remote locations
      if (proc != myproc) {

	j = parts[myrank].NumReqs++;
	MPI_Isend(&(parts[myrank].NumSendPoints[i]), 1, MPI_INT, proc, 
		 ranks[i], comm, &(parts[myrank].Reqs[j]));

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
// returns total number of points received
//
int Lattice4D::ReceiveNeighbors(int myrank, MPI_Comm comm) {

  int ranks[MAX_NEIGHBORS]; // ranks of my neighbors
  int num = 0;
  int proc, myproc;  
  int i, j, k;
  int n;

  MPI_Comm_rank(comm, &myproc);
  GetNeighborRanks(myrank, ranks);

  // for all neighbors
  for (i = 0; i < MAX_NEIGHBORS; i++) {

    // if neighbor exists (not outside of domain)
    if (ranks[i] >= 0) {

      proc = GetProc(ranks[i]);

      // remote: get number of received points and tag
      if (proc != myproc) {
	j = parts[myrank].NumReqs++;	
	MPI_Irecv(&(parts[myrank].NumRecvPoints[i]), 1, MPI_INT, 
		  proc, myrank, comm, &(parts[myrank].Reqs[j]));
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
