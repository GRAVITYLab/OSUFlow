
#include "Lattice4D.h"

//---------------------------------------------------------------------------
//
// constructs and initializes the time-varying lattice class
//
// xlen, ylen, zlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of parition in time 
// myid: rank, process number, thread number, identification of the owner
// default = -1 (can omit if single process sequential program)
// nid: number of processes, threads, owners
// default = -1 (can omit if single process sequential program)
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp, int myid, int nid) {

  int i, j, n;
  volume_bounds_type *vbs; 

  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  ldim = tlen;  // number of time steps

  // spatial domain partitioning first 
  vbs = calc_subvolume(xlen, ylen, zlen, ghost, nsp, idim, jdim, kdim); 
  tdim = ntp; 
  npart = nsp * ntp; 
  flowMatrix = new int[npart * npart]; 
  memset(flowMatrix, '\0', npart * npart * sizeof(int)); 
  vb_list = new volume_bounds_type[npart]; 

  for (int t = 0; t < tdim; t++) {

    for(int k = 0; k < kdim; k++) {

      for (int j = 0; j < jdim; j++) {

	for (int i = 0; i < idim; i++) {

	  int idx = t * idim * jdim * kdim + k * idim * jdim + j * idim + i; 
	  int idx2 = k * idim * jdim + j * idim + i; 
	  vb_list[idx].xmin = vbs[idx2].xmin; 
	  vb_list[idx].xmax = vbs[idx2].xmax; 
	  vb_list[idx].ymin = vbs[idx2].ymin; 
	  vb_list[idx].ymax = vbs[idx2].ymax; 
	  vb_list[idx].zmin = vbs[idx2].zmin; 
	  vb_list[idx].zmax = vbs[idx2].zmax; 
	  vb_list[idx].tmin = tlen / ntp * t - 1; 
	  if (vb_list[idx].tmin < 0)
	    vb_list[idx].tmin = 0; 
	  vb_list[idx].tmax = tlen / ntp * (t + 1) + 1; 
	  if (vb_list[idx].tmax > tlen - 1)
	    vb_list[idx].tmax = tlen - 1; 

	}

      }

    }

  }

  // neighborhood size
  nbhd = 81;

  // table of neighbor ranks
  assert((neighbor_ranks = (int *)malloc(nbhd * sizeof (int))) != NULL);

  // create partition class
  part = new Partition(nsp, ntp);

  delete [] vbs; 

  // ranks of my blocks
  myproc = myid;
  nproc = nid;
  if (myproc >= 0) {
    RoundRobin_proc(nproc);
    n = GetNumPartitions(myproc);
    assert((block_ranks = (int *)malloc(n * sizeof(int))) != NULL);
    GetPartitions(myproc, block_ranks);
  }

}
//---------------------------------------------------------------------------
//
// destructor
//
Lattice4D::~Lattice4D()
{

  if (flowMatrix!=NULL) delete [] flowMatrix; 
  if (vb_list!=NULL) delete [] vb_list; 
  if (seedlists!=NULL) delete [] seedlists; 

  delete part;
  delete block_ranks;

}
//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void Lattice4D::RoundRobin_proc(int nproc) {

  for (int i = 0; i < npart; i++) 
    part->parts[i].Proc = i % nproc; 

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

  return (part->parts[idx].Proc); 

}
//---------------------------------------------------------------------------
//
int Lattice4D::GetProc(int rank) {

  return(part->parts[rank].Proc); 

}
//----------------------------------------------------------------------------
//
// query the paritions that are assigned to processor 'proc' 
//
void Lattice4D::GetPartitions(int proc, int**p_list, int& num) {

  int cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (part->parts[i].Proc == proc)
      cnt++;
  }

  num = cnt; 
  (*p_list) = new int[cnt]; 
  cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (part->parts[i].Proc == proc)
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
//--------------------------------------------------------------------------
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
//
void Lattice4D::InitSeedLists() {

  seedlists = new list<VECTOR4>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------
//
void Lattice4D::ResetSeedLists() {

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------
//
void Lattice4D::ResetSeedLists(int i) {

  seedlists[i].clear(); 

}
//--------------------------------------------------------------------------
//
bool Lattice4D::InsertSeed(int i, int j, int k, int t, VECTOR4 p) {

  int rank = GetRank(i,j,k, t); 

  if (rank ==-1) return(false); 
  else {
    seedlists[rank].push_back(p); 
    return(true); 
  }

}
//--------------------------------------------------------------------------
//
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
//
bool Lattice4D::InsertSeed(int i, VECTOR4 p) {

  if (i>=npart) return(false); 
  else {
    seedlists[i].push_back(p); 
    return(true); 
  }

}
//--------------------------------------------------------------------------
//
bool Lattice4D::InsertSeed(int from_rank, int to_rank, VECTOR4 p) {

  if (from_rank >=npart || to_rank>=npart) return(false); 
  else {
    flowMatrix[from_rank*npart+to_rank]++; 
    seedlists[to_rank].push_back(p); 
  }
  return(true); 

}
//--------------------------------------------------------------------------
//
void Lattice4D::ResetFlowMatrix() {

  if (flowMatrix !=NULL) 
    memset(flowMatrix, '\0', npart*npart*sizeof(int)); 

}
//--------------------------------------------------------------------------
//
// return the number of partitions that are assigned to processor 'proc' 
//
int Lattice4D::GetNumPartitions(int proc) {

  int n = 0; 
  int i;

  for (i = 0; i < npart; i++) {
    if (part->parts[i].Proc == proc)
      n++;
  }

  return n;

}
//---------------------------------------------------------------------------
//
// query the partitions that are assigned to processor 'proc' 
// p_list must be allocated large enough prior to calling
//
void Lattice4D::GetPartitions(int proc, int*p_list) {

  int n = 0; 
  int i;

  for (i = 0; i < npart; i++) {
    if (part->parts[i].Proc == proc)
      p_list[n++] = i;
  }

  return;

}
//---------------------------------------------------------------------------
//
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
// gets ranks of all neighbors
// ranks are the global partition numbers for all neighbors
// the neighbor of an edge gets rank -1
//
void Lattice4D::GetNeighborRanks(int myrank) {

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
// gets local subvolume bounds
//
// block: local block number (0-nblocks)
// min_s, max_s: (output) spatial min and max bounds
// min_t, max_t: (output) temporal min and max bounds
//
void Lattice4D::GetVB(int block, float *min_s, float *max_s, 
		      int *min_t, int *max_t) {

  min_s[0] = vb_list[block_ranks[block]].xmin;
  min_s[1] = vb_list[block_ranks[block]].ymin;
  min_s[2] = vb_list[block_ranks[block]].zmin;
  max_s[0] = vb_list[block_ranks[block]].xmax;
  max_s[1] = vb_list[block_ranks[block]].ymax;
  max_s[2] = vb_list[block_ranks[block]].zmax;
  *min_t = vb_list[block_ranks[block]].tmin;
  *max_t = vb_list[block_ranks[block]].tmax;

}
//---------------------------------------------------------------------------
//
// gets global subvolume bounds
//
// part: global partition number
// min_s, max_s: (output) spatial min and max bounds
// min_t, max_t: (output) temporal min and max bounds
//
void Lattice4D::GetGlobalVB(int part, float *min_s, float *max_s, 
		      int *min_t, int *max_t) {

  min_s[0] = vb_list[part].xmin;
  min_s[1] = vb_list[part].ymin;
  min_s[2] = vb_list[part].zmin;
  max_s[0] = vb_list[part].xmax;
  max_s[1] = vb_list[part].ymax;
  max_s[2] = vb_list[part].zmax;
  *min_t = vb_list[part].tmin;
  *max_t = vb_list[part].tmax;

}
//---------------------------------------------------------------------------
//
// wrappers around blockwise partition methods
//
// block: local block number (0 - nblocks)
// iter: round number
// ls: seeds list
//
//
// sets the request status
//
void Lattice4D::SetReq(int block) {
  part->SetReq(block_ranks[block]);
}
//
// clears the request status
//
void Lattice4D::ClearReq(int block) {
  part->ClearReq(block_ranks[block]);
}
//
// gets the request status
//
int Lattice4D::GetReq(int block) {
  part->GetReq(block_ranks[block]); 
}
//
// sets the load status
//
void Lattice4D::SetLoad(int block) {
  part->SetLoad(block_ranks[block]); 
}
//
// clears the load status
//
void Lattice4D::ClearLoad(int block) { 
  part->ClearLoad(block_ranks[block]);
}
//
// gets the load status
//
int Lattice4D::GetLoad(int block) {
  part->GetLoad(block_ranks[block]); 
}
//
// sets the computed status
//
void Lattice4D::SetComp(int block, int iter) {
  part->SetComp(block_ranks[block], iter); 
}
//
// clears the computed status
//
void Lattice4D::ClearComp(int block) {
  part->ClearComp(block_ranks[block]); 
}
//
// gets the computed status
//
int Lattice4D::GetComp(int block, int iter) {
  part->GetComp(block_ranks[block], iter); 
}
//
// posts a point for sending
//
void Lattice4D::PostPoint(int block, VECTOR4 p) {
  int neighbor = GetNeighbor(block_ranks[block], p[0], p[1], p[2], p[3]);
  if (neighbor >= 0)
    part->PostPoint(block_ranks[block], p, neighbor); 
}
//
// prints the posted points
//
void Lattice4D::PrintPost(int block) {
  part->PrintPost(block_ranks[block]); 
}
//
// prints the received points
//
void Lattice4D::PrintRecv(int block) { 
  part->PrintRecv(block_ranks[block]); 
}
//
// copies the received points to a list
//
void Lattice4D::GetRecvPts(int block, VECTOR4 *ls) { 
  part->GetRecvPts(block_ranks[block], ls); 
}
//
// sends points to all neighbors
//
void Lattice4D::SendNeighbors(int block) { 
  GetNeighborRanks(block_ranks[block]);
  part->SendNeighbors(block_ranks[block], neighbor_ranks);
}
//
// receives points from all neighbors
//
int Lattice4D::ReceiveNeighbors(int block) {
  GetNeighborRanks(block_ranks[block]);
  part->ReceiveNeighbors(block_ranks[block], neighbor_ranks);
}
//---------------------------------------------------------------------------
