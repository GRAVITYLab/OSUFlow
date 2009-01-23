
#include "Lattice4D.h"
//---------------------------------------------------------------------------
//
// Lattice 4D 
//
// constructs and initializes the time-varying lattice class
//
// xlen, ylen, zlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total number of spatial partitions in the lattice
// ntp: total number of parition in time 
//
// Tom Peterka, 12/8/08
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp) {

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
// GetPartitions
//
// query the paritions that are assigned to processor 'proc' 
//
// p_list must be allocated large enough prior to calling
//
// Tom Peterka, 12/5/08
//
void Lattice4D::GetPartitions(int proc, int*p_list, int& num) {

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
int Lattice4D::GetRank(int i, int j, int k, int t) {

  if (i < 0 || i >= idim) 
    return(-1); 
  else if (j < 0 || j >= jdim) 
    return(-1); 
  else if (k < 0 || k >= kdim) 
    return(-1); 
  else if (t < 0 || t >= tdim) 
    return(-1); 

  int idx = t * idim*jdim*kdim + k * idim * jdim + j * idim + i; 

  return (idx); 

}
//----------------------------------------------------------------------------
//
// GetNeighborRanks
//
// gets ranks of 80 neighbors
// the neighbor of an edge gets rank -1
//
// Tom Peterka, 11/26/08
//
void Lattice4D::GetNeighborRanks(int myrank, int *neighbor_ranks) {

  int i, j, k, t; // my lattice coords

  GetIndices(myrank, i, j, k, t);

  int idx = 0; 
  for (int nt=-1; nt<=1; nt++) 
    for (int nk=-1; nk<=1; nk++) 
      for (int nj=-1; nj<=1; nj++) 
	for (int ni=-1; ni<=1; ni++)  {
	  if (nt==0 && nk==0 && nj==0 && ni==0) continue; //self 
	  else if (nt+t<0 || nt+t>tdim-1) neighbor_ranks[idx] = -1; 
	  else if (nk+k<0 || nk+k>kdim-1) neighbor_ranks[idx] = -1; 
	  else if (nj+j<0 || nj+j>jdim-1) neighbor_ranks[idx] = -1; 
	  else if (ni+i<0 || ni+i>idim-1) neighbor_ranks[idx] = -1; 
	  else neighbor_ranks[idx] = GetRank(ni+i, nj+j, nk+k, nt+t); 
	  idx++; 
	}
}
//---------------------------------------------------------------------------
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
//check if the point (x,y,z) is in the lattice element [i,j,k,l]
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
//---------------------------------------------------------------------------
//
// GetNeighbor
//
// returns neighbor number (0-5) of x,y,z point and sets i,j,k of neighbor
// myrank: global partition number
//
// Tom Peterka, 11/26/08
//
int Lattice4D::GetNeighbor(int myrank, float x, float y, float z, float t, 
			   int &ei, int &ej, int &ek, int &et) {

  int neighbor;
  int si, sj, sk;

  neighbor = CheckNeighbor(myrank, x, y, z, t);
  
  GetIndices(neighbor, ei, ej, ek, et); 
  
  return (neighbor); 

}
//----------------------------------------------------------------------------

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

void Lattice4D::InsertSeed(int i, int j, int k, int t, VECTOR4 p) {

  int rank = GetRank(i,j,k, t); 

  seedlists[rank].push_back(p); 

}

void Lattice4D::InsertSeed(int i, VECTOR4 p) {

  seedlists[i].push_back(p); 

}


