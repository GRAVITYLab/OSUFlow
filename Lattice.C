
#include "Lattice.h"

Lattice::Lattice(int xlen, int ylen, int zlen, int ghost, int np) {

  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  vb_list = calc_subvolume(xlen, ylen, zlen, 
			   ghost, np, 
			   idim, jdim, kdim); 
  npart = np; 
  proc_list = new int[np]; 
  for (int i=0; i<np; i++) 
    proc_list[i] = i; // this is a default assignment. the processor which is responsible for 
                      // the i-th vb is just i 
                      // remember proc_list[i] stores the rank of proc that is responsible for i-th vb 
}

// assign the partitions to the processors in a round-robin manner 
void Lattice::RoundRobin_proc(int nproc) {

  for (int i=0; i<npart; i++) 
    proc_list[i] = i % nproc; 



}

// look up the lattice[i,j,k] element to get the processor for the subdomain number (rank)
int Lattice::GetProc(int i, int j, int k) {

  if (i<0 || i >=idim) return(-1); 
  else if (j<0 || j >=jdim) return(-1); 
  else if (k<0 || k >=kdim) return(-1); 

  int idx = k*idim*jdim+j*idim+i; 
  return (proc_list[idx]); 
}

int Lattice::GetProc(int rank) {

  return(proc_list[rank]); 
}

// query the paritions that are assigned to processor 'proc' 
void Lattice::GetPartitions(int proc, int**p_list, int& num) {

  int cnt = 0; 
  for (int i=0; i<npart; i++) 
    if (proc_list[i] == proc) cnt ++; 
  num = cnt; 
  (*p_list) = new int[cnt]; 
  cnt = 0; 
  for (int i=0; i<npart; i++) 
    if (proc_list[i] == proc)  {
      (*p_list)[cnt++] = i; 
    }
}


// look up the lattice[i,j,k] element to get the subdomain number (rank)
int Lattice::GetRank(int i, int j, int k) {

  if (i<0 || i >=idim) return(-1); 
  else if (j<0 || j >=jdim) return(-1); 
  else if (k<0 || k >=kdim) return(-1); 

  int idx = k*idim*jdim+j*idim+i; 
  return (idx); 
}

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

// Find the subdomain that contains the physical location (x,y,z) 
int Lattice::GetRank(float x, float y, float z) {
  
  for (int i=0; i<npart; i++) {
    if (vb_list[i].xmin>x ||vb_list[i].xmax<x) continue; 
    else if (vb_list[i].ymin>y ||vb_list[i].ymax<y) continue; 
    else if (vb_list[i].zmin>z ||vb_list[i].zmax<z) continue; 
    return(i); 
  }
  return(-1); 
}

// find the indices of the lattice element that has its subdomain number = rank 
int Lattice::GetIndices(int rank, int &iidx, int &jidx, int &kidx) {
  
  if (rank <0 || rank >= npart) return(-1); 

  kidx = rank/(idim*jdim) ; 
  jidx = (rank % (idim*jdim))/idim; 
  iidx = rank % idim; 

  return(1); 
}
// find the indices of the lattice element that contanis (x,y,z) 
int Lattice::GetIndices(float x, float y, float z, int &iidx, int &jidx, int &kidx) {

  int rank = GetRank(x,y,z); 
  if (rank!=-1) {
    GetIndices(rank, iidx, jidx, kidx); 
    return(rank); 
  }
  else return(-1); 

}

// look up the volume bounds of lattice[i,j,k] 
int Lattice::GetBounds(int i, int j, int k, volume_bounds_type& vb)  {

  if (i<0 || i >=idim) return(-1); 
  else if (j<0 || j >=jdim) return(-1); 
  else if (k<0 || k >=kdim) return(-1); 
  
  int idx = k*idim*jdim+j*idim+i; 
  vb = vb_list[idx]; 
  return(1); 
}

// look up the volume bounds of the subdomain 'rank' 
int Lattice::GetBounds(int rank, volume_bounds_type &vb) {

  if (rank <0 || rank >= npart) return(-1); 
  vb = vb_list[rank]; 
  return(1); 

}

//check if the point (x,y,z) is in the lattice element [i,j,k]
bool Lattice::isIn(float x, float y, float z, int i, int j, int k) 
{
  if (i<0 || i>idim-1 || j<0 || j>jdim-1 || k<0 || k>kdim-1) return(false); 
  int idx = k *idim*jdim + j*idim + i; 
  volume_bounds_type vb = vb_list[idx]; 
  if (vb.xmin>x || vb.xmax<x) return (false); 
  else if (vb.ymin>y || vb.ymax<y) return (false); 
  else if (vb.zmin>z || vb.zmax<z) return (false); 
  return(true); 
}

// returns neighbor number (0-5) of x,y,z point
int Lattice::CheckNeighbor(int myrank, float x, float y, float z)
{
  int i,j,k; 
  GetIndices(myrank, i, j, k); 
  if (isIn(x,y,z,i-1,j,k)==true) return(0);  //-X
  if (isIn(x,y,z,i+1,j,k)==true) return(1);  //+X
  if (isIn(x,y,z,i,j-1,k)==true) return(2);  //-Y
  if (isIn(x,y,z,i,j+1,k)==true) return(3);  //+Y
  if (isIn(x,y,z,i,j,k-1)==true) return(4);  //-Z
  if (isIn(x,y,z,i,j,k+1)==true) return(5);  //+Z
  return(-1); 
}

// returns neighbor number (0-5) of x,y,z point and sets i,j,k of neighbor
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

void Lattice::InitSeedLists() 
{
  seedlists = new list<VECTOR3>[npart]; 
  for (int i=0; i<npart; i++)
    seedlists[i].clear(); 
}

void Lattice::ResetSeedLists() 
{
  for (int i=0; i<npart; i++)
    seedlists[i].clear(); 
}

void Lattice::InsertSeed(int i, int j, int k, VECTOR3 p) 
{
  int rank = GetRank(i,j,k); 
  //  printf("Lattice::InsertSeed to %d\n", rank); 
  seedlists[rank].push_back(p); 
}

