
#include "Lattice4D.h"

//---------------------------------------------------------------------------
//
// constructs and initializes the time-varying lattice class
// using default round robin partitioning
//
// xlen, ylen, zlen, tlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of parition in time 
// nid: number of processes, threads, owners
// default = 1 (can omit if single process sequential program)
// myid: rank, process number, thread number, identification of the owner
// default = 0 (can omit if single process sequential program)
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, 
		     int nsp, int ntp, int nid, int myid) {

  int i, j;
  volume_bounds_type *vbs; 

  // extents
  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  ldim = tlen;
  min_extent[0] = min_extent[1] = min_extent[2] = min_extent[3] = 0.0;
  max_extent[0] = xlen - 1;
  max_extent[1] = ylen - 1;
  max_extent[2] = zlen - 1;
  max_extent[3] = tlen - 1;

  // spatial domain partitioning
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

  part = new Partition(nsp * ntp, nid);
  delete [] vbs; 
  myproc = myid;
  nproc = nid;

  if (myproc >= 0) {

    // assign partitions
    RoundRobin_proc();
    nb = GetNumPartitions(myproc);

    // allocate table of neighbor ranks for one neighbor initially
    assert((neighbor_ranks = (int **)malloc(nb * sizeof(int *))) != NULL);
    for (i = 0; i < nb; i++)
      assert((neighbor_ranks[i] = (int *)malloc(sizeof(int))) != NULL);

    // ranks of my blocks
    assert((block_ranks = (int *)malloc(nb * sizeof(int))) != NULL);
    GetPartitions(myproc, block_ranks);

    // learn who my neighbors are
    for(i = 0; i < nb; i++)
      GetNeighborRanks(i);

  }

  // average number of neighbors for each block
  int tot_neighbors = 0;
  for (i = 0; i < nb; i++)
    tot_neighbors += part->parts[block_ranks[i]].NumNeighbors;
  avg_neigh = tot_neighbors / nb;

  // more stats
  tot_pts_send = 0;

}
//---------------------------------------------------------------------------
//
// constructs and initializes the time-varying lattice class
// using explicit partitioning
//
// currently only works for one time-step
//
// part_file: partitioning information file
// xlen, ylen, zlen, tlen: total size of the data
// ghost: size of ghost layer per side
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of parition in time 
// nid: number of processes, threads, owners
// default = 1 (can omit if single process sequential program)
// myid: rank, process number, thread number, identification of the owner
// default = 0 (can omit if single process sequential program)
//
Lattice4D::Lattice4D(char *part_file, int xlen, int ylen, int zlen, int tlen,
		     int ghost, int nsp, int ntp, int nid, int myid) {

  int nprocs, nblocks; // number of procs, number of blocks
  int block_size[3]; // block size
  int *block_procs; // process assigned to each block
  float *block_extents; // extent of each block
  FILE *fp;
  int idx; // block index
  int i, t;

  // open the partition file
  assert((fp = fopen(part_file,"rb")) != NULL);

  // global extents
  fread(&xdim, sizeof(int), 1, fp);
#ifdef BYTE_SWAP
  swap4((char *)&xdim);
#endif
  assert(xdim == xlen); // size in partition file must match run script
  fread(&ydim, sizeof(int), 1, fp);
#ifdef BYTE_SWAP
  swap4((char *)&ydim);
#endif
  assert(ydim == ylen);
  fread(&zdim, sizeof(int), 1, fp);
#ifdef BYTE_SWAP
  swap4((char *)&zdim);
#endif
  assert(zdim == zlen);
  tdim = ldim = tlen;
  min_extent[0] = min_extent[1] = min_extent[2] = min_extent[3] = 0.0;
  max_extent[0] = xdim - 1;
  max_extent[1] = ydim - 1;
  max_extent[2] = zdim - 1;
  max_extent[3] = tdim - 1;

  // number of procs, blocks, block size
  fread(&nprocs, sizeof(int), 1, fp);
#ifdef BYTE_SWAP
  swap4((char *)&nprocs);
#endif
  fread(&nblocks, sizeof(int), 1, fp);
#ifdef BYTE_SWAP
  swap4((char *)&nprocs);
#endif
  fread(block_size, sizeof(int), 3, fp);
#ifdef BYTE_SWAP
  for (i = 0; i < 3; i++)
    swap4((char *)&(block_size[i]));
#endif
  assert(nprocs == nid); // # procs in partition file must match run script
  assert(nblocks == nsp); // # blocks in partition file must match run script

  // process assignment and block extents
  assert((block_procs = new int[nblocks]) != NULL);
  assert((block_extents = new float[nblocks * 6]) != NULL);
  fread(block_procs, sizeof(int), nblocks, fp);
#ifdef BYTE_SWAP
  for (i = 0; i < nblocks; i++)
    swap4((char *)&(block_procs[i]));
#endif
  fread(block_extents, sizeof(float), nblocks * 6, fp);
#ifdef BYTE_SWAP
  for (i = 0; i < nblocks * 6; i++)
    swap4((char *)&(block_extents[i]));
#endif

  fclose(fp);

  // add ghost cells
  for (i = 0; i < nblocks; i++) {
    if (block_extents[i * 6] > 0)
      block_extents[i * 6] -= ghost;
    if (block_extents[i * 6 + 1] > 0)
      block_extents[i * 6 + 1] -= ghost;
    if (block_extents[i * 6 + 2] > 0)
      block_extents[i * 6 + 2] -= ghost;
    if (block_extents[i * 6 + 3] < xdim)
      block_extents[i * 6 + 3] += ghost;
    if (block_extents[i * 6 + 4] < ydim)
      block_extents[i * 6 + 4] += ghost;
    if (block_extents[i * 6 + 5] < zdim)
      block_extents[i * 6 + 5] += ghost;
  }

  // print contents
  if (myid == 0) {
    printf("xdim ydim zdim nprocs nblocks = %d %d %d %d %d\n",
	   xdim, ydim, zdim, nprocs, nblocks);
    printf("block size = %d %d %d\n", 
	   block_size[0], block_size[1], block_size[2]);
//     for(i = 0; i < nblocks; i++)
//       printf("block_procs[%d] = %d\n", i, block_procs[i]);
//     for(i = 0; i < nblocks; i++)
//       printf("block %d min = [%.0f %.0f %.0f] max = [%.0f %.0f %.0f]\n",
// 	     i, block_extents[i * 6], block_extents[i * 6 + 1], 
//           block_extents[i * 6 + 2],
// 	     block_extents[i * 6 + 3], 
//           block_extents[i * 6 + 4], block_extents[i * 6 + 5]);
  }

  // spatial domain partitioning
  npart = nblocks * ntp; 
  vb_list = new volume_bounds_type[npart]; 
  VolumeBounds(block_extents, nblocks, block_size, vb_list, idim, jdim, kdim); 

  // temporal domain partitioning
  idx = 0;
  for (t = 0; t < tdim; t++) {
    vb_list[idx].tmin = tdim / ntp * t - 1; 
    if (vb_list[idx].tmin < 0)
      vb_list[idx].tmin = 0; 
    vb_list[idx].tmax = tdim / ntp * (t + 1) + 1; 
    if (vb_list[idx].tmax > tdim - 1)
      vb_list[idx].tmax = tdim - 1; 
    idx++;
  }

  part = new Partition(nblocks * ntp, nid);
  myproc = myid;
  nproc = nid;

  if (myproc >= 0) {

    // assign partitions
    Explicit_proc(block_procs);
    nb = GetNumPartitions(myproc);

    // allocate table of neighbor ranks for one neighbor initially
    assert((neighbor_ranks = (int **)malloc(nb * sizeof(int *))) != NULL);
    for (i = 0; i < nb; i++)
      assert((neighbor_ranks[i] = (int *)malloc(sizeof(int))) != NULL);

    // ranks of my blocks
    assert((block_ranks = (int *)malloc(nb * sizeof(int))) != NULL);
    GetPartitions(myproc, block_ranks);

    // learn who my neighbors are
    for(i = 0; i < nb; i++)
      GetNeighborRanks(i);

  }

  // average number of neighbors for each block
  int tot_neighbors = 0;
  for (i = 0; i < nb; i++)
    tot_neighbors += part->parts[block_ranks[i]].NumNeighbors;
  avg_neigh = tot_neighbors / nb;

  // more stats
  tot_pts_send = 0;

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
//   delete block_ranks;

}
//---------------------------------------------------------------------------
//
// GetExtents
//
// overall min, max in x,y,z,t of the entire data
//
void Lattice4D::GetExtents(float *min, float *max) {

  min[0] = min_extent[0];
  min[1] = min_extent[1];
  min[2] = min_extent[2];
  min[3] = min_extent[3];
  max[0] = max_extent[0];
  max[1] = max_extent[1];
  max[2] = max_extent[2];
  max[3] = max_extent[3];

}
//----------------------------------------------------------------------------
//
// VolumeBounds
//
// populate the volume bounds list and determine size in i, j, k
//
// block_extents: explicit extents of each block
// nblocks: total number of blocks
// block_size: size of block, eg 16x16x16
// vb_list (output): volume bounds list
// idim, jdim, kdim (output): index sizes
//
void Lattice4D::VolumeBounds(float *block_extents, int nblocks, 
			     int *block_size, volume_bounds_type *vb_list, 
			     int &idim, int &jdim, int &kdim){
  int i;

  for (i = 0; i < nblocks; i++) {
    vb_list[i].xmin = block_extents[i * 6];
    vb_list[i].ymin = block_extents[i * 6 + 1];
    vb_list[i].zmin = block_extents[i * 6 + 2];
    vb_list[i].xmax = block_extents[i * 6 + 3];
    vb_list[i].ymax = block_extents[i * 6 + 4];
    vb_list[i].zmax = block_extents[i * 6 + 5];
  }

  idim = ceil((float)xdim / (float)block_size[0]);
  jdim = ceil((float)ydim / (float)block_size[1]);
  kdim = ceil((float)zdim / (float)block_size[2]);

}
//----------------------------------------------------------------------------
//
// assign the partitions to the processors explicitely
//
void Lattice4D::Explicit_proc(int *block_procs) {

  int proc; // process number
  int n; // index into a row of the process table, local block number
  int i; // partition rank

  for (i = 0; i < nproc; i++)
    part->proc_nparts[i] = 0;

  for (i = 0; i < npart; i++) {

    proc = block_procs[i];
    n = part->proc_nparts[proc];
    part->parts[i].Proc = proc; 
    part->proc_parts[proc][n] = i;
    part->proc_nparts[proc]++;

  }

}
//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void Lattice4D::RoundRobin_proc() {

  int proc; // process number
  int n; // index into a row of the process table, local block number
  int i; // partition rank

  for (i = 0; i < npart; i++) {

    proc = i % nproc;
    n = i / nproc;

    // assign the process number
    part->parts[i].Proc = proc; 

    // add entry into the process table
    part->proc_parts[proc][n] = i;
    part->proc_nparts[proc] = n + 1;

  }

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

  int idx = t* idim * jdim * kdim + k * idim * jdim + j * idim + i; 

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

  int neighbor_rank;
  int si, sj, sk;

  neighbor_rank = CheckNeighbor(myrank, x, y, z, t);
  
  GetIndices(neighbor_rank, ei, ej, ek, et); 
  
  return (neighbor_rank); 

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
// return the number of partitions that are assigned to my proc
//
int Lattice4D::GetMyNumPartitions() {

  int n = 0; 
  int i;

  for (i = 0; i < npart; i++) {
    if (part->parts[i].Proc == myproc)
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
// returns neighbor number of neighbor containing point
//
// returns -1 if point is not in one of the neighbors, ie, is out of bounds
// returns my neighbor number if point stayed in my block (valid)
//
int Lattice4D::GetNeighbor(int block, float x, float y, float z, float t) {

  int n; // neighbor number
  int r; // rank of partition containing the point

  r = GetRank(x, y, z, t);

  // for all neighbors
  for (n = 0; n < part->parts[block_ranks[block]].NumNeighbors; n++) {
    if (r == neighbor_ranks[block][n])
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
// block: my local block number
//
void Lattice4D::GetNeighborRanks(int block) {

  int my_i, my_j, my_k, my_l; // my lattice coords
  int i, j, k, l; // lattice coords offset
  int nr; // my neighbor's rank (global partition) number
  int myrank = block_ranks[block];

  GetIndices(myrank, my_i, my_j, my_k, my_l);

  for (i = -1; i <= 1; i++) {
    for (j = -1; j <= 1; j++) {
      for (k = -1; k <= 1; k++) {
	for (l = -1; l <= 1; l++) {
	  nr = GetRank(my_i + i, my_j + j, my_k + k, my_l + l);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	}
      }
    }
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
// adds a neighbor to the end of the neighbor_ranks table and sets its value
// updates the local block neighbors as well as the global partition
// grows local block neighbors if necessary
//
// myblock: my local block number
// neighrank: global partition number of the neighbor
//
void Lattice4D::AddNeighbor(int myblock, int neighrank) {

  int nn; // number of neighbors allocated
  int myrank = block_ranks[myblock];
  int num_neighbors = part->GetNumNeighbors(myrank);

  if ((nn = part->GetAllocNeighbors(myrank)) < num_neighbors + 1)
    assert((neighbor_ranks[myblock] = (int *)realloc(neighbor_ranks[myblock],
				    nn * 2 * sizeof(int))) != NULL);

  neighbor_ranks[myblock][num_neighbors] = neighrank;

  part->AddNeighbor(myrank, myblock, neighrank);

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
  return part->GetReq(block_ranks[block]); 
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
  return part->GetLoad(block_ranks[block]); 
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
  return part->GetComp(block_ranks[block], iter); 
}
//
// posts a point for sending
//
void Lattice4D::PostPoint(int block, VECTOR4 p) {
  int neighbor = GetNeighbor(block, p[0], p[1], p[2], p[3]);
  // only post points that remain inside the overall domain boundary
  if (neighbor >= 0) {
    part->PostPoint(block_ranks[block], p, neighbor); 
    // for performance stats, only count points that leave my proc
    if (part->GetProc(neighbor_ranks[block][neighbor]) != myproc)
      tot_pts_send++;
  }
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

#ifdef _MPI
//
// exchanges points with all neighbors (MPI version)
// returns total number of points received by this process
//
int Lattice4D::ExchangeNeighbors(VECTOR4 **seeds, int *size_seeds,
				 int *num_seeds) { 
  int n;
  comm_time = MPI_Wtime();
  n = part->ExchangeNeighbors(neighbor_ranks, seeds, size_seeds, num_seeds);
  comm_time = MPI_Wtime() - comm_time;
  return n;
}
#endif
//
// exchanges points with all neighbors (serial version)
// returns total number of points received by this process
//
int Lattice4D::SerExchangeNeighbors(VECTOR4 **seeds, int *size_seeds, 
				 int *num_seeds) { 
  int n;
  n = part->SerExchangeNeighbors(neighbor_ranks, seeds, size_seeds, num_seeds);
  return n;
}

//---------------------------------------------------------------------------
