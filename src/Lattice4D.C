//---------------------------------------------------------------------------
//
// 4D lattice
//
// Han-Wei Shen
// The Ohio State University
// Columbus, OH
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

#include "Lattice4D.h"

//---------------------------------------------------------------------------
//
// constructs and initializes the time-varying lattice class
// using default round robin partitioning
//
// xlen, ylen, zlen, tlen: total size of the data
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of partitions in time 
// nid: number of processes, threads, owners
// default = 1 (can omit if single process sequential program)
// myid: rank, process number, thread number, identification of the owner
// default = 0 (can omit if single process sequential program)
// track_ids: should id tracking be turned on. defaults to false
// ghost: size of ghost layer per side
// default = 0
//
Lattice4D::Lattice4D(int xlen, int ylen, int zlen, int tlen, 
		     int nsp, int *ntp, int nid, int myid, 
		     bool track_ids, int ghost)
    : flowMatrix(NULL), alloc_neighbors(NULL), block_ranks(NULL) {

  int lat_dim[3]; // number of blocks in each direction
  int data_dim[4]; // xlen, ylen, zlen
  int i, j;

  myproc = myid;
  tdim = *ntp; 
  npart = nsp * tdim; 

  // extents
  xdim = xlen; 
  ydim = ylen; 
  zdim = zlen; 
  ldim = tlen;
  min_extent[0] = 0; min_extent[1] = 0; 
  min_extent[2] = 0; min_extent[3] = 0;
  max_extent[0] = xlen - 1; max_extent[1] = ylen - 1; 
  max_extent[2] = zlen - 1; max_extent[3] = tlen - 1;

  // compute partitions
  data_dim[0] = xlen; data_dim[1] = ylen; 
  data_dim[2] = zlen; data_dim[3] = tlen;
  ComputePartition(data_dim, ghost, nsp, *ntp, lat_dim);
  idim = lat_dim[0]; jdim = lat_dim[1]; kdim = lat_dim[2];

  // for time-varying flow, number time blocks must be < number time steps
  if (tlen > 1 && *ntp >= tlen)
    *ntp = tlen - 1;

  part = new Partition(npart, nid, myid, track_ids);
  nproc = nid;

  // partition
  if (myproc >= 0) {

    // assign partitions
    RoundRobin_proc();
    nb = part->nb;
    alloc_blocks = nb;

    // init my blocks
    if (alloc_blocks > 0)
      assert((block_ranks = (int *)malloc(nb * sizeof(int))) != NULL);
    GetPartitions(myproc, block_ranks);

    // init my neighbors
    if (alloc_blocks > 0) {
      assert((alloc_neighbors = 
	      (int *)malloc(alloc_blocks * sizeof(int))) != NULL);
      for (j = 0; j < alloc_blocks; j++)
	alloc_neighbors[j] = 0;
      assert((neighbor_ranks = (int **)malloc(nb * sizeof(int *))) != NULL);
      assert((neighbor_procs = (int **)malloc(nb * sizeof(int *))) != NULL);
    }

    // learn who my neighbors are
    for(i = 0; i < nb; i++)
      GetNeighborRanks(i); // gets neighbor_procs too

    // make the parts data structure strictly local from now on
    // by clearing the process number for blocks not belonging to me
    for (i = 0; i < npart; i++) {
      if (part->parts[i].Proc != myproc)
	part->parts[i].Proc = -1;
    }

  }

#ifdef ZOLTAN

  InitRepartition4D(this, MPI_COMM_WORLD);
//   ChangePartition(&nb, &block_ranks, &neighbor_ranks, &neighbor_procs, part,
// 		  NULL, NULL, NULL, &avg_neigh, &alloc_blocks, 
// 		  &alloc_neighbors, MPI_COMM_WORLD, &AddNeighbor);

#endif

  // performance stats
  int tot_neighbors = 0;
  for (i = 0; i < nb; i++)
    tot_neighbors += part->parts[block_ranks[i]].NumNeighbors;
  avg_neigh = (nb > 0 ? tot_neighbors / nb : 0);
  tot_pts_send = 0;

}
//---------------------------------------------------------------------------
//
// constructs and initializes the time-varying lattice class
// using explicit partitioning
//
// part_file: partitioning information file
// xlen, ylen, zlen, tlen: total size of the data
// nsp: total (global) number of spatial partitions in the lattice
// ntp: total (global) number of partitions in time 
// nid: number of processes, threads, owners
// default = 1 (can omit if single process sequential program)
// myid: rank, process number, thread number, identification of the owner
// default = 0 (can omit if single process sequential program)
// track_ids: should id tracking be turned on, defaults to false
// ghost: size of ghost layer per side
// default = 0
//
Lattice4D::Lattice4D(char *part_file, int xlen, int ylen, int zlen, int tlen,
		     int nsp, int *ntp, int nid, int myid, 
		     bool track_ids, int ghost)
    : flowMatrix(NULL), alloc_neighbors(NULL), block_ranks(NULL) {

  int nprocs, nblocks;  // number of procs, number of blocks
  int ntimepart;        // number of time partitions
  int block_size[4];    // block size
  int *block_procs;     // process assigned to each block
  float *block_extents; // extent of each block
  FILE *fp;
  int i;

  // open the partition file
  assert((fp = fopen(part_file,"rb")) != NULL);

  // global extents
  assert(fread(&xdim, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&xdim);
#endif
  assert(xdim == xlen); // size in partition file must match run script

  assert(fread(&ydim, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&ydim);
#endif
  assert(ydim == ylen);

  assert(fread(&zdim, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&zdim);
#endif
  assert(zdim == zlen);

  assert(fread(&ldim, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&ldim);
#endif
  assert(ldim == tlen);

  // for time-varying flow, number time blocks must be < number time steps
  if (tlen > 1 && *ntp >= tlen)
    *ntp = tlen - 1;

  tdim = *ntp;
  ldim = tlen;
  min_extent[0] = 0; min_extent[1] = 0; 
  min_extent[2] = 0; min_extent[3] = 0;
  max_extent[0] = xlen - 1; max_extent[1] = ylen - 1; 
  max_extent[2] = zlen - 1; max_extent[3] = tlen - 1;

  // number of procs, blocks, block size
  assert(fread(&nprocs, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&nprocs);
#endif

  assert(fread(&nblocks, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&nblocks);
#endif

  assert(fread(&ntimepart, sizeof(int), 1, fp) > 0);
#ifdef BYTE_SWAP
  swap4((char *)&ntimepart);
#endif

  assert(fread(block_size, sizeof(int), 4, fp) > 0);
#ifdef BYTE_SWAP
  for (i = 0; i < 4; i++)
    swap4((char *)&(block_size[i]));
#endif

  assert(nprocs == nid);     // # procs in partition file must match run script
  assert(nblocks == nsp);    // # blocks in partition file must match run script
  assert(ntimepart == *ntp); // # of time partition file must match run script

  npart = nblocks * tdim; 

  // process assignment and block extents
  assert((block_procs = new int[npart]) != NULL);
  assert((block_extents = new float[npart * 8]) != NULL);
  assert(fread(block_procs, sizeof(int), npart, fp) > 0);
#ifdef BYTE_SWAP
  for (i = 0; i < npart; i++)
    swap4((char *)&(block_procs[i]));
#endif
  assert(fread(block_extents, sizeof(float), npart * 8, fp) > 0);
#ifdef BYTE_SWAP
  for (i = 0; i < npart * 8; i++)
    swap4((char *)&(block_extents[i]));
#endif

  fclose(fp);

#ifdef DEBUG
  // print contents
  if (myid == 0) {
    fprintf(stderr, "xdim ydim zdim ldim nprocs nblocks ntimepart = %d %d %d %d %d %d %d\n",
	   xdim, ydim, zdim, ldim, nprocs, nblocks, ntimepart);
    fprintf(stderr, "block size = %d %d %d %d\n", 
	    block_size[0], block_size[1], block_size[2], block_size[3]);
    for(i = 0; i < npart; i++)
      fprintf(stderr, "block_procs[%d] = %d\n", i, block_procs[i]);
    for(i = 0; i < npart; i++)
      fprintf(stderr, 
	   "block %d min = [%.0f %.0f %.0f %.0f] max = [%.0f %.0f %.0f %.0f]\n",
	     i, block_extents[i * 8], block_extents[i * 8 + 1], 
	        block_extents[i * 8 + 2], block_extents[i * 8 + 3], 
		block_extents[i * 8 + 4], block_extents[i * 8 + 5], 
		block_extents[i * 8 + 6], block_extents[i * 8 + 7]);
  }
#endif

  // volume and time bounds
  vb_list = new volume_bounds_t[npart];  // volume bounds w/ ghost cells
  vbr_list = new volume_bounds_t[npart]; // volume bounds w/o ghost cells
  for (i = 0; i < npart; i++) {
    vbr_list[i].xmin = vb_list[i].xmin = block_extents[i * 8];
    vbr_list[i].ymin = vb_list[i].ymin = block_extents[i * 8 + 1];
    vbr_list[i].zmin = vb_list[i].zmin = block_extents[i * 8 + 2];
    vbr_list[i].tmin = vb_list[i].tmin = block_extents[i * 8 + 3];
    vbr_list[i].xmax = vb_list[i].xmax = block_extents[i * 8 + 4];
    vbr_list[i].ymax = vb_list[i].ymax = block_extents[i * 8 + 5];
    vbr_list[i].zmax = vb_list[i].zmax = block_extents[i * 8 + 6];
    vbr_list[i].tmax = vb_list[i].tmax = block_extents[i * 8 + 7];
  }

  idim = ceil((float)xdim / (float)block_size[0]);
  jdim = ceil((float)ydim / (float)block_size[1]);
  kdim = ceil((float)zdim / (float)block_size[2]);

  ApplyGhost(ghost);
  FindTimeBounds();

  part = new Partition(npart, nid, myid, track_ids);

  myproc = myid;
  nproc = nid;

  if (myproc >= 0) {

    // assign partitions
    Explicit_proc(block_procs);
    nb = part->nb;
    alloc_blocks = nb;

    // init my blocks
    if (alloc_blocks > 0)
      assert((block_ranks = (int *)malloc(nb * sizeof(int))) != NULL);
    GetPartitions(myproc, block_ranks);

    // init my neighbors
    if (alloc_blocks > 0) {
      assert((alloc_neighbors = 
	      (int *)malloc(alloc_blocks * sizeof(int))) != NULL);
      for (i = 0; i < alloc_blocks; i++)
	alloc_neighbors[i] = 0;
      assert((neighbor_ranks = (int **)malloc(nb * sizeof(int *))) != NULL);
      assert((neighbor_procs = (int **)malloc(nb * sizeof(int *))) != NULL);
    }

    // learn who my neighbors are
    for(i = 0; i < nb; i++)
      GetNeighborRanks(i); // handles neighbor_procs too

    // make the parts data structure strictly local from now on
    // by clearing the process number for blocks not belonging to me
    for (i = 0; i < npart; i++) {
      if (part->parts[i].Proc != myproc)
	part->parts[i].Proc = -1;
    }

  }

  // performance stats
  int tot_neighbors = 0;
  for (i = 0; i < nb; i++)
    tot_neighbors += part->parts[block_ranks[i]].NumNeighbors;
  avg_neigh = (nb > 0 ? tot_neighbors / nb : 0);
  tot_pts_send = 0;
}
//---------------------------------------------------------------------------
//
// destructor
//
Lattice4D::~Lattice4D()
{

  if (alloc_neighbors != NULL) free(alloc_neighbors);
  if (block_ranks != NULL) free(block_ranks);
  if (flowMatrix!=NULL) delete [] flowMatrix; 
//   if (vb_list!=NULL) delete [] vb_list; 

  delete part;
//   delete block_ranks;

}
//---------------------------------------------------------------------------
//
// manually time the boundaries of each time group by looking at every block
//
void Lattice4D::FindTimeBounds()
{
  // find bounds of time partitions
  tb_list = new time_bounds_t[tdim];
  int ti = 0;
  int i, j;
  for(i = 0; i < npart; i++) {

    int tmin = vb_list[i].tmin;
    int tmax = vb_list[i].tmax;

    // see if these time bounds are already accounted
    for(j=0; j<ti; j++) {
      if(tb_list[j].tmin == tmin && tb_list[j].tmax == tmax)
	break;
    }

    if(j == ti) {
      // add new bounds
      tb_list[ti].tmin = tmin;
      tb_list[ti].tmax = tmax;
      ti++;
    }
  }

  // sort list since blocks could have been in any order
  qsort(tb_list, tdim, sizeof(time_bounds_t), Lattice4D::compare_time_bounds);
}
// used to compare and sort time bounds list
int Lattice4D::compare_time_bounds(const void* a, const void* b)
{
  time_bounds_t* aa = (time_bounds_t*)a;
  time_bounds_t* bb = (time_bounds_t*)b;

  if(aa->tmin == bb->tmin)
    return aa->tmax - bb->tmax;
  else
    return aa->tmin - bb->tmin;
}
//---------------------------------------------------------------------------
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
			     int *block_size, volume_bounds_t *vb_list, 
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

  int i;

  part->nb = 0;
  for (i = 0; i < npart; i++) {
    part->parts[i].Proc = block_procs[i]; 
    if (part->parts[i].Proc == myproc)
      part->nb++;
  }

}
//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void Lattice4D::RoundRobin_proc() {

  int i;

  part->nb = 0;
  for (i = 0; i < npart; i++) {
    part->parts[i].Proc = i % nproc; 
    if (part->parts[i].Proc == myproc)
      part->nb++;
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
// query the partitions that are assigned to processor 'proc' 
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
// Takes into account ghost cells
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
//----------------------------------------------------------------------------
//
// Find the subdomain that contains the physical location (x,y,z) 
// Does not take into account ghost cells
//
int Lattice4D::GetRealRank(float x, float y, float z, float t) {
  
  for (int i = 0; i < npart; i++) {

    if (vbr_list[i].xmin > x || vbr_list[i].xmax < x) 
      continue; 
    else if (vbr_list[i].ymin > y || vbr_list[i].ymax < y) 
      continue; 
    else if (vbr_list[i].zmin > z || vbr_list[i].zmax < z) 
      continue; 
    else if (vbr_list[i].tmin > t || vbr_list[i].tmax < t) 
      continue; 
    return(i); 
 
 }

  return(-1); 

}
//---------------------------------------------------------------------------
//
// Find the subdomain that contains the physical location (x,y,z) 
// my hacked version that does space with ghost cells and time without
//
int Lattice4D::MyGetRank(float x, float y, float z, float t) {
  
  for (int i = 0; i < npart; i++) {

    if (vb_list[i].xmin > x || vb_list[i].xmax < x) 
      continue; 
    else if (vb_list[i].ymin > y || vb_list[i].ymax < y) 
      continue; 
    else if (vb_list[i].zmin > z || vb_list[i].zmax < z) 
      continue; 

    // time is a bit of a hack
    // using original volume bounds without ghost (vo_list instead of vb_list)
    else if (vbr_list[i].tmin > t || vbr_list[i].tmax < t) 
      continue; 

    return(i); 
 
 }

  return(-1); 

}
//---------------------------------------------------------------------------
//
// Find the global rank number of the partition
// given my local block number
//
int Lattice4D::GetRank(int block) {

  return block_ranks[block];

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
// find the indices of the lattice element that has its subdomain number = rank 
//
// from old Lattice.C
//
int Lattice4D::GetIndices(int rank, int &iidx, int &jidx, int &kidx) {
  
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
int Lattice4D::GetBounds(int i, int j, int k, int t, volume_bounds_t& vb)  {

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
int Lattice4D::GetBounds(int rank, volume_bounds_t &vb) {

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
  volume_bounds_t vb = vb_list[idx]; 

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
// check if the point (x,y,z) is in the lattice element [i,j,k]
//
// from old Lattice.C
//
bool Lattice4D::isIn(float x, float y, float z, int i, int j, int k) {

  if (i < 0 || i > idim - 1 || j < 0 || j > jdim - 1 || k < 0 || k > kdim - 1) 
    return(false); 

  int idx = k * idim * jdim + j * idim + i; 
  volume_bounds_t vb = vb_list[idx]; 

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
// returns neighbor number (0-5) of x,y,z point
//
// from old Lattice.C
//
int Lattice4D::CheckNeighbor(int myrank, float x, float y, float z) {

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
int Lattice4D::GetNeighbor(int myrank, float x, float y, float z, float t, 
                           int &ei, int &ej, int &ek, int &et) {

  int neighbor_rank;
  neighbor_rank = CheckNeighbor(myrank, x, y, z, t);
  GetIndices(neighbor_rank, ei, ej, ek, et); 
  return (neighbor_rank); 

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

  // do not include ghost cells
  r = GetRealRank(x, y, z, t);

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
	    AddNeighbor(block, nr, part->parts[nr].Proc, block_ranks, 
			alloc_neighbors, &neighbor_ranks,
			&neighbor_procs, part);
	}
      }
    }
  }

}
//---------------------------------------------------------------------------
//
// gets local subvolume bounds (counting ghost cells)
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
// gets local subvolume bounds (not counting ghost cells)
//
// block: local block number (0-nblocks)
// min_s, max_s: (output) spatial min and max bounds
// min_t, max_t: (output) temporal min and max bounds
//
void Lattice4D::GetRealVB(int block, float *min_s, float *max_s, 
			  int *min_t, int *max_t) {

  min_s[0] = vbr_list[block_ranks[block]].xmin;
  min_s[1] = vbr_list[block_ranks[block]].ymin;
  min_s[2] = vbr_list[block_ranks[block]].zmin;
  max_s[0] = vbr_list[block_ranks[block]].xmax;
  max_s[1] = vbr_list[block_ranks[block]].ymax;
  max_s[2] = vbr_list[block_ranks[block]].zmax;
  *min_t = vbr_list[block_ranks[block]].tmin;
  *max_t = vbr_list[block_ranks[block]].tmax;
}
//---------------------------------------------------------------------------
//
// gets local time bounds w/o ghost
//
// block: local block number (0-nblocks)
// min_t, max_t: (output) temporal min and max bounds
//
void Lattice4D::GetTB(int block, int *min_t, int *max_t) {

  *min_t = vbr_list[block_ranks[block]].tmin;
  *max_t = vbr_list[block_ranks[block]].tmax;

}
//---------------------------------------------------------------------------
//
// gets the bounds of a time group
//
// group: the time group
// min_t, max_t: (output) temporal min and max bounds
//
void Lattice4D::GetTimeGroupBounds(int group, int* min_t, int* max_t) {

  *min_t = tb_list[group].tmin;
  *max_t = tb_list[group].tmax;
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
// calculates a default partition of the data
// by factoring the total number of blocks in alternating x, y, z directions
//
// data dim: x, y, z, t data size
// ghost: ghost layer per side
// nsp: number of spatial partitions
// ntp: number of temporal partitions
// lat_dim (output): x, y, z lattice dimensions
// returns: pointer to volume bounds list
//
void Lattice4D::ComputePartition(int *data_dim, int ghost, 
				 int nsp, int ntp, int *lat_dim) { 

  int rem = nsp; // unfactored remaining portion of nsp
  int block_dim[3]; // current block size
  int d[4]; // delta x, y, z, t (block size)
  int i, j, k, l, n;
  int max_dir; // longest remaining direction (0, 1, 2)

  // init
  vbr_list = new volume_bounds_t[nsp * ntp]; // volume bounds w/o ghost
  vb_list = new volume_bounds_t[nsp * ntp];  // volume bounds w/ ghost
  tb_list = new time_bounds_t[ntp];          // time partition bounds w/ ghost
  for (i = 0; i < 3; i++) {
    lat_dim[i] = 1;
    block_dim[i] = data_dim[i];
  }

  // compute factorization of data dimensions into lattice dimensions
  while (1) {

    // find longest division direction max_dir = 0, 1, or 2 (x, y, or z)
    max_dir = 0;
    for(i = 1; i < 3; i++) {
      if (block_dim[i] > block_dim[max_dir])
	max_dir = i;
    }

    // smallest factor remaining gets assigned to this direction
    for (j = 2; j <= rem; j++) {
      if (rem % j == 0) {
	lat_dim[max_dir] *= j;
	block_dim[max_dir] /= j;
	rem /= j;
	break;
      }
    }

    if (rem == 1)
      break;

    if (j > rem)
      fprintf(stderr,"Unable to partition the volume into %d spatial blocks. Please select a different number of spatial blocks and rerun.\n", nsp);
    assert(j <= rem);

  }

  // debug: print the lattice dims
//   fprintf(stderr,"lat_dim = [%d %d %d]\n", lat_dim[0], lat_dim[1], lat_dim[2]);

  // deltas
  for(i = 0; i < 3; i++) // x, y, z
	// MOD-BY-LEETEN 07/01/2011-FROM:
		// d[i] = roundf((float)data_dim[i] / (float)lat_dim[i]);
	// TO:
    d[i] = floorf(0.5f + (float)data_dim[i] / (float)lat_dim[i]);
	// MOD-BY-LEETEN 07/01/2011-END
  d[3] = data_dim[3] / ntp; // t

  // volume bounds in row-major index order (x, y, z, t)
  // x changes fastest, t slowest
  n = 0;
  for (l = 0; l < ntp; l++) {
    for (k = 0; k < lat_dim[2]; k++) {
      for (j = 0; j < lat_dim[1]; j++) {
	for (i = 0; i < lat_dim[0]; i++) {

	  vbr_list[n].xmin = i * d[0];
	  vbr_list[n].xmax = (i == lat_dim[0] - 1 ? data_dim[0] - 1 : 
			     (i + 1) * d[0]);
	  vbr_list[n].ymin = j * d[1];
	  vbr_list[n].ymax = (j == lat_dim[1] - 1 ? data_dim[1] - 1 : 
			     (j + 1) * d[1]);
	  vbr_list[n].zmin = k * d[2];
	  vbr_list[n].zmax = (k == lat_dim[2] - 1 ? data_dim[2] - 1 : 
			     (k + 1) * d[2]);
	  vbr_list[n].tmin = l * d[3];
	  vbr_list[n].tmax = (l == ntp - 1 ? data_dim[3] - 1 : 
			     (l + 1) * d[3]);
	  if (data_dim[3] == 1) // static case
	    vbr_list[n].tmax = 0;

	  vb_list[n].xmin = vbr_list[n].xmin;
	  vb_list[n].xmax = vbr_list[n].xmax;
	  vb_list[n].ymin = vbr_list[n].ymin;
	  vb_list[n].ymax = vbr_list[n].ymax;
	  vb_list[n].zmin = vbr_list[n].zmin;
	  vb_list[n].zmax = vbr_list[n].zmax;
	  vb_list[n].tmin = vbr_list[n].tmin;
	  vb_list[n].tmax = vbr_list[n].tmax;


	  n++;

	}
      }
    }
    tb_list[l].tmin = vb_list[n-1].tmin;
    tb_list[l].tmax = vb_list[n-1].tmax;
  }

  ApplyGhost(ghost);

  // debug: print the vb list
#if 0
  if (myproc == 0) {
    for (i = 0; i < nsp * ntp; i++)
      fprintf(stderr, "vbr_list[%d] min = [%d %d %d %d] max = [%d %d %d %d]\n", i, vbr_list[i].xmin, vbr_list[i].ymin, vbr_list[i].zmin, vbr_list[i].tmin, vbr_list[i].xmax, vbr_list[i].ymax, vbr_list[i].zmax, vbr_list[i].tmax);
    for (i = 0; i < nsp * ntp; i++)
      fprintf(stderr, "vb_list[%d] min = [%d %d %d %d] max = [%d %d %d %d]\n", i, vb_list[i].xmin, vb_list[i].ymin, vb_list[i].zmin, vb_list[i].tmin, vb_list[i].xmax, vb_list[i].ymax, vb_list[i].zmax, vb_list[i].tmax);
 }
#endif

}
//----------------------------------------------------------------------------
//
// adds ghost cell layers to block boundaries
//
// ghost: ghost layer per side
void Lattice4D::ApplyGhost(int ghost)
{
  for (int i = 0; i < npart; i++) {
    // Note: adding ghost cells to tmin is never necessary
    vb_list[i].xmin = max(vbr_list[i].xmin - ghost, 0);
    vb_list[i].ymin = max(vbr_list[i].ymin - ghost, 0);
    vb_list[i].zmin = max(vbr_list[i].zmin - ghost, 0);

    vb_list[i].xmax = min(vbr_list[i].xmax + ghost, xdim - 1);
    vb_list[i].ymax = min(vbr_list[i].ymax + ghost, ydim - 1);
    vb_list[i].zmax = min(vbr_list[i].zmax + ghost, zdim - 1);

    if(ldim != 1)  // don't add ghost cells when doing streamlines
      vb_list[i].tmax = min(vbr_list[i].tmax + ghost, ldim - 1);
  }
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
// void Lattice4D::SetReq(int block) {
//   part->SetReq(block_ranks[block]);
// }
// //
// // clears the request status
// //
// void Lattice4D::ClearReq(int block) {
//   part->ClearReq(block_ranks[block]);
// }
// //
// // gets the request status
// //
// int Lattice4D::GetReq(int block) {
//   return part->GetReq(block_ranks[block]); 
// }
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
// //
// // sets the computed status
// //
// void Lattice4D::SetComp(int block, int iter) {
//   part->SetComp(block_ranks[block], iter); 
// }
// //
// // clears the computed status
// //
// void Lattice4D::ClearComp(int block) {
//   part->ClearComp(block_ranks[block]); 
// }
// //
// // gets the computed status
// //
// int Lattice4D::GetComp(int block, int iter) {
//   return part->GetComp(block_ranks[block], iter); 
// }
//
// posts a point for sending
//
void Lattice4D::PostPoint(int block, float *p, int recirc, int64_t seed_id) {

  // neighbor number of neighbor block
  int neighbor = GetNeighbor(block, p[0], p[1], p[2], p[3]);

  if (neighbor >= 0) {
    int nr = neighbor_ranks[block][neighbor]; // rank of neighbor block
    int i;

    // check if the block is one of mine, and which block number it is
    for (i = 0; i < nb; i++) {
      if (block_ranks[i] == nr)
        break;
    }

    // only post points that remain inside the overall domain boundary
    // and that leave my own block
    if (neighbor >= 0 && (recirc || i != block)) {
      part->PostPoint(block_ranks[block], p, neighbor, seed_id); 
      if (neighbor_procs[block][neighbor] != myproc)
        tot_pts_send++;
    }
  }
}
// //
// // prints the posted points
// //
// void Lattice4D::PrintPost(int block) {
//   part->PrintPost(block_ranks[block]); 
// }

// #ifdef _MPI
// //
// // DEPRECATED
// //
// // exchanges points with all neighbors (old MPI version)
// // returns total number of points received by this process
// //
// int Lattice4D::OldExchangeNeighbors(float ***seeds, int *alloc_seeds,
// 				    int *num_seeds, int64_t **seed_ids) { 
//   int n;
//   comm_time = MPI_Wtime();
//   n = part->OldExchangeNeighbors(block_ranks, neighbor_ranks, neighbor_procs, 
// 			      seeds, alloc_seeds, num_seeds, MPI_COMM_WORLD, seed_ids);
//   comm_time = MPI_Wtime() - comm_time;
//   return n;
// }
// //
// // exchanges points with all neighbors (stable, synchronous MPI version)
// // returns total number of points received by this process
// //
// int Lattice4D::SyncExchangeNeighbors(int *nproc, float ***pts, int ***cts,
// 				    int64_t ***pids) {

//   int n;
//   comm_time = MPI_Wtime();
//   n = part->SyncExchangeNeighbors(block_ranks, neighbor_ranks, 
// 				  neighbor_procs, nproc, pts, cts, pids, 
// 				  MPI_COMM_WORLD);
//   comm_time = MPI_Wtime() - comm_time;
//   return n;
// }
// #endif

//---------------------------------------------------------------------------
//
// AddNeighbor is not a member of Lattice4D so that it can be used
// as a callback function
//
// adds a neighbor to the end of the neighbor_ranks table and sets its value
// updates the local block neighbors as well as the global partition
// grows local block neighbors if necessary
//
// myblock: my local block number
// neighrank: global partition number of the neighbor
// neighproc: process id of neighbor
// block_ranks: global ranks of local blocks
// alloc_neighbors: number of neighbors allocated for each local block
// neighbor_ranks: global ranks of neighbors of local blocks
// neighbor_procs: process ids of neighbors of local blocks
// part: partition data structure
//
void AddNeighbor(int myblock, int neighrank, int neighproc, int *block_ranks,
		 int *alloc_neighbors, int ***neighbor_ranks,
		 int ***neighbor_procs, Partition *part) {

  int myrank = block_ranks[myblock];
  int num_neighbors = part->parts[block_ranks[myblock]].NumNeighbors;
  int old_n = alloc_neighbors[myblock];  // num neighbors currently alloc'd
  int n; // new number of neighbors allocated

  if (old_n < num_neighbors + 1) {

    if (old_n == 0) {
      n = 1;
      assert(((*neighbor_ranks)[myblock] =
	      (int *)malloc(n * sizeof(int))) != NULL);
      assert(((*neighbor_procs)[myblock] = 
	      (int *)malloc(n * sizeof(int))) != NULL);
    }
    else {
      n = old_n * 2;
      assert(((*neighbor_ranks)[myblock] = 
	      (int *)realloc((*neighbor_ranks)[myblock],
			     n * sizeof(int))) != NULL);
      assert(((*neighbor_procs)[myblock] = 
	      (int *)realloc((*neighbor_procs)[myblock],
			     n * sizeof(int))) != NULL);

    }

    alloc_neighbors[myblock] = n;

  }

  (*neighbor_ranks)[myblock][num_neighbors] = neighrank;
  (*neighbor_procs)[myblock][num_neighbors] = neighproc;

  part->AddNeighbor(myrank);

}
//---------------------------------------------------------------------------

// utility functions

//---------------------------------------------------------------------------
//
// swap4(n)
//
// Swaps 4 bytes from 1-2-3-4 to 4-3-2-1 order.
// cast the input as a char and use on any 4 byte variable
//
void Lattice4D::swap4(char *n) {

  char *n1;
  char c;

  n1 = n + 3;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//----------------------------------------------------------------------------
