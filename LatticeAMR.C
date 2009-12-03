//------------------------------------------------------------------------------
//
// AMR lattice
//
// Copyright (c) 2009 Han-Wei Shen and Tom Peterka
//
// Contact:
//
// Han-Wei Shen
// The Ohio State University
// Columbus, OH
//
// Tom Peterka
// MCS Radix Lab
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// All rights reserved. May not be used, modified, or copied
// without permission
//
//--------------------------------------------------------------------------

#include "LatticeAMR.h"
#include "FlashAMR.h"

//--------------------------------------------------------------------------
//
// Lattice AMR
//
// constructs and initializes a time-varying regular structured 
// AMR lattice 
//
// filename: file containing list of timestep files
// tlen: number of timesteps
// vx, vy, vz: names of velocity components in the dataset
// nid: number of processes, threads, owners
// default = 1 (can omit if single process sequential program)
// myid: rank, process number, thread number, identification of the owner
// default = 0 (can omit if single process sequential program)
//
LatticeAMR::LatticeAMR(char *filename, int tlen, char *vx, char *vy, char *vz,
		       int nid, int myid) {

  float min[3], max[3]; // spatial extents
  float blockSize[3]; // physical size of a block
  float levelMinB[3], levelMaxB[3]; // extents in a level
  float *center; // spatial center of block
  int start_block, end_block; // start and end block numbers of my blocks
                              // in a contiguous range distribution
  float *data; // pointer to data
  FlashAMR *amr; // AMR object for one time step
  int level; // current level
  int idx; // index in level
  int rank; // partition rank of one block
  int p; // process id for a block
  int i, j, n, t;
  double t0; // temporary timer for performance instrumentation

  // init AMR
  tamr = new TimeVaryingFlashAMR(myid);
  t0 = MPI_Wtime();
  tamr->LoadMetaData(filename, min, max); 
  io_time = MPI_Wtime() - t0;
  num_levels = tamr->GetNumLevels();
  tamr->GetDims(block_dims); 

  // lattice overall bounds
  xdim = max[0]-min[0]; 
  ydim = max[1]-min[1]; 
  zdim = max[2]-min[2]; 
  tdim = tlen;  // number of time steps

  // the bounds of each level 
  xmin = new float[num_levels];   xmax = new float[num_levels]; 
  ymin = new float[num_levels];   ymax = new float[num_levels]; 
  zmin = new float[num_levels];   zmax = new float[num_levels]; 
  tmin = new int [num_levels];    tmax = new int[num_levels]; 

  // the xyzt lengths of a block in each level, assume the same for all blocks
  xlength = new float[num_levels]; 
  ylength = new float[num_levels]; 
  zlength = new float[num_levels]; 
  tlength = new int[num_levels]; 

  // the data resolution of blocks in each level, same for all blocks
  xres = new int[num_levels]; 
  yres = new int[num_levels]; 
  zres = new int[num_levels]; 
  tres = new int[num_levels]; 

  // the dimensions of the lattice in each level for partitioning purpose 
  idim = new int[num_levels]; 
  jdim = new int[num_levels]; 
  kdim = new int[num_levels]; 
  ldim = new int[num_levels];  

  // data structures for book-keeping 
  has_data = new bool*[num_levels]; // whether the lattice element has data
  has_data_from_merger = new int*[num_levels]; // used for block merger 
  finest_level = new int*[num_levels];   // what is the finest level of data 
                                         // in this spatial region
  data_ptr = new float**[num_levels];    // one float ptr per time step 
  index_to_rank = new int*[num_levels];  // mapping from index to rank 
  index_to_seq = new int*[num_levels];  // mapping from index to block sequence
                                        // number
  nblocks = new int[num_levels];         // how many blocks in each region 
                                         // regardless of empty or not 
  // process info
  myproc = myid;
  nproc = nid;

  // create AMR levels
  for (i = 0; i < num_levels; i++)  {

    tamr->GetLevelBlockSize(i, blockSize); 
    tamr->GetLevelBounds(i, levelMinB, levelMaxB); 

#ifdef DEBUG
    if (myproc == 0)
      fprintf(stderr, "Level %d: physical size of one block in this level = [%.4e %.4e %.4e] min=[%.4e %.4e %.4e %d] max=[%.4e %.4e %.4e %d]\n", 
	      i, blockSize[0], blockSize[1], blockSize[2], 
	      levelMinB[0], levelMinB[1], levelMinB[2], 
	      levelMaxB[0], levelMaxB[1], levelMaxB[2],
	      0, tlen - 1); 
#endif

    CreateLevel(i, blockSize[0], blockSize[1], blockSize[2], 
		block_dims[0], block_dims[1], block_dims[2], 
		levelMinB[0], levelMaxB[0], levelMinB[1], levelMaxB[1], 
		levelMinB[2], levelMaxB[2], 0, tlen - 1); 

  }

  // check in each block, without the data for now
  for (t = 0; t < tlen; t++) {
    amr = tamr->GetTimeStep(t);
    for (i = 0; i < amr->GetNumBlocks(); i++) {
      center = amr->GetBlockCenter(i); 
      CheckIn(amr->GetLevel(i), center[0], center[1], center[2], t, 0);
    }
  }

  // assign a contiguous distribution of blocks to processes
  AssignContig(0, tlen - 1, &start_block, &end_block);

  // final bookkeeping and combining any blocks that are possible
  CompleteLevels(tlen);

  // create partition data structure
  part = new Partition(npart, nproc);

  // read my data blocks
  t0 = MPI_Wtime();
  tamr->LoadData(filename, start_block, end_block, vx, vy, vz);
  io_time += (MPI_Wtime() - t0);

  // check in actual data into my blocks
  // record partition data struture info for all blocks
  n = 0; // block number in the entire sequence so far
  for (t = 0; t < tlen; t++) {

    amr = tamr->GetTimeStep(t);
    for (i = 0; i < amr->GetNumBlocks(); i++) {

      // identify the block by its center -> level and index in level
      center = amr->GetBlockCenter(i); 
      level = amr->GetLevel(i);
      assert((idx = GetIndexinLevel(level, center[0], center[1], center[2],
				    (float)t)) != -1);

      // get and save process id
      rank = index_to_rank[level][idx];
      p = GetProc(level, idx);
      part->parts[rank].Proc = p;

      // update proc_parts and check in data for my own blocks only
      if (n >= start_block && n <= end_block) {

	// check if we have this rank already
	for (j = 0; j < part->proc_nparts[p]; j++) {
	  if (part->proc_parts[p][j] == rank)
	    break;
	}
	// if the block is in a new partition rank, 
	// add it to the proc_parts structure
	if (j == part->proc_nparts[p]) {
	  part->proc_parts[p][part->proc_nparts[p]] = rank;
	  part->proc_nparts[p]++;
	}
	// check in the data of the block
	data_ptr[level][idx] = amr->GetDataPtr(i); 

      }

      n++;

    }

  }

  // proc_parts needs to be sorted in same order as blocks (ascending order)
  for (i = 0; i < nid; i++) {
    if (part->proc_nparts[i] > 0)
      int_sort_list(part->proc_parts[i],part->proc_nparts[i]);
  }

  // ranks of my blocks and my neighbors' blocks
  nb = GetMyNumPartitions(); // my number of blocks

#ifdef DEBUG
  fprintf(stderr,"Process %d has %d partitions\n", myproc, nb);
#endif

  if (myproc >= 0) {

    // allocate table of neighbor ranks for one neighbor initially
    assert((neighbor_ranks = (int **)malloc(nb * sizeof(int *))) != NULL);
    for (i = 0; i < nb; i++)
      assert((neighbor_ranks[i] = (int *)malloc(sizeof(int))) != NULL);

    // ranks of my blocks
    assert((block_ranks = (int *)malloc(nb * sizeof(int))) != NULL);
    GetMyPartitions(myproc, block_ranks);

    // learn who my neighbors are
    for (i = 0; i < nb; i++)
      GetNeighborRanks(i);

  }

  // average number of neighbors for each block
  int tot_neighbors = 0;
  for (i = 0; i < nb; i++)
    tot_neighbors += part->parts[i].NumNeighbors;
  avg_neigh = tot_neighbors / nb;

}
//------------------------------------------------------------------------------
//
//  Destructor
// 
LatticeAMR::~LatticeAMR()
{
  delete []xmin;   delete []ymin;   delete []zmin; delete []tmin;  
  delete []xmax;   delete []ymax;   delete []zmax; delete []tmax; 
  delete []idim;   delete []jdim;   delete []kdim; delete []ldim; 
  delete []xlength;  delete []ylength;  delete []zlength; delete []tlength; 
  delete [] rank_to_index; 
  delete [] nblocks; 

  for (int i=0; i<num_levels; i++) {
    delete [] has_data[i]; 
    delete [] finest_level[i]; 
    delete [] index_to_rank[i]; 
    delete [] index_to_seq[i]; 
  }
  delete [] has_data; 
  delete [] finest_level; 
  delete [] index_to_rank; 
  delete [] index_to_seq; 

  delete part;

  if (vb_list!=NULL) delete [] vb_list; 
  if (seedlists!=NULL) delete [] seedlists; 

  delete tamr;

}
//------------------------------------------------------------------------------
//
// Create a new level of lattice. The level has to be smaller than the 
// largest level set in the constructor 
//
// level: the level id [0..max_level-1]
// {x,y,z}_size: the lengths of the block in this level 
// {x,y,z,t}_res: the data resolution in each dimension 
// {x,y,z,t}_{min,max}: the physical bounds of the level 
// 
bool LatticeAMR::CreateLevel(int level, float x_size, float y_size, 
			     float z_size, 
			     int x_res, int y_res, int z_res, 
			     float x_min, float x_max, 
			     float y_min, float y_max, float z_min, 
			     float z_max, int t_min, int t_max) 
{

  if (level <0 || level >=num_levels) return false; 

  xlength[level]=x_size; ylength[level]=y_size; zlength[level]=z_size;
  tlength[level]= 1;  // one time step in each AMR block 
  
  // the data resolution 
  xres[level] = x_res; yres[level] = y_res; zres[level] = z_res; 
  tres[level] = 1;      // data always come one time step at a time 
 
  // the physical bounds of this level 
  xmin[level] = x_min; xmax[level] = x_max; 
  ymin[level] = y_min; ymax[level] = y_max; 
  zmin[level] = z_min; zmax[level] = z_max;
  tmin[level] = t_min; tmax[level] = t_max;  // the total time range 
  
  //i,j,k,l dim are the number of blocks in each dimension if the whole 
  //space-time domain is filled with blocks from this level 
  // this is the lattice dimensions for this level 
  // 
  idim[level] = (int)((x_max-x_min)/x_size); 
  jdim[level] = (int)((y_max-y_min)/y_size); 
  kdim[level] = (int)((z_max-z_min)/z_size); 
  ldim[level] = t_max-t_min+1; 

  // total number of theoretical blocks at this level
  // (not all have data)
  int size = idim[level]*jdim[level]*kdim[level]*ldim[level];
  nblocks[level] = size;

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"Level %d: theoretical lattice dims %dx%dx%dx%d=%d \n", 
	    level, idim[level], jdim[level], kdim[level], ldim[level], size); 
#endif

  has_data[level] = new bool[size]; 
  has_data_from_merger[level] = new int[size]; 
  finest_level[level] = new int[size]; 
  data_ptr[level] = new float*[size]; 
  index_to_rank[level] = new int[size];
  index_to_seq[level] = new int[size];

  int idx = 0; 
  for (int t=0; t<ldim[level]; t++) 
    for (int k=0; k<kdim[level]; k++) 
      for (int j=0; j<jdim[level]; j++) 
	for (int i=0; i<idim[level]; i++) {
	  has_data[level][idx] = false;    // initial value: no data
	  has_data_from_merger[level][idx] = -1; 
	  finest_level[level][idx] = -1;   // initial value: unknown
	  data_ptr[level][idx] = NULL; 
	  index_to_rank[level][idx] = -1;  // no mapping available 
	  idx++; 
	}
  return true; 
}
//------------------------------------------------------------------------------
//
//  Get the index of the block at the given level that contains (x,y,z,t)
//
int LatticeAMR::GetIndexinLevel(int level, float x, float y, float z, float t) {

  if (level < 0 || level >= num_levels || 
      x < xmin[level] || x > xmax[level] || 
      y < ymin[level] || y > ymax[level] || 
      z < zmin[level] || z > zmax[level] || 
      t < tmin[level] || t > tmax[level])
    return -1; 
  
  int i = (int)((x - xmin[level]) / xlength[level]); 
  int j = (int)((y - ymin[level]) / ylength[level]); 
  int k = (int)((z - zmin[level]) / zlength[level]); 
  int l = (int)((t - tmin[level]) / tlength[level]); 

  // clamp max edge to be inside last block
  i = (i == idim[level] ? i - 1 : i);
  j = (j == jdim[level] ? j - 1 : j);
  k = (k == kdim[level] ? k - 1 : k);
  l = (l == ldim[level] ? l - 1 : l);

  int idx = l * idim[level] * jdim[level] * kdim[level] + 
    k * idim[level] * jdim[level] + j * idim[level] + i; 

  return idx; 

}
//------------------------------------------------------------------------------
//
// AssignContig
//
// assigns a contiguous partitioning of blocks over a time interval
// sets start and end block numbers for a process
//
// start_time, end_time: time interval in question
// start_block, end_block: (output) starting and ending block range of this proc
// returns: total number of blocks
//
int LatticeAMR::AssignContig(int start_time, int end_time,
			     int *start_block, int *end_block) {

  int blocks_per_proc; // number of blocks per process
  FlashAMR *amr; // AMR object for one time step
  int i;

  // total number of blocks
  tot_nblocks = 0;
  for (i = 0; i < end_time - start_time + 1; i++) {
    amr = tamr->GetTimeStep(i);
    tot_nblocks += amr->GetNumBlocks();
  }

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"\nTotal number of leaf blocks = %d\n\n", tot_nblocks);
#endif

  // start and end blocks
  blocks_per_proc = tot_nblocks / nproc;
  *start_block = myproc * blocks_per_proc;
  *end_block = *start_block + blocks_per_proc - 1;
  if (myproc == nproc - 1)
    *end_block = tot_nblocks - 1;

  return tot_nblocks;

}
//------------------------------------------------------------------------------
//
// GetProc
//
// gets the process id. corrsponding to a block level and index
// assumes block has already been checked in
// assumes contiguous block division to map block sequence number to process
//
int LatticeAMR::GetProc(int level, int idx) {

  int n; // sequence number of a block
  int blocks_per_proc; // number of blocks per process
  int p; // process id

  // convert level and idx to the block sequence number
  n = index_to_seq[level][idx];

  // divide the blocks contiguously among processers
  blocks_per_proc = tot_nblocks / nproc;
  p = n / blocks_per_proc;
  if (p == nproc) // edge case, last process has remainder blocks
    p--;

  return p;

}
//------------------------------------------------------------------------------
//
// CheckIn
//
// Indicates that the lattice element containing x, y, z, t
// has data at this level, or will have it before completing the lattice setup
// 
// It also updates the corresponding lattice element in other level 
// regarding which level has the finest level of data 
//
// Each process must check in all blocks in the entire lattice, 
// in order to have a consistent view of the entire structure
// even if the data will actually come from a different process
//
// If data will be provided later or by another process, pass data = 0
//
bool LatticeAMR::CheckIn(int level, float x, float y, float z, int t, 
			 float* data) {

  int idx; // index in level
  static int n; // sequence number of block (in the order CheckIn was called)

  if ((idx = GetIndexinLevel(level, x, y, z, (float)t)) == -1) {
    fprintf(stderr, "panic: idx = -1\n");
    return false;
  }

  has_data[level][idx] = true; 
  data_ptr[level][idx] = data; 
  index_to_seq[level][idx] = n;

  // update the blocks at other levels as well
  for (int i = 0; i < num_levels; i++){
    idx = GetIndexinLevel(i,x,y,z, (float) t); 
    if (idx != -1) 
      if (finest_level[i][idx] < level)
	finest_level[i][idx] = level; 
  }

  n++;
  return true; 

}
//----------------------------------------------------------------------------
//
// CompleteLevels
//
// Final bookkeeping
//
// Breaks the time range into segments of max length t_interval
// A new segment results if the length > than t_interval, or there is no data
// in the corresponding spatial region 
// 
// stores the final number of partitions after simplifying into time segments
// in the class variable nparts
//
void LatticeAMR::CompleteLevels(int t_interval) {

  int idx; // index in level
  int old_proc; // process id for previous block
  int counter;
  int offset;
  int i, j, k, l, t;

  // count the number of blocks before any simplification
  npart = 0; 
  for (i=0; i<num_levels; i++)
    for (j=0; j<nblocks[i]; j++)
      if (has_data[i][j]) npart++; 

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"Before combining any time intervals, npart = %d\n", npart); 
#endif

  // allocate a few data structures for the max number of partitions
  vb_list = new volume_bounds_type_f[npart]; 
  rank_to_index = new int[npart]; 

  // init
  npart = -1; 
  offset = 0; 
  old_proc = -1;

  // loop through all spatial regions in all levels
  for (l = 0; l < num_levels; l++) {
    for (k = 0; k < kdim[l]; k++) 
      for (j = 0; j < jdim[l]; j++) 
	for (i = 0; i < idim[l]; i++) {

	  counter = 0; 

	  for (t = 0; t < ldim[l]; t++) {

	    idx = t * idim[l] * jdim[l] * kdim[l] + k * idim[l] * jdim[l] +
	      j * idim[l] + i; 

	    if (has_data[l][idx]) {

	      // reset the counter. ready to start another time segment 
	      if (counter >= t_interval || GetProc(l, idx) != old_proc)
		counter = 0;

	      if (counter == 0) {    // a new time segment starts

		npart++; 

		vb_list[npart].xmin= xmin[l] +i * xlength[l]; 
		vb_list[npart].xmax= xmin[l] + (i + 1) * xlength[l]; 
		vb_list[npart].ymin= ymin[l] + j * ylength[l]; 
		vb_list[npart].ymax= ymin[l] + (j + 1) * ylength[l]; 
		vb_list[npart].zmin= zmin[l] + k * zlength[l]; 
		vb_list[npart].zmax= zmin[l] + (k + 1) * zlength[l]; 
		vb_list[npart].xdim = xres[l]; 
		vb_list[npart].ydim = yres[l]; 
		vb_list[npart].zdim = zres[l]; 
		vb_list[npart].tmin= tmin[l] + t * tlength[l]; 
		vb_list[npart].tmax= tmin[l] + t * tlength[l]; 
		vb_list[npart].tdim = 1; 

		index_to_rank[l][idx] = npart; 
		rank_to_index[npart] = offset + idx; 
		old_proc = GetProc(l, idx);

		// update the counter for next block
		if (t + 1 == ldim[l])
		  counter = 0; 
		else 
		  counter++; 

	      }

	      else  { // add to existing time segment
		vb_list[npart].tmax = tmin[l] + t * tlength[l]; 
		vb_list[npart].tdim++; 
		index_to_rank[l][idx] = npart; 
		counter++; 
	      }

	    } // has data

	    else  // do not have data
	      counter = 0;

	  }

	}

    offset += nblocks[l]; 

  }

  npart++;

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"After combining any time intervals, final npart = %d\n", npart); 
#endif

}
//----------------------------------------------------------------------------

///////////////////////////////////////////////////////////////////////
//
//  Call this function after all blocks with data have checked in
//  Go through all levels and collect blocks that have data 
//
void LatticeAMR::MergeAndCompleteLevels()
{
  npart = 0; 
  // first check how many non-empty blocks
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data_from_merger[i][j]!=-1) npart++; 
  // next allocate the volume bounds type etc. 
  vb_list = new volume_bounds_type_f[npart]; 
  rank_to_index = new int[npart]; 

  int isize, jsize, ksize; 
  float x_min, x_max, y_min, y_max, z_min, z_max; 
  int t_min, t_max; 
  int offset = 0; 
  int rank = 0; 
  for (int level=0; level<num_levels; level++) {
    isize = idim[level]; 
    jsize = jdim[level]; 
    ksize = kdim[level]; 

    x_min = xmin[level]; x_max = xmax[level]; 
    y_min = ymin[level]; y_max = ymax[level]; 
    z_min = zmin[level]; z_max = zmax[level]; 
    t_min = tmin[level]; t_max = tmax[level]; 
    // j is to index all blocks in this level
    // rank is to index non-empty blocks in all levels 
    for (int j=0; j<nblocks[level]; j++) {
      int iidx, jidx, kidx, tidx; 
      if (has_data_from_merger[level][j]!=-1) {
	tidx = j / (isize*jsize*ksize); 
	int r = j % (isize*jsize*ksize); 
	kidx = r / (isize * jsize) ; 
	jidx = (r % (isize * jsize)) / isize; 
	iidx = r % isize;
	
	vb_list[rank].xmin= x_min+iidx*xlength[level]; 
	vb_list[rank].xmax= x_min+(iidx+1)*xlength[level]; 
	vb_list[rank].ymin= y_min+jidx*ylength[level]; 
	vb_list[rank].ymax= y_min+(jidx+1)*ylength[level]; 
	vb_list[rank].zmin= z_min+kidx*zlength[level]; 
	vb_list[rank].zmax= z_min+(kidx+1)*zlength[level]; 
	vb_list[rank].tmin= t_min+tidx*tlength[level]; 
	vb_list[rank].tmax= t_min+(tidx+1)*tlength[level]; 

	// need to calculate the resolutions and also
	// update the data pointer  

	int DLevel = has_data_from_merger[level][j]; 
	if (DLevel==level){
	  vb_list[rank].xdim = xres[level]; 
	  vb_list[rank].ydim = yres[level]; 
	  vb_list[rank].zdim = zres[level]; 
	  vb_list[rank].tdim = tres[level]; 
	}
	else {

	  int Imin, Imax, Jmin, Jmax, Kmin, Kmax, Tmin, Tmax; 
	  MapCells(iidx, jidx, kidx, tidx, level, DLevel, Imin, Imax, 
		   Jmin, Jmax, Kmin, Kmax, Tmin, Tmax); 
	  int x_res = xres[DLevel]*(Imax-Imin); 
	  int y_res = yres[DLevel]*(Jmax-Jmin); 
	  int z_res = zres[DLevel]*(Kmax-Kmin); 
	  int t_res = tres[DLevel]*(Tmax-Tmin); 
	  
	  float *data = new float[x_res*y_res*z_res*t_res*3]; 

	  if (data == NULL) {
	    printf(" Panic. cannot allocate memory. \n"); 
	    exit(1); 
	  }
	  printf(" *** %d %d %d %d %d %d %d %d\n", Imin, Imax, Jmin, Jmax, 
		 Kmin, Kmax, Tmin, Tmax); 
	  for (int tFor = Tmin; tFor<Tmax; tFor++) 
	    for (int iFor=Imin; iFor<Imax; iFor++)
	      for (int jFor=Jmin; jFor<Jmax; jFor++)
		for (int kFor=Kmin; kFor<Kmax; kFor++) {
		  int idx = tFor*idim[DLevel]*jdim[DLevel]*kdim[DLevel]+kFor*idim[DLevel]*jdim[DLevel]+jFor*idim[DLevel]+iFor; 
		  // now copy data over 
		  if (data_ptr[DLevel][idx] == NULL) {
		    printf(" Panic. level %d [%d %d %d %d] null data ptr\n", DLevel, 
			   iFor, jFor, kFor, tFor); 
		    exit(1); 
		  }
		  if (has_data[DLevel][idx] == false) {
		    printf(" Panic. level %d has No data!\n", DLevel); 
		    exit(1); 
		  }
		  float* from_data = data_ptr[DLevel][idx]; 
		  for(int dtFor=0; dtFor<tres[DLevel]; dtFor++)
		    for(int dkFor=0; dkFor<zres[DLevel]; dkFor++)
		      for(int djFor=0; djFor<yres[DLevel]; djFor++)
			for(int diFor=0; diFor<xres[DLevel];diFor++) {
			  int source = dtFor*zres[DLevel]*yres[DLevel]*xres[DLevel]+dkFor*yres[DLevel]*xres[DLevel]+djFor*xres[DLevel]+diFor; 
			  int target_i = (iFor-Imin)*xres[DLevel]+diFor; 
			  int target_j = (jFor-Jmin)*yres[DLevel]+djFor; 
			  int target_k = (kFor-Kmin)*zres[DLevel]+dkFor; 
			  int target_t = (tFor-Tmin)*tres[DLevel]+dtFor; 
			  int target = target_t*x_res*y_res*z_res+target_k*x_res*y_res+target_j*x_res+target_i; 
			  data[target*3] = from_data[source*3]; 
			  data[target*3+1] = from_data[source*3+1]; 
			  data[target*3+2] = from_data[source*3+2]; 
			}
		  if (from_data!=NULL) 
		    delete [] from_data; // why cann't you free the data? 
		}

	  data_ptr[level][j] = data; 
	  vb_list[rank].xdim = x_res; 
	  vb_list[rank].ydim = y_res; 
	  vb_list[rank].zdim = z_res; 
	  vb_list[rank].tdim = t_res; 
	  }
	index_to_rank[level][j] = rank; 
	rank_to_index[rank] = offset + j; 
	rank++; 
	}
    }
    offset+= nblocks[level]; 
  }
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++) {
      if (has_data_from_merger[i][j]!=-1) has_data[i][j] = true; 
      else has_data[i][j]=false; 
    }

}

/////////////////////////////////////////////////////////////////
// 
//   Get the data pointer for the block of the given rank 
//

float** LatticeAMR::GetDataPtr(int rank) 
{
  if (rank < 0 || rank >=npart) return (NULL); 

  int n_steps = vb_list[rank].tdim; 
  float ** ppvector = new float*[n_steps]; 

  int remainder = rank_to_index[rank]; 
  int l; 
  for (l=0; l<num_levels; l++) {
    if (remainder-nblocks[l] < 0) break; 
    else 
      remainder -=nblocks[l]; 
  }
  if (l == num_levels) return(NULL); //rank is too big 

  int t_jump = idim[l]*jdim[l]*kdim[l]; 
  int index_offset = 0; 
  for (int i=0; i<n_steps; i++) {
    if (has_data[l][remainder+i*t_jump] == false) {
      printf("panic: has_data[%d][%d] is false, rank = %d remainder = %d \n",
	     l, remainder+i*t_jump, rank, remainder); 
      exit(1); 
    }
    ppvector[i] = data_ptr[l][remainder+i*t_jump]; 
    if (ppvector[i]==NULL) 
      printf("panic: data_ptr[%d][%d] is NULL\n",l,remainder+i*t_jump); 
  }
  return (ppvector); 
}


/////////////////////////////////////////////////////////////////
//
// find the block coords (i,j,k,l) at the given level that contains (x,y,z,t)
// This function does not cheeck whether the element has data or not  
//
int LatticeAMR::GetCoordsinLevel(int level, float x, float y, float z, float t, 
				  int &i, int&j, int&k, int&l)
{
  if (level <0 || level >=num_levels) return -1; 
  if (x<xmin[level] || x>xmax[level] || y<ymin[level] || y>ymax[level] || 
      z<zmin[level] || z>zmax[level] || t<tmin[level] || t>ymax[level])
    return -1; 
  
  i = (int)((x-xmin[level])/xlength[level]); 
  j = (int)((y-ymin[level])/ylength[level]); 
  k = (int)((z-zmin[level])/zlength[level]); 
  l = (int)((t-tmin[level])/tlength[level]); 
  int idx = l*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 
  return idx; 
}

//////////////////////////////////////////////////////////////////
//
//    Return the finest level with data that contains (x,y,z,t) 
//
int LatticeAMR::GetFinestLevel(float x, float y, float z, float t) 
{
  int idx; 
  for (int i=num_levels-1; i>=0; i--) {
    idx= GetIndexinLevel(i, x,y,z,t); 
    if (idx != -1) 
      if (has_data[i][idx] == true) return(i); 
  }
  return(-1); // no data at that point 
}

////////////////////////////////////////////////////////////////////
//
// Find the subdomain rank that contains the physical location (x,y,z,t) 
// with the fineset level of data 
//
int LatticeAMR::GetRank(float x, float y, float z, float t) {

  int idx; 
  for (int i=num_levels-1; i>=0; i--) {
    idx= GetIndexinLevel(i, x,y,z,t); 
    if (idx!=-1) {
      if (has_data[i][idx] == true) {
	int rank = index_to_rank[i][idx]; 
	return rank; 
      }
    }
  }
  return(-1); 
}

////////////////////////////////////////////////////////////////////
//
// Find the subdomain rank that contains the physical location (x,y,z,t) 
// at the given level
//
// returns -1 if the level does not contain data
//
int LatticeAMR::GetRank(float x, float y, float z, float t, int level) {

  int idx = GetIndexinLevel(level, x, y, z, t); 
  if (idx != -1 && has_data[level][idx] == true)
    return(index_to_rank[level][idx]);

  return(-1); 

}
//////////////////////////////////////////////////////////////////////
//
//  Get the subdomain rank for the lattice[i,j,k] element at 
//  the given level.
//  note: the lattice element may not contain data. In that case, 
//  the return value will be -1. 
//
int LatticeAMR::GetRank(int i, int j, int k, int t, int level) {

  if (level <0 || level >=num_levels) return -1; 
  if (i<0 || j<0 || k<0 || t<0 ||
      i>=idim[level] || j>=jdim[level] || k>=kdim[level] || t>=ldim[level]) 
    return -1; 
  
  int idx = t*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 

  return index_to_rank[level][idx]; 
}
////////////////////////////////////////////////////////////////////////
//
// find the indices and level of the lattice element that has its 
// subdomain number = rank 
//
int LatticeAMR::GetCoords(int rank, int &iidx, int &jidx, int &kidx, int& tidx, int& level) 
{
  if (rank < 0 || rank >=npart) return (-1); 

  int remainder = rank_to_index[rank]; 
  int l; 
  for (l=0; l<num_levels; l++) {
    if (remainder-nblocks[l] < 0) break; 
    else 
      remainder -=nblocks[l]; 
  }
  if (l == num_levels) return(-1); //rank is too big 

  int isize = idim[l]; 
  int jsize = jdim[l]; 
  int ksize = kdim[l]; 
  int tsize = ldim[l]; 

  tidx = remainder / (isize*jsize*ksize); 
  int r = remainder % (isize*jsize*ksize); 
  kidx = r / (isize * jsize) ; 
  jidx = (r % (isize * jsize)) / isize; 
  iidx = r % isize; 
  level = l;

  return(1); 
}

///////////////////////////////////////////////////////////////
//
// find the indices and level of the lattice element that has 
// the finest resolution of data at (x,y,z,t)
// it also returns the rank of the lattice element 
//
int LatticeAMR::GetCoords(float x, float y, float z, float t, 
			  int &iidx, int &jidx, int &kidx, int&tidx, int&level) 
{
  int idx, rank; 
  for (level=num_levels-1; level>=0; level--) {
    idx= GetCoordsinLevel(level, x,y,z,t, iidx, jidx, kidx,tidx); 
    if (idx !=-1) 
      if (has_data[level][idx] == true)  {
	rank = index_to_rank[level][idx]; 
	return rank; 
      }
  }
  return(-1); // no data at that point 
}

//////////////////////////////////////////////////////////////
//
//check if the point (x,y,z,t) is in the lattice element [i,j,k,l] 
// at the given level 
//
bool LatticeAMR::isIn(float x, float y, float z, float t, int i, int j, int k, 
		      int l, int level) {

  if (level <0 || level >= num_levels) return (-1); 
  if (i < 0 || i > idim[level] - 1 || j < 0 || j > jdim[level] - 1
      || k < 0 || k > kdim[level] - 1 || l <0 || l> ldim[level]-1) 
    return(false); 

  float min_x, max_x, min_y, max_y, min_z, max_z, min_t, max_t; 

  min_x = xmin[level] + i*xlength[level]; 
  max_x = xmin[level] + (i+1)*xlength[level]; 
  min_y = ymin[level] + i*ylength[level]; 
  max_y = ymin[level] + (i+1)*ylength[level]; 
  min_z = zmin[level] + i*zlength[level]; 
  max_z = zmin[level] + (i+1)*zlength[level]; 
  min_t = tmin[level] + i*tlength[level]; 
  max_t = tmin[level] + (i+1)*tlength[level]; 

  if (min_x > x || max_x < x) 
    return (false); 
  else if (min_y > y || max_y < y) 
    return (false); 
  else if (min_z > z || max_z < z) 
    return (false); 
  else if (min_t > t || max_t < t) 
    return(false); 

  return(true); 

}

///////////////////////////////////////////////////////////
//
// look up the volume bounds of lattice[i,j,k,t] at the given level 
//
int LatticeAMR::GetBounds(int i, int j, int k, int t, int level, volume_bounds_type_f& vb)  {

  if (level <0 || level >=num_levels) return (-1); 

  if (i < 0 || i >= idim[level]) 
    return(-1); 
  else if (j < 0 || j >= jdim[level]) 
    return(-1); 
  else if (k < 0 || k >= kdim[level]) 
    return(-1); 
  else if (t < 0 || t >=ldim[level]) 
    return(-1); 

  vb.xmin = xmin[level] + i*xlength[level]; 
  vb.xmax = xmin[level] + (i+1)*xlength[level]; 
  vb.ymin = ymin[level] + j*ylength[level]; 
  vb.ymax = ymin[level] + (j+1)*ylength[level]; 
  vb.zmin = zmin[level] + k*zlength[level]; 
  vb.zmax = zmin[level] + (k+1)*zlength[level]; 
  vb.tmin = tmin[level] + t*tlength[level]; 
  vb.tmax = tmin[level] + (t+1)*tlength[level]; 

  return(1); 
}

//////////////////////////////////////////////////////////////
//
// look up the volume bounds of the subdomain 'rank' 
//
int LatticeAMR::GetBounds(int rank, volume_bounds_type_f &vb)
{
  if (rank < 0 || rank >=npart) return(-1); 

  vb = vb_list[rank];
  return(1); 

}

///////////////////////////////////////////////////////////
//
// returns neighbor rank where x,y,z,t is in 
//
int LatticeAMR::CheckNeighbor(int myrank, float x, float y, float z, float t) {

  return GetRank(x,y,z,t); 

}
////////////////////////////////////////////////////////
//
// GetNeighbor
//
// Returns the indices and level of the neighbor that has data and contains the point (x,y,z,t)
// Also returns the rank of the element (-1 means no data found at the point in all levels)
// myrank: global subdomain id (not used any more) 
//
//
int LatticeAMR::GetNeighbor(int myrank, float x, float y, float z, float t, 
			   int &ei, int &ej, int &ek, int &et, int &level) {


  int idx = GetCoords(x,y,z,t,ei, ej, ek, et, level); 

  return idx; 

}

/////////////////////////////////////////////////////////////////
//
//   Map a cell from one level to the cell(s) in the over level 
// 
bool LatticeAMR::MapCells(int fromI, int fromJ, int fromK, int fromT, 
		 int fromLevel, int toLevel, int& toImin, int& toImax, 
		 int& toJmin, int& toJmax, int& toKmin, int& toKmax, 
		 int& toTmin, int& toTmax)
{
  float minX = xmin[fromLevel]+fromI*xlength[fromLevel]; 
  float maxX = xmin[fromLevel]+(fromI+1)*xlength[fromLevel]; 
  float minY = ymin[fromLevel]+fromJ*ylength[fromLevel]; 
  float maxY = ymin[fromLevel]+(fromJ+1)*ylength[fromLevel]; 
  float minZ = zmin[fromLevel]+fromK*zlength[fromLevel]; 
  float maxZ = zmin[fromLevel]+(fromK+1)*zlength[fromLevel]; 
  float minT = tmin[fromLevel]+fromT*tlength[fromLevel]; 
  float maxT = tmin[fromLevel]+(fromT+1)*tlength[fromLevel]; 

  toImin = (int)((minX - xmin[toLevel])/(float)xlength[toLevel]); 
  if (toImin<0 || toImin>idim[toLevel]) return false; 
  toImax = (int)((maxX - xmin[toLevel])/(float)xlength[toLevel]); 
  if (toImax<0 || toImax>idim[toLevel]) return false; 
  toJmin = (int)((minY - ymin[toLevel])/(float)ylength[toLevel]); 
  if (toJmin<0 || toJmin>jdim[toLevel]) return false; 
  toJmax = (int)((maxY - ymin[toLevel])/(float)ylength[toLevel]); 
  if (toJmax<0 || toJmax>jdim[toLevel]) return false; 
  toKmin = (int)((minZ - zmin[toLevel])/(float)zlength[toLevel]); 
  if (toKmin<0 || toKmin>kdim[toLevel]) return false; 
  toKmax = (int)((maxZ - zmin[toLevel])/(float)zlength[toLevel]); 
  if (toKmax<0 || toKmax>kdim[toLevel]) return false; 
  toTmin = (int)((minT - tmin[toLevel])/(float)tlength[toLevel]); 
  if (toTmin<0 || toTmin>ldim[toLevel]) return false; 
  toTmax = (int)((maxT - tmin[toLevel])/(float)tlength[toLevel]); 
  if (toTmax<0 || toTmax>ldim[toLevel]) return false; 
  return true; 
}

//////////////////////////////////////////////////////////////////////
//
//    To see if data from blocks of a higher LOD can be 
//    merged into the block (i,j,k,t) at this level
//
//    The purpose of merger is to reduce the number of blocks in the 
//    data set for efficiency 
//    
//
bool LatticeAMR::Mergeable(int i, int j, int k, int t, int level, 
			   int& mergeLevel)
{
  int toImin, toImax, toJmin, toJmax, toKmin, toKmax, toTmin, toTmax; 

  int idx = t*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 

  // Blocks from which level (as deep as possible) can be merged to this level? 
  mergeLevel = level; 
  if (has_data[level][idx]==true) { //the block has data. no merging needed 
    has_data_from_merger[level][idx] = level; //default: its own level 
    return true; 
  }
  else // check finer levels
    if (level+1 >= num_levels) // this is the finest level 
      return false; // this means in this region no data is available 
    else {
      // Can the space of this block be filled by blocks in the next level? 
      // if yes, what are the indices of those next level blocks 
      bool mappable = MapCells(i,j,k,t,level,level+1, toImin, toImax, toJmin, toJmax, 
		toKmin, toKmax, toTmin, toTmax); 
      if (mappable == false) {
	return false; // out of next level's bound 
      }

      bool first = true; 
      for (int tFor = toTmin; tFor<toTmax; tFor++) 
	for (int kFor = toKmin; kFor<toKmax; kFor++)
	  for (int jFor = toJmin; jFor<toJmax; jFor++)
	    for (int iFor = toImin; iFor<toImax; iFor++) {
	      // recursively call the next level 
	      int mLevel; 
	      if (Mergeable(iFor, jFor, kFor, tFor, level+1, 
			    mLevel)== false) {
		return(false); 
	      }
	      else {
		if (first == true)  {
		  mergeLevel = mLevel; 
		  first = false; 
		} 
		else  // all the merge levels from children have to be the same 
		  if (mLevel != mergeLevel) {
		    return(false); 
		  }
	      }
	    }
      // mergeable 
      // reset the children's record 
      for (int tFor = toTmin; tFor<toTmax; tFor++) 
	for (int iFor = toImin; iFor<toImax; iFor++)
	  for (int jFor = toJmin; jFor<toJmax; jFor++)
	    for (int kFor = toKmin; kFor<toKmax; kFor++){
	      int idx2 = tFor*idim[level+1]*jdim[level+1]*kdim[level+1]+kFor*idim[level+1]*jdim[level+1]+jFor*idim[level+1]+iFor; 
	      has_data_from_merger[level+1][idx2] = -1; 
	    }
      has_data_from_merger[level][idx]=mergeLevel; 
      return(true); 
    }
}
////////////////////////////////////////////////////////////
//
//  Merge smaller blocks of the same resolutions into a larger 
//  block whenever possible 
//

void LatticeAMR::MergeBlocks()
{
  int l = 0; // 0-th level;  coarsest resolution 
             // other levels will be done through recursive calls 
  for (int tFor=0; tFor<ldim[l]; tFor++)
    for (int kFor=0; kFor<kdim[l]; kFor++)
      for (int jFor=0; jFor<jdim[l]; jFor++)
	for (int iFor=0; iFor<idim[l]; iFor++) {
	  int mergeLevel; 
	  bool mg = Mergeable(iFor, jFor, kFor, tFor, l, mergeLevel); 
	  if (mg == true) 
	    printf("%d %d %d %d got level %d merger\n", iFor, jFor, kFor, tFor, 
		   mergeLevel); 
	}

  int cnt = 0; 
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data_from_merger[i][j]!=-1) cnt++; 
  printf(" cnt = %d \n", cnt); 

}

//----------------------------------------------------------------------------

void LatticeAMR::InitSeedLists() {

  seedlists = new list<VECTOR4>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

void LatticeAMR::ResetSeedLists() {

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}

void LatticeAMR::ResetSeedLists(int i) {

  seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

bool LatticeAMR::InsertSeed(int i, int j, int k, int t, int level, VECTOR4 p) {

  int rank = GetRank(i,j,k, t, level); 

  if (rank ==-1) return(false); 
  else {
    seedlists[rank].push_back(p); 
    return(true); 
  }
}


bool LatticeAMR::InsertSeed(int i, VECTOR4 p) {

  if (i<0 || i>=npart) return(false); 
  else {
    seedlists[i].push_back(p); 
    return(true); 
  }
}


//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void LatticeAMR::RoundRobin_proc() {

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
// assign the partitions to the processors in contiguous ranges
//
void LatticeAMR::ContigRange_proc() {

  int proc; // process number
  int n = 0; // index into a row of the process table, local block number
  int i; // partition rank
  int groupsize = npart / nproc; // number of blocks per process
                           // plus remainder in last process

  for (i = 0; i < npart; i++) {

    proc = i / groupsize;
    if (proc >= nproc) {// remainder: clamp extra blocks to last process
      proc = nproc - 1;
      n++;
    }
    else
      n = i % groupsize;

    // assign the process number
    part->parts[i].Proc = proc; 

    // add entry into the process table
    part->proc_parts[proc][n] = i;
    part->proc_nparts[proc] = n + 1;

  }

}
//---------------------------------------------------------------------------
//
// return the process that owns the partition rank
//
int LatticeAMR::GetProc(int rank) {

  return(part->parts[rank].Proc); 

}

//---------------------------------------------------------------------------
//
// look up the lattice[i,j,k] element at the given level 
//  for the processor of the subdomain
//
int LatticeAMR::GetProc(int i, int j, int k, int t, int level) 
{
  if (level <0 || level >=num_levels) return (-1); 

  if (i < 0 || i >= idim[level])
    return(-1); 
  else if (j < 0 || j >= jdim[level])
    return(-1); 
  else if (k < 0 || k >= kdim[level])
    return(-1); 
  else if (t < 0 || t >= ldim[level])
    return(-1); 

  int rank = GetRank(i,j,k,t,level); 
  if (rank == -1) return (-1); 
  else 
    return (part->parts[rank].Proc); 
}

//----------------------------------------------------------------------------
//
// query the paritions that are assigned to processor 'proc' 
//
void LatticeAMR::GetPartitions(int proc, int**p_list, int& num) {

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
//---------------------------------------------------------------------------
//
// get the number of voxels per block (eg 16x16x16)
// assumed constant for all blocks
//
void LatticeAMR::GetBlockDims(int *dims) {

  int i;

  for (i = 0; i < 3; i++)
    dims[i] = block_dims[i];

}
//----------------------------------------------------------------------------
//
// query the partitions that are assigned to processor 'proc' 
// p_list must be allocated large enough prior to calling
//
void LatticeAMR::GetMyPartitions(int proc, int* p_list) {

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
int LatticeAMR::GetMyNumPartitions(int proc) {

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
int LatticeAMR::GetMyNumPartitions() {

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
// gets local subvolume bounds
//
// block: local block number (0-nblocks)
// min_s, max_s: (output) spatial min and max bounds
// min_t, max_t: (output) temporal min and max bounds
//
void LatticeAMR::GetVB(int block, float *min_s, float *max_s, 
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
void LatticeAMR::GetGlobalVB(int part, float *min_s, float *max_s, 
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
// returns neighbor number of neighbor containing point
//
// returns -1 if point is not in one of the neighbors
//
// this can happen in two cases:
// the point remained in the current block, or
// the point left the overall data boundary
//
int LatticeAMR::GetNeighbor(int block, float x, float y, float z, float t) {

  int n; // neighbor number
  int r; // rank of partition containing the point

  r = GetRank(x, y, z, t);

  // check if the point never left the current block
  if (r == block_ranks[block])
    return -1;

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
// the neighbor of an edge gets rank -1
//
void LatticeAMR::GetNeighborRanks(int block) {

  // levels coarser than mine are handled slightly separately
  GetCoarseNeighborRanks(block);
  GetFineNeighborRanks(block);

}
//---------------------------------------------------------------------------
//
// gets ranks of all neighbors at a coarser level than mine
// the neighbor of an edge gets rank -1
//
void LatticeAMR::GetCoarseNeighborRanks(int block) {

  float x, y, z, t; // point inside of partition
  float x0, y0, z0, t0; // reference point for point calculation
  float dx, dy, dz, dt; // delta x,y,z,t
  int i, j, k, l; // coords along an end face
  int l_max; // maximum value l (index in time dimension)
  int level; // level number
  int nr; // neighbor rank
  int old_nr[81]; // previous neighbor ranks
  int num_onr; // number of previous neighbor ranks to check
  int myrank = block_ranks[block];
  int n;

  // get minimum corner point and initial block size
  x0 = vb_list[myrank].xmin;
  y0 = vb_list[myrank].ymin;
  z0 = vb_list[myrank].zmin;
  t0 = vb_list[myrank].tmin;
  dx = vb_list[myrank].xmax - vb_list[myrank].xmin;
  dy = vb_list[myrank].ymax - vb_list[myrank].ymin;
  dz = vb_list[myrank].zmax - vb_list[myrank].zmin;
  dt = vb_list[myrank].tmax - vb_list[myrank].tmin;

  // get level number of block
  float cx = x0 + dx / 2.0f;
  float cy = y0 + dy / 2.0f;
  float cz = z0 + dz / 2.0f;
  float ct = t0 + dt / 2.0f;
  int block_level = GetFinestLevel(cx, cy, cz, ct);
  assert(block_level >= 0);

  // for levels coarser than the block level
  for (level = 0; level < block_level; level++) {

    // set iteration params
    l_max = dt > 0.0 ? 2 : 0; // extend

    // -x end face
    num_onr = 0;
    x = vb_list[myrank].xmin - 0.5 * dx;
    for (j = -1; j < 2; j++) { // y direction
      y = y0 + (j + 0.5) * dy;
      for (k = -1; k < 2; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 && 
	      vb_list[nr].xmax == vb_list[myrank].xmin) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // z direction
    } // y direction

    // +x end face
    num_onr = 0;
    x = vb_list[myrank].xmax + 0.5 * dx;
    for (j = -1; j < 2; j++) { // y direction
      y = y0 + (j + 0.5) * dy;
      for (k = -1; k < 2; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 &&
	      vb_list[nr].xmin == vb_list[myrank].xmax) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // z direction
    } // y direction

    // -y end face
    num_onr = 0;
    y = vb_list[myrank].ymin - 0.5 * dy;
    for (i = -1; i < 2; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (k = -1; k < 2; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 &&
	      vb_list[nr].ymax == vb_list[myrank].ymin && 
	      vb_list[nr].xmax != vb_list[myrank].xmin && 
	      vb_list[nr].xmin != vb_list[myrank].xmax) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // z direction
    } // x direction

    // +y end face
    num_onr = 0;
    y = vb_list[myrank].ymax + 0.5 * dy;
    for (i = -1; i < 2; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (k = -1; k < 2; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 &&
	      vb_list[nr].ymin == vb_list[myrank].ymax && 
	      vb_list[nr].xmax != vb_list[myrank].xmin && 
	      vb_list[nr].xmin != vb_list[myrank].xmax) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // z direction
    } // x direction

    // -z end face
    num_onr = 0;
    z = vb_list[myrank].zmin - 0.5 * dz;
    for (i = -1; i < 2; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (j = -1; j < 2; j++) { // y direction
	y = y0 + (j + 0.5) * dy;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 &&
	      vb_list[nr].zmax == vb_list[myrank].zmin && 
	      vb_list[nr].xmax != vb_list[myrank].xmin && 
	      vb_list[nr].xmin != vb_list[myrank].xmax &&
	      vb_list[nr].ymax != vb_list[myrank].ymin && 
	      vb_list[nr].ymin != vb_list[myrank].ymax) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // y direction
    } // x direction

    // +z end face
    num_onr = 0;
    z = vb_list[myrank].zmax + 0.5 * dz;
    for (i = -1; i < 2; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (j = -1; j < 2; j++) { // y direction
	y = y0 + (j + 0.5) * dy;
	for (l = -1; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0 &&
	      vb_list[nr].zmin == vb_list[myrank].zmax && 
	      vb_list[nr].xmax != vb_list[myrank].xmin && 
	      vb_list[nr].xmin != vb_list[myrank].xmax &&
	      vb_list[nr].ymax != vb_list[myrank].ymin && 
	      vb_list[nr].ymin != vb_list[myrank].ymax) {
	    for (n = 0; n < num_onr; n++) {
	      if (nr == old_nr[n])
		break;
	    }
	    if (n == num_onr) {
	      AddNeighbor(block, nr);
	      old_nr[num_onr++] = nr;
	    }
	  }
	} // t direction
      } // y direction
    } // x direction

    // neighbors in the time dimension exist only if the block has some
    // "thickness" in the time dimension (dt > 0)
    if (dt > 0.0) {

      // -t end face
      num_onr = 0;
      t = vb_list[myrank].tmin - 0.5 * dt;
      for (i = -1; i < 2; i++) { // x direction along face
	x = x0 + (i + 0.5) * dx;
	for (j = -1; j < 2; j++) { // y direction
	  y = y0 + (j + 0.5) * dy;
	  for (k = -1; k < 2; k++) { // z direction
	    z = z0 + (k + 0.5) * dz;
	    nr = GetRank(x, y, z, t, level);
	    if (nr >= 0 &&
		vb_list[nr].tmax == vb_list[myrank].tmin && 
		vb_list[nr].xmax != vb_list[myrank].xmin && 
		vb_list[nr].xmin != vb_list[myrank].xmax &&
		vb_list[nr].ymax != vb_list[myrank].ymin && 
		vb_list[nr].ymin != vb_list[myrank].ymax &&
		vb_list[nr].zmax != vb_list[myrank].zmin && 
		vb_list[nr].zmin != vb_list[myrank].zmax) {
	      for (n = 0; n < num_onr; n++) {
		if (nr == old_nr[n])
		  break;
	      }
	      if (n == num_onr) {
		AddNeighbor(block, nr);
		old_nr[num_onr++] = nr;
	      }
	    }
	  } // z direction
	} // y direction
      } // x direction
    
      // +t end face
      num_onr = 0;
      t = vb_list[myrank].tmax + 0.5 * dt;
      for (i = -1; i < 2; i++) { // x direction along face
	x = x0 + (i + 0.5) * dx;
	for (j = -1; j < 2; j++) { // y direction
	  y = y0 + (j + 0.5) * dy;
	  for (k = -1; k < 2; k++) { // z direction
	    z = z0 + (k + 0.5) * dz;
	    nr = GetRank(x, y, z, t, level);
	    if (nr >= 0 &&
		vb_list[nr].tmin == vb_list[myrank].tmax && 
		vb_list[nr].xmax != vb_list[myrank].xmin && 
		vb_list[nr].xmin != vb_list[myrank].xmax &&
		vb_list[nr].ymax != vb_list[myrank].ymin && 
		vb_list[nr].ymin != vb_list[myrank].ymax &&
		vb_list[nr].zmax != vb_list[myrank].zmin && 
		vb_list[nr].zmin != vb_list[myrank].zmax) {
	      for (n = 0; n < num_onr; n++) {
		if (nr == old_nr[n])
		  break;
	      }
	      if (n == num_onr) {
		AddNeighbor(block, nr);
		old_nr[num_onr++] = nr;
	      }
	    }
	  } // z direction
	} // y direction
      } // x direction

    } // dt > 0

  } // for all coarser levels

}
//---------------------------------------------------------------------------
//
// gets ranks of all neighbors at the same or finer level than mine
// the neighbor of an edge gets rank -1
//
void LatticeAMR::GetFineNeighborRanks(int block) {

  float x, y, z, t; // point inside of partition
  float x0, y0, z0, t0; // reference point for point calculation
  float dx, dy, dz, dt; // delta x,y,z,t
  int i, j, k, l; // coords along an end face
  int i_min, j_min, k_min, l_min; // minimum value of i, j, k, l
  int i_max, j_max, k_max, l_max; // maximum value of i, j, k, l
  int level; // level number
  int n = 1; // number of (smaller) partitions in the end face
  int nr; // neighbor rank
  int myrank = block_ranks[block];

  // get minimum corner point and initial block size
  x0 = vb_list[myrank].xmin;
  y0 = vb_list[myrank].ymin;
  z0 = vb_list[myrank].zmin;
  t0 = vb_list[myrank].tmin;
  dx = vb_list[myrank].xmax - vb_list[myrank].xmin;
  dy = vb_list[myrank].ymax - vb_list[myrank].ymin;
  dz = vb_list[myrank].zmax - vb_list[myrank].zmin;
  dt = vb_list[myrank].tmax - vb_list[myrank].tmin;

  // get level number of block
  float cx = x0 + dx / 2.0f;
  float cy = y0 + dy / 2.0f;
  float cz = z0 + dz / 2.0f;
  float ct = t0 + dt / 2.0f;
  int block_level = GetFinestLevel(cx, cy, cz, ct);
  assert(block_level >= 0);

  // for levels starting at the the block level and finer
  for (level = block_level; level < num_levels; level++) {

    // set iteration params
    j_min = k_min = l_min = -1;
    j_max = dy > 0.0 ? n + 1 : 0; // extend
    k_max = dz > 0.0 ? n + 1 : 0; // extend
    l_max = dt > 0.0 ? n + 1 : 0; // extend

    // -x end face
    x = vb_list[myrank].xmin - 0.5 * dx;
    for (j = j_min; j < j_max; j++) { // y direction
      y = y0 + (j + 0.5) * dy;
      for (k = k_min; k < k_max; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // z direction
    } // y direction

    // +x end face
    x = vb_list[myrank].xmax + 0.5 * dx;
    for (j = j_min; j < j_max; j++) { // y direction
      y = y0 + (j + 0.5) * dy;
      for (k = k_min; k < k_max; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // z direction
    } // y direction

    // -y end face
    i_min = 0;
    i_max = n; // same size
    y = vb_list[myrank].ymin - 0.5 * dy;
    for (i = i_min; i < i_max; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (k = k_min; k < k_max; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // z direction
    } // x direction

    // +y end face
    y = vb_list[myrank].ymax + 0.5 * dy;
    for (i = i_min; i < i_max; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (k = k_min; k < k_max; k++) { // z direction
	z = z0 + (k + 0.5) * dz;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // z direction
    } // x direction

    // -z end face
    j_min = 0;
    j_max = n; // same size
    z = vb_list[myrank].zmin - 0.5 * dz;
    for (i = i_min; i < i_max; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (j = j_min; j < j_max; j++) { // y direction
	y = y0 + (j + 0.5) * dy;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // y direction
    } // x direction

    // +z end face
    z = vb_list[myrank].zmax + 0.5 * dz;
    for (i = i_min; i < i_max; i++) { // x direction along face
      x = x0 + (i + 0.5) * dx;
      for (j = j_min; j < j_max; j++) { // y direction
	y = y0 + (j + 0.5) * dy;
	for (l = l_min; l < l_max; l++) { // t direction
	  t = t0 + (l + 0.5) * dt;
	  nr = GetRank(x, y, z, t, level);
	  if (nr >= 0)
	    AddNeighbor(block, nr);
	} // t direction
      } // y direction
    } // x direction

    // neighbors in the time dimension exist only if the block has some
    // "thickness" in the time dimension (dt > 0)
    if (dt > 0.0) {

      // -t end face
      k_min = 0;
      k_max = n; // same size
      t = vb_list[myrank].tmin - 0.5 * dt;
      for (i = i_min; i < i_max; i++) { // x direction along face
	x = x0 + (i + 0.5) * dx;
	for (j = j_min; j < j_max; j++) { // y direction
	  y = y0 + (j + 0.5) * dy;
	  for (k = k_min; k < k_max; k++) { // z direction
	    z = z0 + (k + 0.5) * dz;
	    nr = GetRank(x, y, z, t, level);
	    if (nr >= 0)
	      AddNeighbor(block, nr);
	  } // z direction
	} // y direction
      } // x direction
    
      // +t end face
      t = vb_list[myrank].tmax + 0.5 * dt;
      for (i = i_min; i < i_max; i++) { // x direction along face
	x = x0 + (i + 0.5) * dx;
	for (j = j_min; j < j_max; j++) { // y direction
	  y = y0 + (j + 0.5) * dy;
	  for (k = k_min; k < k_max; k++) { // z direction
	    z = z0 + (k + 0.5) * dz;
	    nr = GetRank(x, y, z, t, level);
	    if (nr >= 0)
	      AddNeighbor(block, nr);
	  } // z direction
	} // y direction
      } // x direction

    } // dt > 0

    n *= 2;
    dx *= 0.5f;
    dy *= 0.5f;
    dz *= 0.5f;
    dt *= 0.5f;

  } // for all levels

  // don't forget to include myself
  AddNeighbor(block, myrank);

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
void LatticeAMR::AddNeighbor(int myblock, int neighrank) {

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
// gets the data pointer
//
float** LatticeAMR::GetData(int block) {
  return GetDataPtr(block_ranks[block]);
}
//
// sets the request status
//
void LatticeAMR::SetReq(int block) {
  part->SetReq(block_ranks[block]);
}
//
// clears the request status
//
void LatticeAMR::ClearReq(int block) {
  part->ClearReq(block_ranks[block]);
}
//
// gets the request status
//
int LatticeAMR::GetReq(int block) {
  return part->GetReq(block_ranks[block]); 
}
//
// sets the load status
//
void LatticeAMR::SetLoad(int block) {
  part->SetLoad(block_ranks[block]); 
}
//
// clears the load status
//
void LatticeAMR::ClearLoad(int block) { 
  part->ClearLoad(block_ranks[block]);
}
//
// gets the load status
//
int LatticeAMR::GetLoad(int block) {
  return part->GetLoad(block_ranks[block]); 
}
//
// sets the computed status
//
void LatticeAMR::SetComp(int block, int iter) {
  part->SetComp(block_ranks[block], iter); 
}
//
// clears the computed status
//
void LatticeAMR::ClearComp(int block) {
  part->ClearComp(block_ranks[block]); 
}
//
// gets the computed status
//
int LatticeAMR::GetComp(int block, int iter) {
  return part->GetComp(block_ranks[block], iter); 
}
//
// posts a point for sending
//
void LatticeAMR::PostPoint(int block, VECTOR4 p) {
  int neighbor = GetNeighbor(block, p[0], p[1], p[2], p[3]);
  // only post points that move out of the current block and
  // remain inside the overall domain boundary
  if (neighbor >= 0)
    part->PostPoint(block_ranks[block], p, neighbor); 
}
//
// prints the posted points
//
void LatticeAMR::PrintPost(int block) {
  part->PrintPost(block_ranks[block]); 
}
//
// prints the received points
//
void LatticeAMR::PrintRecv(int block) { 
  part->PrintRecv(block_ranks[block]); 
}
//
// exchanges points with all neighbors
// returns total number of points received by this process
//
int LatticeAMR::ExchangeNeighbors(VECTOR4 **seeds, int *size_seeds) { 
  int n;
  comm_time = MPI_Wtime();
  n = part->ExchangeNeighbors(neighbor_ranks, seeds, size_seeds);
  comm_time = MPI_Wtime() - comm_time;
  return n;
}
//---------------------------------------------------------------------------
