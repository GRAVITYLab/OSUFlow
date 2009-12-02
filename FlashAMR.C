//------------------------------------------------------------------------------
//
// static and time-varying FLASH AMR data
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

#include "flashhdf5_float.h"
#include "FlashAMR.h" 

#ifdef _MPI
#include <mpi.h>
#endif


//-----------------------------------------------------------------------
//
// utility functions, not part of any class
//
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//
// int_sort_list (integer version)
//
// sorts a list containing duplicate entries into distinct elements
// only works for nonnegative values
//
// list: input / output list
// cnt: input / output count of elements
// cnt will be updated to the number of distinct elements
// the first cnt distinct elements will be placed at the start of list
//
void int_sort_list(int* list, int& cnt)  {

  int* tmp_list = new int[cnt]; 

  for (int i = 0; i < cnt; i++)
    tmp_list[i] = list[i]; 


  for (int i = 0; i< cnt - 1; i++) {
    int largest = -1; 
    int idx = i; 
    for (int j = i; j < cnt; j++) {
      if (tmp_list[j] > largest) {
	idx = j;
	largest = tmp_list[j];
      }
    }
    int tmp = (int)tmp_list[i]; 
    tmp_list[i] = tmp_list[idx]; 
    tmp_list[idx] = tmp; 
  }

  list[0] = tmp_list[cnt - 1]; 
  int count = 1; 

  for (int i = cnt - 2; i >= 0; i--) {
    if (tmp_list[i] > list[count-1]) {
      list[count] = tmp_list[i]; 
      count++; 
    }
  }

  cnt = count; 

}
//-----------------------------------------------------------------------
//
// float_sort_list (float version)
//
// sorts a list containing duplicate entries into distinct elements
// only works for nonnegative values
//
// list: input / output list
// cnt: input / output count of elements
// cnt will be updated to the number of distinct elements
// the first cnt distinct elements will be placed at the start of list
//
void float_sort_list(float* list, int& cnt)  {

  float* tmp_list = new float[cnt]; 

  for (int i = 0; i < cnt; i++)
    tmp_list[i] = list[i]; 


  for (int i = 0; i< cnt - 1; i++) {
    float largest = -1; 
    int idx = i; 
    for (int j = i; j < cnt; j++) {
      if (tmp_list[j] > largest) {
	idx = j;
	largest = tmp_list[j];
      }
    }
    int tmp = (int)tmp_list[i]; 
    tmp_list[i] = tmp_list[idx]; 
    tmp_list[idx] = tmp; 
  }

  list[0] = tmp_list[cnt - 1]; 
  int count = 1; 

  for (int i = cnt - 2; i >= 0; i--) {
    if (tmp_list[i] > list[count-1]) {
      list[count] = tmp_list[i]; 
      count++; 
    }
  }

  cnt = count; 

}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//
// FlashAMR functions, static (single timestep) AMR
//
//-----------------------------------------------------------------------

#ifdef _MPI

//-----------------------------------------------------------------------
//
// ParallelLoadHDF5MetaData
//
// parallel reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
//
// fname: file name
// min, max: (output) global min, max corners
// comm: MPI communicator
//
int FlashAMR::ParallelLoadHDF5MetaData(char* fname, float min[3], float max[3],
			       MPI_Comm comm) {

  fdf = new FlashHDFFile(fname, comm); // flash file object
  LoadHDF5MetaData(min, max);

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// SerialLoadHDF5MetaData
//
// serial reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
// currently exists in a serial version only
//
// fname: file name
// min, max: (output) global min, max corners
//
void FlashAMR::SerialLoadHDF5MetaData(char* fname, float min[3], float max[3]) {

  fdf = new FlashHDFFile(fname); // flash file object
  LoadHDF5MetaData(min, max);

}
//-----------------------------------------------------------------------
//
// LoadHDF5MetaData
//
// reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
// currently exists in a serial version only
//
// min, max: (output) global min, max corners
//
void FlashAMR::LoadHDF5MetaData(float min[3], float max[3]) {

  float *bounds; // bounding box: min, max corners
  int *init; // initial flag for each level
  int *lmap; // map of levels in file to position in block_level list
  int *local_block_level; // local version of block_level
  int level; // current level
  int i, j;

  // number of blocks and block dimensions
  nb = fdf->GetNumberOfBlocks(); // includes non-leaf blocks
  file_nb = nb; // save the total number of file block (leaf + nonleaf)
  fdf->GetCellDimensions(block_dims);

  // allocate memory for leaf and non-leaf blocks
  block_level = new int[nb]; // level that each block belongs to 
  local_block_level = new int[nb]; // local copy of block_level
  block_center = new float[nb * 3]; // center x y and z of each block 
  block_length = new float[nb * 3]; // physical length (x, y, z) of each block 
  block_minB = new float[nb * 3]; // min corner (x, y, z) of each blk
  block_maxB = new float[nb * 3]; // max corner (x, y, z) of each blk
  bounds = new float[nb * 6]; // bounding box of each blk

  // center, length, level, extents of each block
  nb = fdf->GetFloatVecBlocks((char *)"coordinates", 0, file_nb, 3, 1, 
			      block_center);
  fdf->GetFloatVecBlocks((char *)"block size", 0, file_nb, 3, 1, block_length);
  fdf->GetIntVecBlocks((char *)"refine level", 0, file_nb, 1, 1, block_level);
  fdf->GetFloatMatBlocks((char *)"bounding box", 0, file_nb, 3, 2, 1, bounds);

  for (i = 0; i < nb; i++) {
    local_block_level[i] = block_level[i];
    for (j = 0; j < 3; j++) {
      block_minB[3 * i + j] = bounds[2 * (3 * i + j)];
      block_maxB[3 * i + j] = bounds[2 * (3 * i + j) + 1];
    }
  }

  // overall bounds of entire dataset
  for (i = 0; i < nb; i++) {
    if (i == 0)
      for (j = 0; j < 3; j++) {
	min[j] = block_minB[j];
	max[j] = block_maxB[j];
      }
    else
      for (j = 0; j < 3; j++) {
	if (block_minB[3 * i + j] < min[j]) 
	  min[j] = block_minB[3 * i + j]; 
	if (block_maxB[3 * i + j] > max[j]) 
	  max[j] = block_maxB[3 * i + j]; 
      }
  }

  // list of distinct levels
  num_levels = nb;
  int_sort_list(local_block_level, num_levels);
  max_level_value = local_block_level[num_levels - 1]; // needed later for
                                                // level map in this time step
                                                //  and in merged time steps

  // extents and block size for each level
  // block size is physical size (length)
  level_minB = new float[num_levels * 3];
  level_maxB = new float[num_levels * 3]; 
  block_xsize_inLevel = new float[num_levels];
  block_ysize_inLevel = new float[num_levels]; 
  block_zsize_inLevel = new float[num_levels]; 
  init = new int[num_levels];
  for (i = 0; i < num_levels; i++)
    init[i] = 1;
  lmap = new int[max_level_value + 1];
  for (i = 0; i < num_levels; i++)
    lmap[local_block_level[i]] = i;

  for (i = 0; i < nb; i++) {

    // level = position in block_level list of 
    // refinement level in file
    level = lmap[block_level[i]];

    // level block size
    block_xsize_inLevel[level] = block_length[3 * i]; 
    block_ysize_inLevel[level] = block_length[3 * i + 1]; 
    block_zsize_inLevel[level] = block_length[3 * i + 2]; 

    // level extents
    if (init[level]) {
      for (j = 0; j < 3; j++) {
	level_minB[3 * level + j] = block_minB[3 * i + j];
	level_maxB[3 * level + j] = block_maxB[3 * i + j];
      }
      init[level] = 0;
    }
    else {
      for (j = 0; j < 3; j++) {
	if (block_minB[3 * i + j] < level_minB[3 * level + j]) 
	  level_minB[3 * level + j] = block_minB[3 * i + j]; 
	if (block_maxB[3 * i + j] > level_maxB[3 * level + j]) 
	  level_maxB[3 * level + j] = block_maxB[3 * i + j]; 
      }
    }

  }

  delete [] init;
  delete [] bounds;

#ifdef DEBUG
  if (myproc == 0) {
    fprintf(stderr,"Number of leaf + nonleaf blocks = %d\n", file_nb); 
    fprintf(stderr,"Number of leaf blocks = %d\n", nb); 
    fprintf(stderr,"Each block is %dx%dx%d cells\n",
	    block_dims[0], block_dims[1], block_dims[2]);
  }
#endif

}
//-----------------------------------------------------------------------
//
// LoadHDF5Data
//
// reader for HDF5 Flash vector data from 
// reads start_block to end_block (inclusive)
//
// start_block, end_block: contiguous range of blocks to read
// vx, vy, vz: names of three velocity components in the dataset
//
int FlashAMR::LoadHDF5Data(int start_block, int end_block,
				 char * vx, char * vy, char *vz) {

  int size = block_dims[0] * block_dims[1] * block_dims[2]; // total voxels
  float *comp; // one component of vector data
  int i, j;

  int nblocks = end_block - start_block + 1;

  // allocate memory for leaf + nonleaf blocks
  vectors = new float*[nb]; // the vector data
  comp = new float[file_nb * size]; // don't know how much data will need to be
                                    //  read in order to get nblocks leaf blocks
                                    // todo: find a tighter bound for comp?

  for (i = 0; i < nb; i++)
    vectors[i] = new float[size * 3];

  // x component
  fdf->GetFloatVolBlocks(vx, start_block, nblocks, block_dims[0], 
			     block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {

    vectors[start_block + i] = new float[size * 3]; // large enough for all
                                                    // three components
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 0] = comp[i * size + j]; 
  }

  // y component
  fdf->GetFloatVolBlocks(vy, start_block, nblocks, block_dims[0], 
			     block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 1] = comp[i * size + j]; 
  }

  // z component
  fdf->GetFloatVolBlocks(vz, start_block, nblocks, block_dims[0], 
			     block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 2] = comp[i * size + j]; 
  }

  delete [] comp; 

  return(1); 

}
//-----------------------------------------------------------------------
//
// GetLevelBlockSize
//
// physical size of a block in a level
//
void FlashAMR::GetLevelBlockSize(int level, float size[3]) {

  size[0] = block_xsize_inLevel[level]; 
  size[1] = block_ysize_inLevel[level]; 
  size[2] = block_zsize_inLevel[level]; 

}
//-----------------------------------------------------------------------
//
// GetLevelBounds
//
// extents of all the blocks in a level
//
void FlashAMR::GetLevelBounds(int level, float minB[3], float maxB[3]){

  minB[0] = level_minB[level*3]; 
  minB[1] = level_minB[level*3+1]; 
  minB[2] = level_minB[level*3+2]; 

  maxB[0] = level_maxB[level*3]; 
  maxB[1] = level_maxB[level*3+1]; 
  maxB[2] = level_maxB[level*3+2]; 

}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//
// TimeVaryingFlashAMR functions, time-varying (multiple timesteps) AMR
//
//-----------------------------------------------------------------------

#ifdef _MPI

//-----------------------------------------------------------------------
//
// LoadMetaData
//
// reads metadata collectively for a time-varying flash AMR dataset
//
// fname: file containing list of timestep files
// min, max: (output) global min, max extent over entire dataset (space+time)
//
void TimeVaryingFlashAMR::LoadMetaData(char* fname, float min[3], float max[3],
				      MPI_Comm comm) {

  FILE *fIn; 
  char filename[300]; 
  float pmin[3], pmax[3]; 
  int i, j;

  // get the number of timesteps, allocate amr objects
  assert((fIn = fopen(fname, "r")) != NULL);
  fscanf(fIn, "%d", &num_timesteps); 
  amr_list = new FlashAMR*[num_timesteps]; 

  // load the metadata (extents, level info)
  for (i = 0; i < num_timesteps; i++) {

    // read
    fscanf(fIn, "%s", filename); 

#ifdef DEBUG
    if (myproc == 0)
      fprintf(stderr,"Reading metadata from %s ...\n\n", filename); 
#endif

    amr_list[i] = new FlashAMR(myproc);
//     amr_list[i]->ParallelLoadHDF5MetaData(filename, pmin, pmax, comm); 
    amr_list[i]->SerialLoadHDF5MetaData(filename, pmin, pmax); 

    // compute data extents over all the timesteps
    if (i == 0) {
      for (j = 0; j < 3; j++) {
	min[j] = pmin[j];
	max[j] = pmax[j];
      }
    }
    else {
      for (j = 0; j < 3; j++) {
	if (pmin[j] < min[j])
	  min[j] = pmin[j]; 
	if (pmax[j] > max[j])
	  max[j] = pmax[j]; 
      }
    }

  }

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"Overall volume bounds: min=[%.4e %.4e %.4e]\
   max=[%.4e %.4e %.4e]\n\n", min[0], min[1], min[2], max[0], max[1], max[2]);
#endif
 
  amr_list[0]->GetDims(block_dims); // assume all time steps are the same 
  level_mapping = new int*[num_timesteps]; 

  // reconcile metadata over all timesteps
  MatchTimestepLevels(); 

}
//-----------------------------------------------------------------------

#else

//-----------------------------------------------------------------------
//
// LoadMetaData
//
// reads metadata independently for a time-varying flash AMR dataset
//
// fname: file containing list of timestep files
// min, max: (output) global min, max extent over entire dataset (space+time)
//
void TimeVaryingFlashAMR::LoadMetaData(char* fname, float min[3], float max[3]) {

  FILE *fIn; 
  char filename[300]; 
  float pmin[3], pmax[3]; 
  int i, j;

  // get the number of timesteps, allocate amr objects
  assert((fIn = fopen(fname, "r")) != NULL);
  fscanf(fIn, "%d", &num_timesteps); 
  amr_list = new FlashAMR*[num_timesteps]; 

  // load the metadata (extents, level info)
  for (i = 0; i < num_timesteps; i++) {

    // read
    fscanf(fIn, "%s", filename); 

#ifdef DEBUG
    if (myproc == 0)
      fprintf(stderr,"Reading metadata from %s ...\n\n", filename); 
#endif

    amr_list[i] = new FlashAMR(myproc);
    amr_list[i]->SerialLoadHDF5MetaData(filename, pmin, pmax); 

    // compute data extents over all the timesteps
    if (i == 0) {
      for (j = 0; j < 3; j++) {
	min[j] = pmin[j];
	max[j] = pmax[j];
      }
    }
    else {
      for (j = 0; j < 3; j++) {
	if (pmin[j] < min[j])
	  min[j] = pmin[j]; 
	if (pmax[j] > max[j])
	  max[j] = pmax[j]; 
      }
    }

  }

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"Overall volume bounds: min=[%.4e %.4e %.4e]\
   max=[%.4e %.4e %.4e]\n\n", min[0], min[1], min[2], max[0], max[1], max[2]);
#endif
 
  amr_list[0]->GetDims(block_dims); // assume all time steps are the same 
  level_mapping = new int*[num_timesteps]; 

  // reconcile metadata over all timesteps
  MatchTimestepLevels(); 

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// LoadData
//
// reads vector data independently for a time-varying flash AMR dataset
//
// fname: file containing list of timestep files
// start_block, end_block: contiguous range of blocks to read
// vx, vy, vz: names of three velocity components in the dataset
//
void TimeVaryingFlashAMR::LoadData(char* fname, int start_block, int end_block,
				  char *vx, char *vy, char *vz) {

  FILE *fIn; 
  char filename[300]; 
  int fo; // offset of start of current time step in the total block range
  int fs; // starting block number in current time step
  int fe; // ending block number in current time step
  int nb; // number of blocks in current time step
  int i;

  assert((fIn = fopen(fname, "r")) != NULL);
  fscanf(fIn, "%d", &num_timesteps); 

  fo = 0;

  for (i = 0; i < num_timesteps && start_block <= end_block; i++) {

    // start end end block within the current time step
    fs = start_block - fo;
    fe = end_block - fo;
    nb = amr_list[i]->GetNumBlocks();

    // block range does not include this file
    if (fs < 0 || fs >= nb) {
      fo += nb;
      continue;
    }

    // block range extends beyond end of this time step
    if (fe >= nb)
      fe = nb - 1;

    // read the blocks in the current time step
    fscanf(fIn, "%s", filename); 

#ifdef DEBUG
    if (myproc == 0)
      fprintf(stderr,"Reading vector data from %s ...\n\n", filename); 
#endif

    amr_list[i]->LoadHDF5Data(fs, fe, vx, vy, vz);

    // update for next time step
    start_block += (fe - fs + 1);
    fo += nb;

  }

}
//-----------------------------------------------------------------------
//
// MatchTimestepLevels
//
// Match levels from all timesteps (levels can change over time)
//
void TimeVaryingFlashAMR::MatchTimestepLevels() {

  int cnt = 0; // sum of levels in timesteps
  float *sizes; // x block size in each level
  float lsize[3]; // block size in level
  int *init; // initialization flag for each level
  float minB[3], maxB[3]; // extents of a block
  int i, j, k, l, n;

  for(i = 0; i < num_timesteps; i++) {
    n = amr_list[i]->GetMaxLevelValue();
    l = amr_list[i]->GetNumLevels();
    level_mapping[i] = new int[n + 1];
    cnt += l; 
  }

  sizes = new float[cnt]; 

  n = 0; 
  for (i = 0; i < num_timesteps; i++)  {
    for (j = 0; j < amr_list[i]->GetNumLevels(); j++) {
      amr_list[i]->GetLevelBlockSize(j, lsize); 
      sizes[n++] = lsize[0]; 
    }
  }
  
  float_sort_list(sizes, cnt); // cnt will be updated
  num_levels = cnt; // total different levels in all timesteps

#ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"There are total %d different levels \n", cnt); 
#endif

  // (re)allocate level size and extents for all levels  
  block_xsize_inLevel = new float[num_levels]; 
  block_ysize_inLevel = new float[num_levels]; 
  block_zsize_inLevel = new float[num_levels]; 
  level_minB = new float[num_levels * 3]; 
  level_maxB = new float[num_levels * 3]; 
  init = new int[num_levels]; // initialization flag for each level

  // match the levels 
  for (i = 0; i < num_timesteps; i++) {

    for (j = 0; j < amr_list[i]->GetNumLevels(); j++) {
 
	amr_list[i]->GetLevelBlockSize(j, lsize); 

      for (k = 0; k < num_levels; k++) {

	if (lsize[0] == sizes[k])  {
	  l = num_levels - k + 1;
	  level_mapping[i][l] = j; 
	  block_xsize_inLevel[level_mapping[i][l]] = lsize[0]; 
	  block_ysize_inLevel[level_mapping[i][l]] = lsize[1]; 
	  block_zsize_inLevel[level_mapping[i][l]] = lsize[2]; 
	  break; 
	}

      }
    }
  }

  // fix the levels to make them consistent 
  for (i = 0; i < num_timesteps; i++) {
    for (j = 0; j < amr_list[i]->GetNumBlocks(); j++)
      amr_list[i]->SetLevel(j, level_mapping[i][amr_list[i]->GetLevel(j)]); 
  }

  // set initial flag for each level
  for (i = 0; i < num_levels; i++)
    init[i] = 1;

  // extents for each level
  for (i = 0; i < num_timesteps; i++) {

    for (j = 0; j < amr_list[i]->GetNumLevels(); j++) {

      amr_list[i]->GetLevelBounds(j, minB, maxB); 

      if (init[j]) {
	level_minB[j * 3]     = minB[0]; 
	level_minB[j * 3 + 1] = minB[1]; 
	level_minB[j * 3 + 2] = minB[2]; 
	level_maxB[j * 3]     = maxB[0]; 
	level_maxB[j * 3 + 1] = maxB[1]; 
	level_maxB[j * 3 + 2] = maxB[2]; 
	init[j] = 0;
      }

      else {
	if (minB[0] < level_minB[j * 3])
	  level_minB[j * 3] = minB[0]; 
	if (minB[1] < level_minB[j * 3 + 1])
	  level_minB[j * 3 + 1] = minB[1]; 
	if (minB[2] < level_minB[j * 3 + 2])
	  level_minB[j * 3 + 2] = minB[2]; 
	if (maxB[0] > level_maxB[j * 3])
	  level_maxB[j * 3] = maxB[0]; 
	if (maxB[1] > level_maxB[j * 3 + 1])
	  level_maxB[j * 3 + 1] = maxB[1]; 
	if (maxB[2] > level_maxB[j * 3 + 2])
	  level_maxB[j * 3 + 2] = maxB[2]; 
      }

    }
  }

}
//-----------------------------------------------------------------------
//
// GetLevelBlockSize
//
// physical size of one block in a given level
//
void TimeVaryingFlashAMR::GetLevelBlockSize(int level, float size[3]) {

  size[0] = block_xsize_inLevel[level]; 
  size[1] = block_ysize_inLevel[level]; 
  size[2] = block_zsize_inLevel[level]; 

}
//-----------------------------------------------------------------------
//
// GetLevelBounds
//
// physical bounds of all blocks of a given level
//
void TimeVaryingFlashAMR::GetLevelBounds(int level, float minB[3], float maxB[3]) {

  minB[0] = level_minB[level*3]; 
  minB[1] = level_minB[level*3+1]; 
  minB[2] = level_minB[level*3+2]; 

  maxB[0] = level_maxB[level*3]; 
  maxB[1] = level_maxB[level*3+1]; 
  maxB[2] = level_maxB[level*3+2]; 

}
//-----------------------------------------------------------------------
