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
// sorts a list of integers containing duplicate entries into distinct elements
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
// sorts a list of floats containing duplicate entries into distinct elements
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
// parallel reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
//
// fname: file name
// min, max: (output) global min, max corners
// block_size: (output) voxels in a block eg. 16x16x16
// comm: MPI communicator
//
void FlashAMR::ParallelLoadHDF5MetaData(char* fname, float min[3], float max[3],
					int *block_size, DataMode dm, 
					float scale, MPI_Comm comm) {

  fdf = new FlashHDFFile(fname, scale, comm); // flash file object
  LoadHDF5MetaData(min, max, block_size, dm);

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// serial reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
// currently exists in a serial version only
//
// fname: file name
// min, max: (output) global min, max corners
// block_size: (output) voxels in a block eg. 16x16x16
//
void FlashAMR::SerialLoadHDF5MetaData(char* fname, float min[3], float max[3],
				      int *block_size, DataMode dm, 
				      float scale) {

  fdf = new FlashHDFFile(fname, scale); // flash file object
  LoadHDF5MetaData(min, max, block_size, dm);

}
//-----------------------------------------------------------------------
//
// reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
// currently exists in a serial version only
//
// min, max: (output) global min, max corners
// block_size: (output) voxels in a block eg. 16x16x16
//
void FlashAMR::LoadHDF5MetaData(float min[3], float max[3], int *block_size,
				DataMode dm) {

  float *bounds; // bounding box: min, max corners
  int *init; // initial flag for each level
  int level; // current level
  int max_level; // maximum level
  int i, j;

  // we only support HDF5 AMR data
  assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);

  // number of blocks and block dimensions
  nb = fdf->GetNumberOfBlocks(); // includes non-leaf blocks
  file_nb = nb; // save the total number of file block (leaf + nonleaf)
  fdf->GetCellDimensions(block_dims);
  block_size[0] = block_dims[0];
  block_size[1] = block_dims[1];
  block_size[2] = block_dims[2];

  // allocate memory for leaf and non-leaf blocks
  block_level = new int[nb]; // level that each block belongs to 
  orig_level = new int[nb]; // local copy of block_level
  block_center = new float[nb * 3]; // center x y and z of each block 
  block_length = new float[nb * 3]; // physical length (x, y, z) of each block 
  block_minB = new float[nb * 3]; // min corner (x, y, z) of each blk
  block_maxB = new float[nb * 3]; // max corner (x, y, z) of each blk
  bounds = new float[nb * 6]; // bounding box of each blk

  // center, length, level, extents of each block
  if (dm == HDF_FLOAT) {
    nb = fdf->GetFloatVecBlocks((char *)"coordinates", 0, file_nb, 3, 1, 
				block_center);
    fdf->GetFloatVecBlocks((char *)"block size", 0, file_nb, 3, 1, block_length);
    fdf->GetIntVecBlocks((char *)"refine level", 0, file_nb, 1, 1, block_level);
    fdf->GetFloatMatBlocks((char *)"bounding box", 0, file_nb, 3, 2, 1, bounds);
  } else {
    nb = fdf->GetDoubleVecBlocks((char *)"coordinates", 0, file_nb, 3, 1, 
				block_center);
    fdf->GetDoubleVecBlocks((char *)"block size", 0, file_nb, 3, 1, block_length);
    fdf->GetIntVecBlocks((char *)"refine level", 0, file_nb, 1, 1, block_level);
    fdf->GetDoubleMatBlocks((char *)"bounding box", 0, file_nb, 3, 2, 1, bounds);
  }

  for (i = 0; i < nb; i++) {
    orig_level[i] = block_level[i];
    for (j = 0; j < 3; j++) {
      block_minB[3 * i + j] = bounds[2 * (3 * i + j)];
      block_maxB[3 * i + j] = bounds[2 * (3 * i + j) + 1];
    }
//     fprintf(stderr,"min = [%.2e %.2e %.2e] max = [%.2e %.2e %.2e]\n",block_minB[3*i],block_minB[3*i+1],block_minB[3*i+2],block_maxB[3*i],block_maxB[3*i+1],block_maxB[3*i+2]);
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
  int_sort_list(orig_level, num_levels);
  max_level = orig_level[num_levels - 1];

  // extents and block size for each level
  // block size is physical size (length)
  init = new int[max_level + 1];
  for (i = 0; i < max_level + 1; i++)
    init[i] = 1;
  level_minB = new float[(max_level + 1) * 3];
  level_maxB = new float[(max_level + 1) * 3]; 
  block_xsize_inLevel = new float[max_level + 1];
  block_ysize_inLevel = new float[max_level + 1]; 
  block_zsize_inLevel = new float[max_level + 1]; 
  for (i = 0; i < nb; i++) {

    level = block_level[i];

    block_xsize_inLevel[level] = block_length[3 * i]; 
    block_ysize_inLevel[level] = block_length[3 * i + 1]; 
    block_zsize_inLevel[level] = block_length[3 * i + 2]; 

    if (init[level]) {
      for (j = 0; j < 3; j++) {
	level_minB[3 * level + j] = block_minB[3 * i + j];
	level_maxB[3 * level + j] = block_maxB[3 * i + j];
      }
      init[level] = 0;
    }
    else{
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
// reader for HDF5 Flash vector data from 
// reads start_block to end_block (inclusive)
//
// start_block, end_block: contiguous range of blocks to read
// vx, vy, vz: names of three velocity components in the dataset
// mpi: whether mpi is enabled
//
// returns: I/O bandwidth in MB/s
//
int FlashAMR::LoadHDF5Data(int start_block, int end_block,
			   char * vx, char * vy, char *vz, 
			   DataMode dm, int mpi) {

  int size = block_dims[0] * block_dims[1] * block_dims[2]; // total voxels
  float *comp; // one component of vector data
  int i, j;
  double t = 0.0; // timer

  // we only support HDF5 AMR data
  assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);

  int nblocks = end_block - start_block + 1;

  // allocate memory for leaf + nonleaf blocks
  vectors = new float*[nb]; // the vector data
  comp = new float[file_nb * size]; // don't know how much data will need to be
                                    //  read in order to get nblocks leaf blocks
                                    // todo: find a tighter bound for comp?

  for (i = 0; i < nb; i++)
    vectors[i] = new float[size * 3];

#ifdef _MPI
  if (mpi)
    t = MPI_Wtime();
#endif

  // x component
  if (dm == HDF_FLOAT)
    fdf->GetFloatVolBlocks(vx, start_block, nblocks, block_dims[0], 
			   block_dims[1], block_dims[2], 1, comp);
  else
    fdf->GetDoubleVolBlocks(vx, start_block, nblocks, block_dims[0], 
			    block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {

    vectors[start_block + i] = new float[size * 3]; // large enough for all
                                                    // three components
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 0] = comp[i * size + j]; 
  }

  // y component
  if (dm == HDF_FLOAT)
    fdf->GetFloatVolBlocks(vy, start_block, nblocks, block_dims[0], 
			   block_dims[1], block_dims[2], 1, comp);
  else
    fdf->GetDoubleVolBlocks(vy, start_block, nblocks, block_dims[0], 
			    block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 1] = comp[i * size + j]; 
  }

  // z component
  if (dm == HDF_FLOAT)
    fdf->GetFloatVolBlocks(vz, start_block, nblocks, block_dims[0], 
			   block_dims[1], block_dims[2], 1, comp);
  else
    fdf->GetDoubleVolBlocks(vz, start_block, nblocks, block_dims[0], 
			    block_dims[1], block_dims[2], 1, comp);
  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < size; j++)
      vectors[start_block + i][j * 3 + 2] = comp[i * size + j]; 
  }

  delete [] comp; 
  fdf->Close();

  // return I/O bandwidth
#ifdef _MPI
  if(mpi)
    t = MPI_Wtime() - t;
  return(nblocks * size * 12 / t / 1048576); // 3 components * 4 bytes = 12
#else
  return(0); // no bandwidth to report without MPI_Wtime
#endif

}
//-----------------------------------------------------------------------
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
// reads metadata collectively for a time-varying flash AMR dataset
//
// files: timestep file names
// min, max: (output) global min, max extent over entire dataset (space+time)
// tsize: number of timesteps
// block_size: (output) voxels in a block eg. 16x16x16
// dm: data mode
//
// side effects: allocates an amr object for each timestep
//
void TimeVaryingFlashAMR::LoadMetaData(char** files, float *min, float *max,
				       int tsize, int *block_size,
				       DataMode dm, float scale, 
				       MPI_Comm comm) {

  float pmin[3], pmax[3]; 
  int i, j;

  // we only support HDF5 AMR data
  assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);

  // get the number of timesteps, allocate amr objects
  num_timesteps = tsize;
  amr_list = new FlashAMR*[num_timesteps]; 

  // load the metadata (extents, level info)
  for (i = 0; i < num_timesteps; i++) {

    amr_list[i] = new FlashAMR(myproc);
    amr_list[i]->ParallelLoadHDF5MetaData(files[i], pmin, pmax, block_size,
					  dm, scale, comm); 

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

// #ifdef DEBUG
  if (myproc == 0)
    fprintf(stderr,"Overall volume bounds: min=[%.4e %.4e %.4e]\
   max=[%.4e %.4e %.4e]\n\n", min[0], min[1], min[2], max[0], max[1], max[2]);
// #endif
 
  amr_list[0]->GetDims(block_dims); // assume all time steps are the same 

  // reconcile metadata over all timesteps
  MatchTimestepLevels(); 

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// reads metadata for a time-varying flash AMR dataset
//
// files: timestep file names
// min, max: (output) global min, max extent over entire dataset (space+time)
// tsize: number of timesteps
// block_size: (output) voxels in a block eg. 16x16x16
// dm: data_mode
//
// side effects: allocates an amr object for each timestep
//
void TimeVaryingFlashAMR::SerialLoadMetaData(char** files, float min[3], 
					     float max[3], int tsize,
					     int *block_size, DataMode dm, 
					     float scale) {

  float pmin[3], pmax[3]; 
  int unused;
  int i, j;

  // we only support HDF5 AMR data
  assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);

  // get the number of timesteps, allocate amr objects
  num_timesteps = tsize;
  amr_list = new FlashAMR*[num_timesteps]; 

  // load the metadata (extents, level info)
  for (i = 0; i < num_timesteps; i++) {

    amr_list[i] = new FlashAMR(myproc);
    amr_list[i]->SerialLoadHDF5MetaData(files[i], pmin, pmax, block_size, 
					dm, scale); 

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

  // reconcile metadata over all timesteps
  MatchTimestepLevels(); 

}
//-----------------------------------------------------------------------
//
// DEPRECATED
//
// //
// // reads vector data for a time-varying flash AMR dataset
// //
// // start_block, end_block: contiguous range of blocks to read
// // vx, vy, vz: names of three velocity components in the dataset
// // mpi: whether mpi is enabled
// //
// // returns: average I/O bandwidth in MB/s across all time steps
// //
// int TimeVaryingFlashAMR::LoadData(int start_block, int end_block,
// 				  char *vx, char *vy, char *vz, 
// 				  DataMode dm, int mpi) {

//   int fo; // offset of start of current time step in the total block range
//   int fs; // starting block number in current time step
//   int fe; // ending block number in current time step
//   int nb; // number of blocks in current time step
//   int i;
//   int io_bw = 0; // I/O bandwidth in MB/s
//   int n = 0; // number of reads performed

//   // we only support HDF5 AMR data
//   assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);

//   fo = 0;

//   // num_timesteps was assigned earlier during LoadMetaData
//   for (i = 0; i < num_timesteps && start_block <= end_block; i++) {

//     // start end end block within the current time step
//     fs = start_block - fo;
//     fe = end_block - fo;
//     nb = amr_list[i]->GetNumBlocks();

//     // block range does not include this file
//     if (fs < 0 || fs >= nb) {
//       fo += nb;
//       continue;
//     }

//     // block range extends beyond end of this time step
//     if (fe >= nb)
//       fe = nb - 1;

//     io_bw += amr_list[i]->LoadHDF5Data(fs, fe, vx, vy, vz, dm, mpi);
//     n++;

//     // update for next time step
//     start_block += (fe - fs + 1);
//     fo += nb;

//   }

//   // return average I/O bandwidth across all time steps
//   return(io_bw / n);

// }
//-----------------------------------------------------------------------
//
// Match levels from all timesteps (levels can change over time)
//
void TimeVaryingFlashAMR::MatchTimestepLevels() {

  int cnt; // sum of levels in timesteps
  float lsize[3]; // block size in level
  float minB[3], maxB[3]; // extents of a block
  int i, j, k, n;
  int *orig_level; // original levels in all time steps
  int max_level; // maximum original level

  // space for total number of levels across all timesteps, including duplicates
  cnt = 0;
  for(i = 0; i < num_timesteps; i++)
    cnt += amr_list[i]->GetNumLevels();
  orig_level = new int[cnt]; 

  n = 0; 
  for (i = 0; i < num_timesteps; i++)  {
    for (j = 0; j < amr_list[i]->GetNumLevels(); j++)
      orig_level[n++] = amr_list[i]->orig_level[j];
  }
  
  int_sort_list(orig_level, cnt); // cnt will be updated
  num_levels = cnt; // total different levels in all timesteps
  max_level = orig_level[num_levels - 1];
  lmap = new int[max_level + 1];

  // renumber levels canonically, 0 (coarse) to num_levels - 1 (fine)
  for (i = 0; i < num_levels; i++)
    lmap[orig_level[i]] = i;
  for (i = 0; i < num_timesteps; i++) {
    for (j = 0; j < amr_list[i]->GetNumBlocks(); j++)
      amr_list[i]->SetLevel(j, lmap[amr_list[i]->GetLevel(j)]); 
  }

  // level sizes and extents
  block_xsize_inLevel = new float[num_levels]; 
  block_ysize_inLevel = new float[num_levels]; 
  block_zsize_inLevel = new float[num_levels]; 
  level_minB = new float[num_levels * 3]; 
  level_maxB = new float[num_levels * 3]; 
  for (i = 0; i < num_timesteps; i++) {
    for (j = 0; j < amr_list[i]->GetNumLevels(); j++) {
      amr_list[i]->GetLevelBlockSize(amr_list[i]->orig_level[j], lsize); 
      k = lmap[amr_list[i]->orig_level[j]];
      block_xsize_inLevel[k] = lsize[0]; 
      block_ysize_inLevel[k] = lsize[1]; 
      block_zsize_inLevel[k] = lsize[2]; 
      amr_list[i]->GetLevelBounds(amr_list[i]->orig_level[j], minB, maxB); 
      level_minB[k * 3]     = minB[0]; 
      level_minB[k * 3 + 1] = minB[1]; 
      level_minB[k * 3 + 2] = minB[2]; 
      level_maxB[k * 3]     = maxB[0]; 
      level_maxB[k * 3 + 1] = maxB[1]; 
      level_maxB[k * 3 + 2] = maxB[2]; 
    }
  }

}
//-----------------------------------------------------------------------
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
