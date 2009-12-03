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

#ifndef _FLASHAMR_H_ 
#define _FLASHAMR_H_ 

#include <stdio.h> 
#include <stdlib.h> 
#include <assert.h>
#include "flashhdf5_float.h"

#ifdef _MPI
#include <mpi.h>
#endif

//--------------------------------------------------------------------------
//
// utility functions
//
void int_sort_list(int* list, int& cnt);
void float_sort_list(float* list, int& cnt);

//--------------------------------------------------------------------------
//
// Static FlashAMR class
//
class FlashAMR {

 public: 

  FlashAMR(int myid = 0){ myproc = myid; }
  ~FlashAMR(){} 
  int GetNumBlocks() {return nb;}
  int GetLevel(int i) {return block_level[i];}
  void SetLevel(int i, int level) {block_level[i] = level;}
  float* GetBlockCenter(int i) {return &(block_center[i*3]); }
  float* GetBlockLengths(int i) {return &(block_length[i*3]); }
  int GetNumLevels() {return num_levels; }
  float *GetDataPtr(int i) {return vectors[i]; }
  void GetLevelBlockSize(int level, float size[3]); 
  void GetLevelBounds(int level, float minB[3], float maxB[3]); 
  void GetDims(int dims[3]) {dims[0] = block_dims[0]; dims[1] = block_dims[1];                               dims[2] = block_dims[2]; }
  int GetMaxLevelValue() { return max_level_value; };

#ifdef _MPI
  void ParallelLoadHDF5MetaData(char* fname, float* min, float* max, 
			       MPI_Comm comm); 
#endif
  void SerialLoadHDF5MetaData(char* fname, float* min, float* max); 

  void LoadHDF5MetaData(float* min, float* max); 
  int LoadHDF5Data(int start_block, int end_block, 
		   char *vx, char *vy, char *vz); 

 private: 

  FlashHDFFile *fdf; // flash data file
  int nb; // number of leaf blocks that are kept
  int file_nb; // total number of leaf + nonleaf blocks in the file
  int num_levels; 
  int max_level_value; // maximum level number in file

  int block_dims[3]; 
  int *block_index; 
  int *block_level; 
  float *block_center; 
  float *block_length; 
  float *block_minB, *block_maxB; 
  float *level_minB; 
  float *level_maxB; 

  float *block_xsize_inLevel; 
  float *block_ysize_inLevel; 
  float *block_zsize_inLevel; 

  float **vectors; 

  int myproc; // my process or thread number (0 if serial)

} ; 
//--------------------------------------------------------------------------
//
// Static FlashAMR class
//
class TimeVaryingFlashAMR {

 public: 

  TimeVaryingFlashAMR(int myid = 0){ myproc = myid; }
  ~TimeVaryingFlashAMR(){} 
  
#ifdef _MPI
  void LoadMetaData(char *fname, float*, float*, MPI_Comm = MPI_COMM_WORLD); 
#else
  void LoadMetaData(char *fname, float*, float*); 
#endif
  int LoadData(char *fname, int start_block, int end_block, char *vx, char *vy,
	       char *vz); 

  int GetNumTimeSteps() {return num_timesteps; }
  int GetNumLevels() {return num_levels; } 
  FlashAMR* GetTimeStep(int t) {if (t>=0 && t<num_timesteps) return amr_list[t]; 
                                else return (NULL); } 

  void GetLevelBlockSize(int level, float size[3]); 
  void GetLevelBounds(int level, float minB[3], float maxB[3]); 
  void GetDims(int dims[3]) {dims[0] = block_dims[0]; dims[1] = block_dims[1];                               dims[2] = block_dims[2]; }


 private: 

  int block_dims[3]; 

  int num_timesteps; 
  FlashAMR** amr_list; 
  int num_levels; 

  float *level_minB, *level_maxB; 
  int **level_mapping; 
  
  float *block_xsize_inLevel; 
  float *block_ysize_inLevel; 
  float *block_zsize_inLevel; 
  
  void MatchTimestepLevels(); 

  int myproc; // my process or thread number (0 if serial)

}; 
//--------------------------------------------------------------------------

#endif 
