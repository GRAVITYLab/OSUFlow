
#include "flashhdf5_float.h"
#include "FlashAMR.h" 

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
//
// LoadHDF5MetaData
//
// reader for HDF5 Flash metadata
// computes level info and extents for a single timestep
// currently exists in a serial version only
//
// fname: file name
// min, max: (output) global min, max corners
//
int FlashAMR::LoadHDF5MetaData(char* fname, float min[3], float max[3]) {

  FlashHDFFile fdf(fname); // flash file object
  float bounds[6]; // bounding box: min, max corners
  int n = 0; // total number of (leaf) blocks read and kept
  int *init; // initial flag for each level
  int *lmap; // map of levels in file to position in block_level list
  int *local_block_level; // local version of block_level
  int level; // current level
  int i, j;

  // number of blocks and block dimensions
  nb = fdf.GetNumberOfBlocks(); // includes non-leaf blocks
  file_nb = nb; // save the total number of file block (leaf + nonleaf)
  fdf.GetCellDimensions(block_dims);
  if (myproc == 0) {
    fprintf(stderr,"Total number of leaf + nonleaf blocks = %d\n", nb); 
    fprintf(stderr,"Each block is %dx%dx%d cells\n",
	    block_dims[0], block_dims[1], block_dims[2]);
  }

  // allocate memory  
  // allocates more than necessary for now, some blocks are non-leaf
  // and not stored; no memory needed for them
  block_level = new int[nb]; // level that each block belongs to 
  local_block_level = new int[nb]; // local copy of block_level
  block_center = new float[nb * 3]; // center x y and z of each block 
  block_length = new float[nb * 3]; // physical length (x, y, z) of each block 
  block_minB = new float[nb * 3]; // min corner (x, y, z) of each blk
  block_maxB = new float[nb * 3]; // max corner (x, y, z) of each blk

  // center, length, level, extents of each block
  // reading one block at a time is inefficient--
  // todo: change to reading the entire dataset at once
  for (i = 0; i < nb; i++) {

    // skip non-leaf blocks
    if (fdf.GetNodeType(i) != 1)
      continue;

    fdf.Get3dCoordinate(i, &(block_center[3 * n]));
    fdf.Get3dBlockSize(i, &(block_length[3 * n]));
    block_level[n] = fdf.GetRefinementLevel(i);
    local_block_level[n] = block_level[n];
    fdf.Get3dBoundingBox(i, bounds);

    for (j = 0; j < 3; j++) {
      block_minB[3 * n + j] = bounds[2 * j];
      block_maxB[3 * n + j] = bounds[2 * j + 1];
    }

    n++;

  }

  // list of distinct levels
  num_levels = n;
  int_sort_list(local_block_level, num_levels);
  max_level_value = local_block_level[num_levels - 1]; // needed later for
                                                // level map in this time step
                                                //  and in merged time steps

  // overall bounds of the entire data
  fdf.GetCoordinateRangeEntireDataset(bounds);
  for (i = 0; i < 3; i++) {
    min[i] = bounds[2 * i];
    max[i] = bounds[2 * i + 1];
  }

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

  n = 0;

  for (i = 0; i < nb; i++) {

    // skip non-leaf blocks
    if (fdf.GetNodeType(i) != 1)
      continue;

    // level = position in block_level list of 
    // refinement level in file
    level = lmap[block_level[n]];

    // level block size
    block_xsize_inLevel[level] = block_length[3 * n]; 
    block_ysize_inLevel[level] = block_length[3 * n + 1]; 
    block_zsize_inLevel[level] = block_length[3 * n + 2]; 

    // level extents
    if (init[level]) {
      for (j = 0; j < 3; j++) {
	level_minB[3 * level + j] = block_minB[3 * n + j];
	level_maxB[3 * level + j] = block_maxB[3 * n + j];
      }
      init[level] = 0;
    }
    else {
      for (j = 0; j < 3; j++) {
	if (block_minB[3 * n + j] < level_minB[3 * level + j]) 
	  level_minB[3 * level + j] = block_minB[3 * n + j]; 
	if (block_maxB[3 * n + j] > level_maxB[3 * level + j]) 
	  level_maxB[3 * level + j] = block_maxB[3 * n + j]; 
      }
    }

    n++;

  }

  nb = n; // store the actual number of leaf blocks that are kept

  delete [] init;

}
//-----------------------------------------------------------------------
//
// SerialLoadHDF5Data
//
// serial reader for HDF5 Flash vector data
//
// fname: file name
//
int FlashAMR::SerialLoadHDF5Data(char* fname) {

  int size = block_dims[0] * block_dims[1] * block_dims[2]; // total voxels
  int n = 0; // total number of (leaf) blocks read and kept
  float *xcomp, *ycomp, *zcomp; // x,y,z components of vector data
  FlashHDFFile fdf(fname);
  char scalar1[] = "velx"; // hard-coded variable names for now
  char scalar2[] = "vely";
  char scalar3[] = "velz";
  int i, j;

  // allocate memory  
  // allocates more than necessary for now, some blocks are non-leaf
  // and not stored; no memory needed for them
  vectors = new float*[nb]; // the vector data
  xcomp = new float[size]; 
  ycomp = new float[size]; 
  zcomp = new float[size]; 

  // read all blocks
  for (i = 0; i < file_nb; i++) {

    // skip non-leaf blocks
    if (fdf.GetNodeType(i) != 1)
      continue;

    vectors[n] = new float[size * 3]; 
    fdf.GetScalarVariable(scalar1, i, xcomp);
    fdf.GetScalarVariable(scalar2, i, ycomp);
    fdf.GetScalarVariable(scalar3, i, zcomp);
    for (j = 0; j < size; j++) {
      vectors[n][j * 3 + 0] = xcomp[j]; 
      vectors[n][j * 3 + 1] = ycomp[j]; 
      vectors[n][j * 3 + 2] = zcomp[j]; 
    }

    n++;

  }

  delete [] xcomp; 
  delete [] ycomp; 
  delete [] zcomp; 

  return(1); 

}
//-----------------------------------------------------------------------
//
// ParallelLoadHDF5Data
//
// parallel reader for HDF5 Flash vector data
//
// fname: file name
// start_block, end_block: contiguous range of blocks to read
//
int FlashAMR::ParallelLoadHDF5Data(char* fname, int start_block, 
				   int end_block) {

  int size = block_dims[0] * block_dims[1] * block_dims[2]; // total voxels
  int n = start_block; // index of (leaf) blocks read and kept
  float *xcomp, *ycomp, *zcomp; // x,y,z components of vector data
  FlashHDFFile fdf(fname);
  char scalar1[] = "velx"; // hard-coded variable names for now
  char scalar2[] = "vely";
  char scalar3[] = "velz";
  int start_file_block, end_file_block; // file blocks (leaf + nonleaf)
                                        // corresponding to start/end blocks
  int num_file_block; // number of blocks to read from the file (leaf + nonleaf)
  int i, j;

  // allocate memory  
  // allocates more than necessary for now, some blocks are non-leaf
  // and not stored; no memory needed for them
  vectors = new float*[nb]; // the vector data
  xcomp = new float[size]; 
  ycomp = new float[size]; 
  zcomp = new float[size]; 

  // find file blocks (leaf + nonleaf) corresponding to start/end leaf blocks
  // todo: store this in memory to avoid another file access
  j = 0;
  for (i = 0; i < file_nb; i++) {
    if (fdf.GetNodeType(i) != 1)
      continue;
    if (j == start_block)
      break;
    j++;
  }
  start_file_block = i;
  for ( ; i < file_nb; i++) {
    if (fdf.GetNodeType(i) != 1)
      continue;
    if (j == end_block)
      break;
    j++;
  }
  end_file_block = i;
  num_file_block = end_file_block - start_file_block + 1;

  // read my blocks
  for (i = start_file_block; i <= end_file_block; i++) {

    // skip non-leaf blocks
    if (fdf.GetNodeType(i) != 1)
      continue;

    vectors[n] = new float[size * 3]; 
    fdf.GetScalarVariable(scalar1, i, xcomp);
    fdf.GetScalarVariable(scalar2, i, ycomp);
    fdf.GetScalarVariable(scalar3, i, zcomp);
    for (j = 0; j < size; j++) {
      vectors[n][j * 3 + 0] = xcomp[j]; 
      vectors[n][j * 3 + 1] = ycomp[j]; 
      vectors[n][j * 3 + 2] = zcomp[j]; 
    }

    n++;

  }

  delete [] xcomp; 
  delete [] ycomp; 
  delete [] zcomp; 

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
//
// LoadMetaData
//
// reads metadata for a time-varying flash AMR dataset
//
// fname: file containing list of timestep files
// min, max: (output) global min, max extent over entire dataset (space+time)
//
int TimeVaryingFlashAMR::LoadMetaData(char* fname, float min[3], float max[3]) {

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
    if (myproc == 0)
      fprintf(stderr,"Reading metadata from %s ...\n\n", filename); 
    amr_list[i] = new FlashAMR(myproc);
    amr_list[i]->LoadHDF5MetaData(filename, pmin, pmax); 

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

  if (myproc == 0)
    fprintf(stderr,"Overall volume bounds: min=[%.4e %.4e %.4e]\
   max=[%.4e %.4e %.4e]\n\n", min[0], min[1], min[2], max[0], max[1], max[2]);
 
  amr_list[0]->GetDims(block_dims); // assume all time steps are the same 
  level_mapping = new int*[num_timesteps]; 

  // reconcile metadata over all timesteps
  MatchTimestepLevels(); 

}
//-----------------------------------------------------------------------
//
// LoadData
//
// reads vector data for a time-varying flash AMR dataset
//
// fname: file containing list of timestep files
// start_block, end_block: contiguous range of blocks to read
//
int TimeVaryingFlashAMR::LoadData(char* fname, int start_block,
				  int end_block) {

  FILE *fIn; 
  char filename[300]; 
  int i;

  assert((fIn = fopen(fname, "r")) != NULL);
  fscanf(fIn, "%d", &num_timesteps); 

  for (i = 0; i < num_timesteps; i++) {
    fscanf(fIn, "%s", filename); 
    if (myproc == 0)
      fprintf(stderr,"Reading vector data from %s ...\n\n", filename); 
    amr_list[i]->ParallelLoadHDF5Data(filename, start_block, end_block);
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
  if (myproc == 0)
    fprintf(stderr,"There are total %d different levels \n", cnt); 

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
