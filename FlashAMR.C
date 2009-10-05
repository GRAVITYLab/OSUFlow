

#include "FlashAMR.h" 

void sort_list(float* list, int& cnt) 
{
  float* tmp_list = new float[cnt]; 
  for (int i=0; i<cnt; i++) tmp_list[i] = list[i]; 
  for (int i=0; i<cnt-1; i++) {
    float largest = -1; 
    int idx = i; 
    for (int j = i; j<cnt; j++) {
      if (tmp_list[j]>largest) {idx = j; largest = tmp_list[j]; }
    }
    int tmp = tmp_list[i]; 
    tmp_list[i] = tmp_list[idx]; 
    tmp_list[idx] =tmp; 
  }
  list[0] = tmp_list[cnt-1]; 
  int count = 1; 
  for (int i=cnt-2; i>=0; i--) {
    if (tmp_list[i]>list[count-1]) {
      list[count] = tmp_list[i]; 
      count++; 
    }
  }
  cnt = count; 
}

FlashAMR::FlashAMR() 
{

}

FlashAMR::~FlashAMR() 
{

}

int FlashAMR::LoadData(char* fname, float min[3], float max[3]) 
{

  float *sizes; 
  FILE * in = fopen(fname, "r"); 

  if (in == NULL) return(0); 

  fread(&nb, sizeof(int), 1, in); // number of blocks
  printf(" %d blocks.\n", nb); 

  fread(block_dims, sizeof(int), 3, in); // dims of block, should be the same for all levels
  printf(" block dims %d %d %d\n", block_dims[0], block_dims[1], block_dims[2]); //16x16x16 for example 

  int size = block_dims[0]*block_dims[1]*block_dims[2];  // number of voxels in each block 
  
  block_level = new int[nb];        // the level that each block belongs to 
  block_index = new int[nb];        // the index of each block (defined by application)
  block_center = new float[nb*3];   // the center x y and z of each block 
  block_length = new float[nb*3];   // the physical length (in x y and z) of each block 
  block_minB = new float[nb*3];     // the min corner (min_x, min_y, min_z) of each block 
  block_maxB = new float[nb*3];     // the max corner (max_x, max_y, max_z) of each block 
  voxel_xsize = new float[nb]; // the physical length in x of a voxel in each block  
  voxel_ysize = new float[nb]; // the physical length in y of a voxel in each block 
  voxel_zsize = new float[nb]; // the physical length in z of a voxel in each block 

  vectors = new float*[nb]; 

  sizes = new float[nb];         // use one of the x y z legnths to determin block level 

  for (int i=0; i<nb; i++) {

    int idx = i*3; 

    fread(&(block_index[i]), sizeof(int), 1, in); 
    fread(&(block_center[idx]), sizeof(float), 3, in); 
    fread(&(block_length[idx]), sizeof(float), 3, in); 

    float bounds[6]; 
    fread(bounds, sizeof(float), 6, in);   // the two opposite corners of each block 
    block_minB[idx]   = bounds[0];  block_maxB[idx]   = bounds[1]; 
    block_minB[idx+1] = bounds[2];  block_maxB[idx+1] = bounds[3]; 
    block_minB[idx+2] = bounds[4];  block_maxB[idx+2] = bounds[5]; 

    if (i==0) {
      min[0] = block_minB[0]; min[1] = block_minB[1]; min[2] = block_minB[2]; 
      max[0] = block_maxB[0]; max[1] = block_maxB[1]; max[2] = block_maxB[2]; 
    }
    else {   // find the global min/max cornders 
      if (block_minB[idx]<min[0]) min[0] = block_minB[idx]; 
      if (block_minB[idx+1]<min[1]) min[1] = block_minB[idx+1];      
      if (block_minB[idx+2]<min[2]) min[2] = block_minB[idx+2]; 
      if (block_maxB[idx]>max[0]) max[0] = block_maxB[idx]; 
      if (block_maxB[idx+1]>max[1]) max[1] = block_maxB[idx+1]; 
      if (block_maxB[idx+2]>max[2]) max[2] = block_maxB[idx+2]; 
    }
    vectors[i] = new float[size*3]; 

    float *xcomp = new float[size]; 
    float *ycomp = new float[size]; 
    float *zcomp = new float[size]; 

    fread(xcomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(ycomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(zcomp, sizeof(float), size, in);  // read the data (vectors) 

    for (int c=0; c<size; c++) {
      vectors[i][c*3]   =xcomp[c]; 
      vectors[i][c*3+1] =ycomp[c]; 
      vectors[i][c*3+2] =zcomp[c]; 
    }
    delete [] xcomp; 
    delete [] ycomp; 
    delete [] zcomp; 

    // compute the length of each cell in the block 
    // will be used to determine the level which the block belongs to 
    sizes[i] = voxel_xsize[i] = block_length[idx]/(float)block_dims[0]; 
    voxel_ysize[i] = block_length[idx+1]/(float)block_dims[1]; 
    voxel_zsize[i] = block_length[idx+2]/(float)block_dims[2]; 
  }
  int count = nb; 
  sort_list(sizes, count); 
  num_levels = count; // the total levels in this data 
  level_minB = new float[num_levels*3];  // now determin the min/max corners of each level 
  level_maxB = new float[num_levels*3]; 

  block_xsize_inLevel = new float[num_levels]; //  the physical block_length in each level 
  block_ysize_inLevel = new float[num_levels]; 
  block_zsize_inLevel = new float[num_levels]; 

  for (int i=0; i<num_levels; i++) {
    int idx = i*3; 
    level_minB[idx] = level_minB[idx+1] = level_minB[idx+2] = 999999999; 
    level_maxB[idx] = level_maxB[idx+1] = level_maxB[idx+2] = -999999999; 
  }
  // write the level info for each block 
  // also find the block lengths in x y z for all levels, and 
  // the min max corners for all levels 

  for (int i=0; i<nb; i++) {
    for (int j=0; j<num_levels; j++) {
      if (voxel_xsize[i] == sizes[j]) {
	int level = num_levels-1-j; 
	block_level[i]= level; 
	block_xsize_inLevel[level] = block_length[i*3]; 
	block_ysize_inLevel[level] = block_length[i*3+1]; 
	block_zsize_inLevel[level] = block_length[i*3+2]; 

	int idx = i*3; 
	// update the level min corner 
	if (block_minB[idx]<level_minB[level*3]) 
	  level_minB[level*3] = block_minB[idx]; 
	if (block_minB[idx+1]<level_minB[level*3+1]) 
	  level_minB[level*3+1] = block_minB[idx+1]; 
	if (block_minB[idx+2]<level_minB[level*3+2]) 
	  level_minB[level*3+2] = block_minB[idx+2]; 

	// update the level max corner 
	if (block_maxB[idx]>level_maxB[level*3]) 
	  level_maxB[level*3] = block_maxB[idx]; 
	if (block_maxB[idx+1]>level_maxB[level*3+1]) 
	  level_maxB[level*3+1] = block_maxB[idx+1]; 
	if (block_maxB[idx+2]>level_maxB[level*3+2]) 
	  level_maxB[level*3+2] = block_maxB[idx+2]; 
      }
    }
  }
  return(1); 
}

void FlashAMR::GetLevelBlockSize(int level, float size[3]) 
{
  size[0] = block_xsize_inLevel[level]; 
  size[1] = block_ysize_inLevel[level]; 
  size[2] = block_zsize_inLevel[level]; 
}

void FlashAMR::GetLevelBounds(int level, float minB[3], float maxB[3])
{
  minB[0] = level_minB[level*3]; 
  minB[1] = level_minB[level*3+1]; 
  minB[2] = level_minB[level*3+2]; 

  maxB[0] = level_maxB[level*3]; 
  maxB[1] = level_maxB[level*3+1]; 
  maxB[2] = level_maxB[level*3+2]; 
}



/////////////////////////////////////////////////////////////// 
/////////////////////////////////////////////////////////////// 

TimeVaryingFlashAMR::TimeVaryingFlashAMR()
{

}

TimeVaryingFlashAMR::~TimeVaryingFlashAMR()
{
}

int TimeVaryingFlashAMR::LoadData(char* fname, float min[3], float max[3])
{
  FILE *fIn; 
  char filename[300]; 
  float pmin[3], pmax[3]; 
  
  fIn = fopen(fname, "r"); 
  assert(fIn !=NULL); 
  fscanf(fIn, "%d", &num_timesteps); 

  amr_list = new FlashAMR*[num_timesteps]; 

  for (int i = 0; i<num_timesteps; i++) {
    fscanf(fIn, "%s", filename); 
    printf(" to read %s ....\n", filename); 
    amr_list[i] = new FlashAMR; 
    amr_list[i]->LoadData(filename, pmin, pmax); 

    if (i==0) {
      min[0] = pmin[0]; min[1] = pmin[1]; min[2] = pmin[2]; 
      max[0] = pmax[0]; max[1] = pmax[1]; max[2] = pmax[2]; 
    }
    else {
      if (pmin[0] < min[0]) min[0] = pmin[0]; 
      if (pmin[1] < min[1]) min[1] = pmin[1]; 
      if (pmin[2] < min[2]) min[2] = pmin[2]; 

      if (pmax[0] > max[0]) max[0] = pmax[0]; 
      if (pmax[1] > max[1]) max[1] = pmax[1]; 
      if (pmax[2] > max[2]) max[2] = pmax[2]; 
    }
  }

  amr_list[0]->GetDims(block_dims); // assume all time steps are the same 
  level_mapping = new int*[num_timesteps]; 
  match_all_timesteps(); 
}


/////////////////////////////////////////////////////////////////////
//
//  since the grid is refined over time,  we need to match levels from all time
//  time steps 

void TimeVaryingFlashAMR::match_all_timesteps() {

  int cnt = 0; 

  for(int i=0; i<num_timesteps; i++) {
    int nlevels = amr_list[i]->GetNumLevels(); 
    level_mapping[i] = new int[nlevels]; 
    cnt += nlevels; 
  }
  float* sizes = new float[cnt]; 

  int idx = 0; 
  for (int i=0; i<num_timesteps; i++)  {
    int nlevels = amr_list[i]->GetNumLevels();     
    for (int j=0; j<nlevels; j++) {
      float lsizes[3]; 
      amr_list[i]->GetLevelBlockSize(j, lsizes); 
      sizes[idx++] = lsizes[0]; 
    }
  }
  
  sort_list(sizes, cnt); // cnt will be changed in the function. 
  printf(" there are total %d different levels \n", cnt); 
  
  num_levels = cnt; 
  
  block_xsize_inLevel = new float[num_levels]; 
  block_ysize_inLevel = new float[num_levels]; 
  block_zsize_inLevel = new float[num_levels]; 
  level_minB = new float[num_levels*3]; 
  level_maxB = new float[num_levels*3]; 

  // now match the levels 

  for (int i=0; i<num_timesteps; i++) {

    int nl = amr_list[i]->GetNumLevels(); 

    for (int j=0; j<nl; j++) 
      for (int l = 0; l<num_levels; l++) {
	float bsize[3]; 
	amr_list[i]->GetLevelBlockSize(j, bsize); 
	if (bsize[0] == sizes[l])  {
	  level_mapping[i][j] = num_levels-1-l; 
	  block_xsize_inLevel[num_levels-1-l] = bsize[0]; 
	  block_ysize_inLevel[num_levels-1-l] = bsize[1]; 
	  block_zsize_inLevel[num_levels-1-l] = bsize[2]; 
	  break; 
	}
      }
  }

  // now we need to fix the levels to make them consistent 

  for (int i=0; i<num_timesteps; i++) {
    int nblocks = amr_list[i]->GetNumBlocks(); 
    for (int j=0; j<nblocks; j++) {
      int lvl = amr_list[i]->GetLevel(j); 
      amr_list[i]->SetLevel(j, level_mapping[i][lvl]); 
    }
  }

  for (int i=0; i<num_levels; i++) {
    level_minB[i*3] = level_minB[i*3+1] = level_minB[i*3+2] = 999999999; 
    level_maxB[i*3] = level_maxB[i*3+1] = level_maxB[i*3+2] = -999999999; 
  }

  for (int i=0; i <num_timesteps; i++) {
    int nl = amr_list[i]->GetNumLevels(); 
    for (int j=0; j<nl; j++) {
      int level = level_mapping[i][j]; 
      float minB[3], maxB[3]; 
      amr_list[i]->GetLevelBounds(j, minB, maxB); 

      if (minB[0] < level_minB[level*3])
	level_minB[level*3] = minB[0]; 
      if (minB[1] < level_minB[level*3+1])
	level_minB[level*3+1] = minB[1]; 
      if (minB[2] < level_minB[level*3+2])
	level_minB[level*3+2] = minB[2]; 

      if (maxB[0] > level_maxB[level*3])
	level_maxB[level*3] = maxB[0]; 
      if (maxB[1] > level_maxB[level*3+1])
	level_maxB[level*3+1] = maxB[1]; 
      if (maxB[2] > level_maxB[level*3+2])
	level_maxB[level*3+2] = maxB[2]; 
    }
  }
  printf(" done.\n"); 
}


void TimeVaryingFlashAMR::GetLevelBlockSize(int level, float size[3]) 
{
  size[0] = block_xsize_inLevel[level]; 
  size[1] = block_ysize_inLevel[level]; 
  size[2] = block_zsize_inLevel[level]; 
}

void TimeVaryingFlashAMR::GetLevelBounds(int level, float minB[3], float maxB[3])
{
  minB[0] = level_minB[level*3]; 
  minB[1] = level_minB[level*3+1]; 
  minB[2] = level_minB[level*3+2]; 

  maxB[0] = level_maxB[level*3]; 
  maxB[1] = level_maxB[level*3+1]; 
  maxB[2] = level_maxB[level*3+2]; 
}
