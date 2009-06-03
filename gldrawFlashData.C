
#include <stdio.h>
#include <stdlib.h> 

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#else
#include <GL/glut.h> 
#include <GL/gl.h>
#endif 

#include "OSUFlow.h"
#include "LatticeAMR.h" 

///////////////////////////////
//
// AMR book-keeping variables 
//
int num_timesteps;     // number of time steps 

float domain_lengths[3];      // x y z lengths of the entire domain 

float domain_center[3];       // center of the entire domain 


//***********************************************
//    Per time step information 
//

int   *num_levels;       // total number of levels in each time step 

int   *nb;               // nb[t], total number of blocks in each time step 

////////////////////////////////
//  Per block information 

int   block_dims[3];     // dimensions of a data block; assumes all levels all time steps are the same 
int   **block_index;     // block id 
int   **block_level;     // the level that each block belongs to 
float **block_center;    // center of each block 
float **block_length;    // physical block_length (xyz) of each block
float **block_minB;      // min/max corners for every block 
float **block_maxB;        

float **voxel_xsize_inBlock;     // x length of a voxel in a block for every block
float **voxel_ysize_inBlock;     // y length of a voxel in a block for every block
float **voxel_zsize_inBlock;     // z length of a voxel in a block for every block

/////////////////////////////////
//  per level information 
//
float **block_xsize_inLevel;     // the physical length of a block in each level at every time step
float **block_ysize_inLevel;     // all blocks in a level should have the same length 
float **block_zsize_inLevel; 

float ***level_minB;             // the min corner for each level at every time step 
float ***level_maxB;             // the max corner for each level at every time step 

float ***vectors;                 // data for blocks in all time steps

//***********************************************
// Aggregate information across all time steps

int    total_levels;      // how many levels in total across all time steps 
int    **level_mapping;   // a mapping from level number in each time step to the global level number 

float  *g_block_xsize_inLevel; // the physical length of a block in each level 
float  *g_block_ysize_inLevel; 
float  *g_block_zsize_inLevel; 

float  **g_level_minB;    // the coordinates of min corner of each level over all time steps 
float  **g_level_maxB;    // the coordinates of max corner of each level over all time steps 

//////////////////////////////////////////////////////
//  Lattice data structure built on top of the grid 

LatticeAMR * latticeAMR; 

//////////////////////////////////////////////////////
//  The subdomain parition 
//
volume_bounds_type_f* vb_list; 


//////////////////////////////////////////////////////
OSUFlow **osuflow_list; 

list<vtListSeedTrace*> *sl_list; 
VECTOR3 **osuflow_seeds; 
int *osuflow_num_seeds; 


int total_seeds = 5000; 
int npart; 
int first_frame = 1; 

////////////////////////////////////////////////////////


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

/////////////////////////////////////////////////////////////////////
//
//   A routine to read one time step of data from the Flash application generated 
//   by Mike Papka's group. This format is not the original Flash format
//   but just an in-house format for quick development 
// 
void LoadAMR(char* fname, float* min, float *max, int t)
{
  float *sizes; 
  FILE * in = fopen(fname, "r"); 

  fread(&nb[t], sizeof(int), 1, in); // number of blocks
  printf(" %d blocks.\n", nb[t]); 

  fread(block_dims, sizeof(int), 3, in); // dims of block, should be the same for all levels
  printf(" block dims %d %d %d\n", block_dims[0], block_dims[1], block_dims[2]); //16x16x16 for example 

  int size = block_dims[0]*block_dims[1]*block_dims[2];  // number of voxels in each block 
  
  block_level[t] = new int[nb[t]];        // the level that each block belongs to 
  block_index[t] = new int[nb[t]];        // the index of each block (defined by application)
  block_center[t] = new float[nb[t]*3];   // the center x y and z of each block 
  block_length[t] = new float[nb[t]*3];   // the physical length (in x y and z) of each block 
  block_minB[t] = new float[nb[t]*3];     // the min corner (min_x, min_y, min_z) of each block 
  block_maxB[t] = new float[nb[t]*3];     // the max corner (max_x, max_y, max_z) of each block 
  voxel_xsize_inBlock[t] = new float[nb[t]]; // the physical length in x of a voxel in each block  
  voxel_ysize_inBlock[t] = new float[nb[t]]; // the physical length in y of a voxel in each block 
  voxel_zsize_inBlock[t] = new float[nb[t]]; // the physical length in z of a voxel in each block 

  vectors[t] = new float*[nb[t]]; 

  sizes = new float[nb[t]];         // use one of the x y z legnths to determin block level 

  for (int i=0; i<nb[t]; i++) {

    int idx = i*3; 

    fread(&(block_index[t][i]), sizeof(int), 1, in); 
    //    printf(" read block id: %d\n", block_index[t][i]); 

    fread(&(block_center[t][idx]), sizeof(float), 3, in); 
    //    printf(" block %d center: [%f %f %f]\n", block_index[t][i], 
    //	   block_center[t][idx], block_center[t][idx+1], block_center[t][idx+2]); 

    fread(&(block_length[t][idx]), sizeof(float), 3, in); 
    //    printf(" block %d length: [%f %f %f]\n", block_index[t][i], 
    //	   block_length[t][idx], block_length[t][idx+1], block_length[t][idx+2]); 

    float bounds[6]; 
    fread(bounds, sizeof(float), 6, in);   // the two opposite corners of each block 
    block_minB[t][idx]   = bounds[0];  block_maxB[t][idx]   = bounds[1]; 
    block_minB[t][idx+1] = bounds[2];  block_maxB[t][idx+1] = bounds[3]; 
    block_minB[t][idx+2] = bounds[4];  block_maxB[t][idx+2] = bounds[5]; 

    if (i==0) {
      min[0] = block_minB[t][0]; min[1] = block_minB[t][1]; min[2] = block_minB[t][2]; 
      max[0] = block_maxB[t][0]; max[1] = block_maxB[t][1]; max[2] = block_maxB[t][2]; 
    }
    else {   // find the global min/max cornders 
      if (block_minB[t][idx]<min[0]) min[0] = block_minB[t][idx]; 
      if (block_minB[t][idx+1]<min[1]) min[1] = block_minB[t][idx+1];      
      if (block_minB[t][idx+2]<min[2]) min[2] = block_minB[t][idx+2]; 
      if (block_maxB[t][idx]>max[0]) max[0] = block_maxB[t][idx]; 
      if (block_maxB[t][idx+1]>max[1]) max[1] = block_maxB[t][idx+1]; 
      if (block_maxB[t][idx+2]>max[2]) max[2] = block_maxB[t][idx+2]; 
    }
    vectors[t][i] = new float[size*3]; 

    float *xcomp = new float[size]; 
    float *ycomp = new float[size]; 
    float *zcomp = new float[size]; 

    fread(xcomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(ycomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(zcomp, sizeof(float), size, in);  // read the data (vectors) 

    for (int c=0; c<size; c++) {
      vectors[t][i][c*3]   =xcomp[c]; 
      vectors[t][i][c*3+1] =ycomp[c]; 
      vectors[t][i][c*3+2] =zcomp[c]; 
    }
    delete [] xcomp; 
    delete [] ycomp; 
    delete [] zcomp; 

    // compute the length of each cell in the block 
    // will be used to determine the level which the block belongs to 
    sizes[i] = voxel_xsize_inBlock[t][i] = block_length[t][idx]/(float)block_dims[0]; 
    voxel_ysize_inBlock[t][i] = block_length[t][idx+1]/(float)block_dims[1]; 
    voxel_zsize_inBlock[t][i] = block_length[t][idx+2]/(float)block_dims[2]; 
  }
  int count = nb[t]; 
  sort_list(sizes, count); 
  printf(" there are %d different sizes\n", count); 
  num_levels[t] = count; // the total levels in this data 

  level_minB[t] = new float*[num_levels[t]];  // now determin the min/max corners of each level 
  level_maxB[t] = new float*[num_levels[t]]; 
  block_xsize_inLevel[t] = new float[num_levels[t]]; //  the physical block_length in each level 
  block_ysize_inLevel[t] = new float[num_levels[t]]; 
  block_zsize_inLevel[t] = new float[num_levels[t]]; 
  for (int i=0; i<num_levels[t]; i++) {
    level_minB[t][i] = new float[3]; 
    level_maxB[t][i] = new float[3]; 
    level_minB[t][i][0] = level_minB[t][i][1] = level_minB[t][i][2] = 999999999; 
    level_maxB[t][i][0] = level_maxB[t][i][1] = level_maxB[t][i][2] = -999999999; 
  }
  // write the level info for each block 
  // also find the block lengths in x y z for all levels, and 
  // the min max corners for all levels 

  for (int i=0; i<nb[t]; i++) {
    for (int j=0; j<num_levels[t]; j++) {
      if (voxel_xsize_inBlock[t][i] == sizes[j]) {
	int level = num_levels[t]-1-j; 
	block_level[t][i]= level; 
	block_xsize_inLevel[t][level] = block_length[t][i*3]; 
	block_ysize_inLevel[t][level] = block_length[t][i*3+1]; 
	block_zsize_inLevel[t][level] = block_length[t][i*3+2]; 
	int idx = i*3; 

	// update the level min corner 
	if (block_minB[t][idx]<level_minB[t][level][0]) 
	  level_minB[t][level][0] = block_minB[t][idx]; 
	if (block_minB[t][idx+1]<level_minB[t][level][1]) 
	  level_minB[t][level][1] = block_minB[t][idx+1]; 
	if (block_minB[t][idx+2]<level_minB[t][level][2]) 
	  level_minB[t][level][2] = block_minB[t][idx+2]; 

	// update the level max corner 
	if (block_maxB[t][idx]>level_maxB[t][level][0]) 
	  level_maxB[t][level][0] = block_maxB[t][idx]; 
	if (block_maxB[t][idx+1]>level_maxB[t][level][1]) 
	  level_maxB[t][level][1] = block_maxB[t][idx+1]; 
	if (block_maxB[t][idx+2]>level_maxB[t][level][2]) 
	  level_maxB[t][level][2] = block_maxB[t][idx+2]; 
      }
    }
  }


}


/////////////////////////////////////////////////////////////////////
//
//  since the grid is refined over time,  we need to match levels from all time
//  time steps 

void match_all_levels_over_time() {

  int cnt = 0; 

  for(int i=0; i<num_timesteps; i++) {
    level_mapping[i] = new int[num_levels[i]]; 
    cnt += num_levels[i]; 
  }
  float* sizes = new float[cnt]; 

  int idx = 0; 
  for (int i=0; i<num_timesteps; i++) 
    for (int j=0; j<num_levels[i]; j++) 
      sizes[idx++] = block_xsize_inLevel[i][j]; 
  
  sort_list(sizes, cnt); 
  printf(" there are total %d different levels \n", cnt); 
  
  total_levels = cnt; 
  
  g_block_xsize_inLevel = new float[total_levels]; 
  g_block_ysize_inLevel = new float[total_levels]; 
  g_block_zsize_inLevel = new float[total_levels]; 
  g_level_minB = new float*[total_levels]; 
  g_level_maxB = new float*[total_levels]; 

  for (int i=0; i<total_levels; i++) {
    g_level_minB[i] = new float[3]; 
    g_level_maxB[i] = new float[3]; 
  }
  
  // now match the levels 

  for (int i=0; i<num_timesteps; i++) 
    for (int j=0; j<num_levels[i]; j++) 
      for (int l = 0; l<total_levels; l++) {
	if (block_xsize_inLevel[i][j] == sizes[l])  {

	  level_mapping[i][j] = total_levels-1-l; 
	  g_block_xsize_inLevel[total_levels-1-l] = block_xsize_inLevel[i][j]; 
	  g_block_ysize_inLevel[total_levels-1-l] = block_ysize_inLevel[i][j]; 
	  g_block_zsize_inLevel[total_levels-1-l] = block_zsize_inLevel[i][j]; 

	  printf(" time step %d level %d is mapped to global level %d\n", i, j, level_mapping[i][j]); 
	  break; 
	}
      }

  // now we need to fix the levels to make them consistent 

  for (int i=0; i<num_timesteps; i++) 
    for (int j=0; j<nb[i]; j++) {
      int lvl = block_level[i][j]; 
      block_level[i][j] = level_mapping[i][lvl]; 
    }

  for (int i=0; i<total_levels; i++) {
    g_level_minB[i][0] = g_level_minB[i][1] = g_level_minB[i][2] = 999999999; 
    g_level_maxB[i][0] = g_level_maxB[i][1] = g_level_maxB[i][2] = -999999999; 
  }

  for (int i=0; i <num_timesteps; i++) {
    for (int j=0; j<num_levels[i]; j++) {
      int level = level_mapping[i][j]; 
      if (level_minB[i][j][0] < g_level_minB[level][0])
	g_level_minB[level][0] = level_minB[i][j][0]; 
      if (level_minB[i][j][1] < g_level_minB[level][1])
	g_level_minB[level][1] = level_minB[i][j][1]; 
      if (level_minB[i][j][2] < g_level_minB[level][2])
	g_level_minB[level][2] = level_minB[i][j][2]; 

      if (level_maxB[i][j][0] > g_level_maxB[level][0])
	g_level_maxB[level][0] = level_maxB[i][j][0]; 
      if (level_maxB[i][j][1] > g_level_maxB[level][1])
	g_level_maxB[level][1] = level_maxB[i][j][1]; 
      if (level_maxB[i][j][2] > g_level_maxB[level][2])
	g_level_maxB[level][2] = level_maxB[i][j][2]; 

    }
  }

}


/////////////////////////////////////////////////////////////////////////
//
//   Load a time-varying AMR data set 
// 

void LoadTimeVaryingAMR(char* fname, float* min, float* max)  {

  FILE *fIn; 
  char filename[300]; 
  float pmin[3], pmax[3]; 
  
  fIn = fopen(fname, "r"); 
  assert(fIn !=NULL); 
  fscanf(fIn, "%d", &num_timesteps); 

  // allocate memory for every time steps 

  nb = new int [num_timesteps];        // number of blocks in each time steps

  num_levels = new int[num_timesteps]; // no of levels
  level_mapping = new int*[num_timesteps]; 

  block_center = new float*[num_timesteps]; // block_center of blocks in the domain 
  block_length = new float*[num_timesteps]; // physical dimensions 
  block_index = new int*[num_timesteps]; 

  block_minB = new float*[num_timesteps]; 
  block_maxB = new float*[num_timesteps]; 

  voxel_xsize_inBlock = new float*[num_timesteps];  
  voxel_ysize_inBlock = new float*[num_timesteps]; 
  voxel_zsize_inBlock = new float*[num_timesteps]; 

  block_xsize_inLevel = new float*[num_timesteps]; 
  block_ysize_inLevel = new float*[num_timesteps]; 
  block_zsize_inLevel = new float*[num_timesteps]; 

  block_level = new int*[num_timesteps]; 

  level_minB = new float**[num_timesteps]; 
  level_maxB = new float**[num_timesteps]; 

  vectors = new float**[num_timesteps]; // the actual data 

  for (int i = 0; i<num_timesteps; i++) {
    fscanf(fIn, "%s", filename); 
    printf(" to read %s ....\n", filename); 
    LoadAMR(filename, pmin, pmax, i); 
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

  match_all_levels_over_time(); 
}


/////////////////////////////////////////////////////////////////////

main(int argc, char* argv[]) 
{

  float min[3], max[3]; 
  LoadTimeVaryingAMR(argv[1], min, max); 

  domain_lengths[0] = max[0]-min[0]; 
  domain_lengths[1] = max[1]-min[1]; 
  domain_lengths[2] = max[2]-min[2]; 
  domain_center[0] = (min[0]+max[0])/2.0; 
  domain_center[1] = (min[1]+max[1])/2.0; 
  domain_center[2] = (min[2]+max[2])/2.0; 

  latticeAMR = new LatticeAMR(domain_lengths[0], domain_lengths[1], domain_lengths[2], 
			      num_timesteps, total_levels); 


  for (int i=0; i<total_levels; i++)  {
    printf(" level size: [%f %f %f] min corner [%f %f %f] max corner [%f %f %f]\n", 
	   g_block_xsize_inLevel[i], g_block_ysize_inLevel[i], g_block_zsize_inLevel[i], 
	   g_level_minB[i][0], g_level_minB[i][1], g_level_minB[i][2], 
	   g_level_maxB[i][0], g_level_maxB[i][1], g_level_maxB[i][2]); 

    latticeAMR->CreateLevel(i, g_block_xsize_inLevel[i], g_block_ysize_inLevel[i], 
			    g_block_zsize_inLevel[i],
			    block_dims[0], block_dims[1], block_dims[2],
			    g_level_minB[i][0], g_level_maxB[i][0], g_level_minB[i][1], 
			    g_level_maxB[i][1], g_level_minB[i][2], g_level_maxB[i][2], 
			    0, num_timesteps-1); 

  }

  for (int t=0; t<num_timesteps; t++) {
    for (int i=0; i<nb[t]; i++) {
      latticeAMR->CheckIn(block_level[t][i], block_center[t][i*3], block_center[t][i*3+1], block_center[t][i*3+2], 
			  (float) t, vectors[t][i]); 
    }
  }

  latticeAMR->CompleteLevels(2); // this call is very important. It finishes up all the 
                                // book-keeping and complete the lattice setup

  vb_list = latticeAMR->GetBoundsList(npart); 

  // now initialize one osuflow object for each partition 

  osuflow_list = new OSUFlow*[npart]; 

  for(int i=0; i<npart; i++) {

    osuflow_list[i] = new OSUFlow(); 
    float **data = latticeAMR->GetDataPtr(i); 
    if (data == NULL) {
      printf(" panic\n"); 
      exit(1); 
    }
    float minB[3], maxB[3]; 
    float min_t, max_t; 
    minB[0] = vb_list[i].xmin; maxB[0] = vb_list[i].xmax;     
    minB[1] = vb_list[i].ymin; maxB[1] = vb_list[i].ymax;
    minB[2] = vb_list[i].zmin; maxB[2] = vb_list[i].zmax;
    min_t = vb_list[i].tmin; max_t = vb_list[i].tmax; 

    // dims[0/1/2] are the x y z dat dimensions 
    // minB and maxB are the physical space min/max corners 

    osuflow_list[i]->CreateTimeVaryingFlowField(data, block_dims[0], block_dims[1], 
						block_dims[2], minB, maxB, 
						min_t, max_t);     

  }

}
