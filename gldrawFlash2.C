
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

int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0; 
float scale_size = 1;

int xform_mode = 0; 

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

///////////////////////////////
//
// AMR book-keeping variables 
//
int nb;                // number of blocks
int num_levels;        // the total number of level 
int dims[3];           // resolutions of each block
float *mins, *maxs;    // two global min/max corners 
float *centers;        // center of each block 
float *lengths;        // block dimensions
int *indices;          // block id 
float **vectors;        // cell center block data 
float len[3];          // global physical x y z lengths
float center[3];       // center of the domain 
float *element_xsize;  // voxel x length in each block 
float *element_ysize;  // voxel y length
float *element_zsize;  // voxel z length
float **level_minB;    // the coordinates of min corner of each level
float **level_maxB;    // the coordinates of max corner of each level 
float *level_xsize, *level_ysize, *level_zsize; // the physical length of each level 
int *its_level;        // the level which each block belongs to 

LatticeAMR * latticeAMR; 
volume_bounds_type_f* vb_list; 

OSUFlow **osuflow_list; 

list<vtListSeedTrace*> *sl_list; 
VECTOR3 **osuflow_seeds; 
int *osuflow_num_seeds; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 

int num_frames = 20; 
int current_frame = 0; 

int draw_level=0; 

bool to_draw_bbx = false; 
bool to_draw_center = false; 
bool to_draw_vector = false; 

int total_seeds = 5000; 
int npart; 
int first_frame = 1; 
////////////////////////////////////////////////////////

void compute_streamlines() {
  
  float from[3], to[3]; 

      for (int i=0; i<npart; i++) {

    VECTOR3 minB, maxB; 
    from[0] = vb_list[i].xmin;  
    from[1] = vb_list[i].ymin;     
    from[2] = vb_list[i].zmin; 
    to[0] = vb_list[i].xmax;  
    to[1] = vb_list[i].ymax;     
    to[2] = vb_list[i].zmax; 

    printf(" from %f %f %f to %f %f %f \n", from[0], from[1], from[2], 
	   to[0], to[1], to[2]); 
    osuflow_list[i]->SetRandomSeedPoints(from, to, 1); // set range for seed locations
    int num; 
    osuflow_seeds[i] = osuflow_list[i]->GetSeeds(num); 
    osuflow_num_seeds[i] = num; 
    sl_list[i].clear(); 
  }
  
  // Now begin to perform particle tracing in all subdomains
  bool has_seeds = true;      // initially we always have seeds
  int num_seeds_left = 20*npart; 

    for (int i=0; i<npart; i++) {
  //    for (int i=0; i<40; i++) {
      list<vtListSeedTrace*> list; 
      osuflow_list[i]->SetIntegrationParams(1, 5); 
      osuflow_list[i]->GenStreamLines(osuflow_seeds[i], FORWARD_DIR, 
      				      osuflow_num_seeds[i], 50, list); 
      printf("domain %d done integrations", i); 
      printf(" %d streamlines. \n", list.size()); 

      std::list<vtListSeedTrace*>::iterator pIter; 

      //looping through and store the trace points into sl_list[i] 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
        vtListSeedTrace *trace = *pIter; 
	printf(" trace length %d   ", trace->size()); 
	sl_list[i].push_back(trace); 
      }
  }
}



void draw_streamlines() {

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  std::list<vtListSeedTrace*>::iterator pIter; 

  for (int i=0; i<npart; i++) {   // looping through all subdomains 

    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 
    /*
    glColor3f(1,1,1); 
    glBegin(GL_POINTS); 
    for (int j=0; j<nSeeds; j++) {
      printf(" Seed: %f %f %f \n", seeds[j][0], seeds[j][1], seeds[j][2]); 
      glVertex3f(seeds[j][0], seeds[j][1], seeds[j][2]); 
    }
    glEnd(); 
    */

    pIter = sl_list[i].begin(); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      vtListSeedTrace *trace = *pIter; 
      std::list<VECTOR3*>::iterator pnIter; 
      pnIter = trace->begin(); 
      glColor3f(1,1,0); 
      glBegin(GL_LINE_STRIP); 
      //      glBegin(GL_POINTS); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
	//	printf(" %f %f %f \n", p[0], p[1], p[2]); 
	glVertex3f(p[0], p[1], p[2]); 
      }
      glEnd(); 
    }
  }
  glPopMatrix(); 
}

void animate_streamlines() {

  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  std::list<VECTOR3*>::iterator pnIter; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  float min_time = current_frame/(float) num_frames; 
  float max_time = (current_frame+1)/(float) num_frames; 

  glColor3f(1,1,0); 
  for (int i=0; i<npart; i++) {
    pIter = sl_list[i].begin(); 
    glBegin(GL_POINTS); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      trace = *pIter; 
      pnIter = trace->begin(); 
      int cnt=0; 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
	float ratio = cnt/(float)trace->size(); 
	if (ratio>= min_time && ratio < max_time) {
	  float x = p[0]; 
	  float y = p[1]; 
	  float z = p[2]; 
	  glVertex3f(x,y,z); 
	}
	cnt++; 
      }
    }
    glEnd(); 
  }
  glPopMatrix(); 
  current_frame = (current_frame+1) % num_frames; 
}

////////////////////////////////////////////// 

void draw_cube(float r, float g, float b)
{
  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube 
}

///////////////////////////////////////////////

void draw_bbx(int i) {


  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  float xsize, ysize, zsize; 
  if (its_level[i]%3 == 0) glColor3f(1,0,0); 
  if (its_level[i]%3 == 1) glColor3f(0,1,0); 
  if (its_level[i]%3 == 2) glColor3f(0,0,1); 

  /* 
  // test the centers
  glPushMatrix(); 
  glTranslatef(centers[i*3], centers[i*3+1], centers[i*3+2]); 
  glScalef(lengths[i*3], lengths[i*3+1], lengths[i*3+2]); 

  glutWireCube(1.0); 
  glPopMatrix(); 
  */

  // test the volume bounds 

  int rank = latticeAMR->GetRank(centers[i*3], centers[i*3+1], centers[i*3+2], 0); 

  float xmin = vb_list[rank].xmin; 
  float ymin = vb_list[rank].ymin; 
  float zmin = vb_list[rank].zmin; 
  float xmax = vb_list[rank].xmax; 
  float ymax = vb_list[rank].ymax; 
  float zmax = vb_list[rank].zmax; 
  
  glBegin(GL_LINE_STRIP); 
  glVertex3f(xmin, ymin, zmin); 
  glVertex3f(xmax, ymin, zmin); 
  glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmin, ymin, zmax); 
  glVertex3f(xmin, ymin, zmin); 
  glEnd(); 

  glPopMatrix(); 

}

void draw_vector(int i) 
{
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glBegin(GL_LINES); 
  glColor3f(1,0,0); 
  glVertex3f(centers[i*3], centers[i*3+1], centers[i*3+2]); 
  glColor3f(0,0,1); 
  glVertex3f(centers[i*3]+ vectors[i][0]*0.3, 
	     centers[i*3+1]+vectors[i][1]*0.3, 
	     centers[i*3+2]+vectors[i][2]*0.3); 
  glEnd(); 

  glPopMatrix(); 

}

/////////////////////////////////////////////////////
void draw_center(int i) 
{
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glPointSize(2); 

  glColor3f(1,1,1); 

  int level = latticeAMR->GetFinestLevel(centers[i*3], centers[i*3+1], centers[i*3+2],0);

  if (level %3 == 0) glColor3f(1,0,0); 
  if (level %3 == 1) glColor3f(0,1,0); 
  if (level %3 == 2) glColor3f(0,0,1); 

  glBegin(GL_POINTS); 
  glVertex3f(centers[i*3], centers[i*3+1], centers[i*3+2]); 
  glEnd(); 

  glPopMatrix(); 

}
///////////////////////////////////////////////////////////

void mymouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    press_x = x; press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE; 
	 else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode = XFORM_SCALE; 
  }
  else if (state == GLUT_UP) {
	  xform_mode = XFORM_NONE; 
  }
}

/////////////////////////////////////////////////////////

void mymotion(int x, int y)
{
    if (xform_mode==XFORM_ROTATE) {
      x_angle += (x - press_x)/5.0; 
      if (x_angle > 180) x_angle -= 360; 
      else if (x_angle <-180) x_angle += 360; 
      press_x = x; 
	   
      y_angle += (y - press_y)/5.0; 
      if (y_angle > 180) y_angle -= 360; 
      else if (y_angle <-180) y_angle += 360; 
      press_y = y; 
    }
	else if (xform_mode == XFORM_SCALE){
      float old_size = scale_size;
      scale_size *= (1+ (y - press_y)/60.0); 
      if (scale_size <0) scale_size = old_size; 
      press_y = y; 
    }
	glutPostRedisplay(); 
}


//////////////////////////////////////////////////////

void display()
{
  glEnable(GL_DEPTH_TEST); 
  glClearColor(0,0,0,1); 
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 
  gluPerspective(60, 1, .1, 100); 

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(0,0,5,0,0,0,0,1,0); 

  glRotatef(x_angle, 0, 1,0); 
  glRotatef(y_angle, 1,0,0); 
  glScalef(scale_size, scale_size, scale_size); 

  
  //  printf(" len %f %f %f\n", len[0], len[1], len[2]); 
  glPushMatrix(); 
  //  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 

  if (toggle_draw_streamlines == true)
    draw_streamlines(); 
  else if (toggle_animate_streamlines == true)
    animate_streamlines(); 

  for (int i=0; i<nb; i++) {
    if (its_level[i]>=draw_level) {
      if (to_draw_bbx) draw_bbx(i); 
      if (to_draw_center) draw_center(i); 
      if (to_draw_vector) draw_vector(i); 
    }
  }
    
  //  glColor3f(1,1,1); 
  glBegin(GL_LINES); 
  glColor3f(1,0,0); 
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0); 
  glColor3f(0,1, 0); 
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1); 
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 

  glutSwapBuffers(); 
}

void timer(int val) {
  if (toggle_animate_streamlines == true) {
    glutPostRedisplay(); 
  }
  glutTimerFunc(200, timer, 0); 
}

///////////////////////////////////////////////////////////////

void mykey(unsigned char key, int x, int y)
{
        switch(key) {
	case 'q': exit(1);
	  break; 
	case 's': compute_streamlines(); 
	  glutPostRedisplay(); 
	  break; 
	case 'l': 
	  draw_level++; 
	  draw_level = draw_level % (num_levels); 
	  break; 
	case 'b':
	  to_draw_bbx = !to_draw_bbx; 
	  break; 
	case 'c': 
	  to_draw_center = !to_draw_center; 
	  break; 
	case 'v': 
	  to_draw_vector = !to_draw_vector; 
	  break; 
	case 'd': 
	  toggle_draw_streamlines = !toggle_draw_streamlines; 
	  toggle_animate_streamlines = false; 
	  break; 
	case'a': 
	  toggle_animate_streamlines = !toggle_animate_streamlines; 
	  toggle_draw_streamlines = false; 
	  first_frame = 1; 
	}
	glutPostRedisplay(); 
}

void sort_list(float* list, int& cnt) 
{
  float* tmp_list = new float[cnt]; 
  for (int i=0; i<cnt; i++) { 
    tmp_list[i] = list[i]; 
  }
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
//   A routine to read the data from the Flash application generated 
//   by Mike Papka's group. This format is not the original Flash format
//   but just an in-house format for quick development 
// 
void LoadAMR(char* fname, float* min, float *max)
{

  float *sizes; 
  FILE * in = fopen(fname, "r"); 
  fread(&nb, sizeof(int), 1, in); // number of blocks
  printf(" %d blocks.\n", nb); 
  fread(dims, sizeof(int), 3, in); // dims per block, same for all level 
  printf(" block dims %d %d %d\n", dims[0], dims[1], dims[2]); //16x16x16 for example 

  int size = dims[0]*dims[1]*dims[2];  // number of elements in each block 

  its_level = new int[nb];      // the level that each block belongs to 
  indices = new int[nb];        // the index of each block (defined by application)
  centers = new float[nb*3];    // the center x y and z of each block 
  lengths = new float[nb*3];    // the physical lengths (in x y and z) of each block 
  mins = new float[nb*3];       // the min corner (min_x, min_y, min_z) of each block 
  maxs = new float[nb*3];       // the max corner (max_x, max_y, max_z) of each block 
  element_xsize = new float[nb]; // the physical length in x of each block  
  element_ysize = new float[nb]; // the physical length in y of each block 
  element_zsize = new float[nb]; // the physical length in z of each block 
  sizes = new float[nb];         // use one of the x y z legnths to determin block level 
  vectors = new float*[nb]; 

  for (int i=0; i<nb; i++) {

    int idx = i*3; 

    fread(&(indices[i]), sizeof(int), 1, in); 
    //    printf(" read block id: %d\n", indices[i]); 

    fread(&(centers[idx]), sizeof(float), 3, in); 
    //    printf(" block %d center: [%f %f %f]\n", indices[i], 
    //	   centers[idx], centers[idx+1], centers[idx+2]); 

    fread(&(lengths[idx]), sizeof(float), 3, in); 
    //    printf(" block %d length: [%f %f %f]\n", indices[i], 
    //	   lengths[idx], lengths[idx+1], lengths[idx+2]); 

    float bounds[6]; 
    fread(bounds, sizeof(float), 6, in);   // the two opposite corners of each block 
    mins[idx]   = bounds[0];  maxs[idx]   = bounds[1]; 
    mins[idx+1] = bounds[2];  maxs[idx+1] = bounds[3]; 
    mins[idx+2] = bounds[4];  maxs[idx+2] = bounds[5]; 

    if (i==0) {
      min[0] = mins[0]; min[1] = mins[1]; min[2] = mins[2]; 
      max[0] = maxs[0]; max[1] = maxs[1]; max[2] = maxs[2]; 
    }
    else {   // find the global min/max cornders 
      if (mins[idx]<min[0]) min[0] = mins[idx]; 
      if (mins[idx+1]<min[1]) min[1] = mins[idx+1];      
      if (mins[idx+2]<min[2]) min[2] = mins[idx+2]; 
      if (maxs[idx]>max[0]) max[0] = maxs[idx]; 
      if (maxs[idx+1]>max[1]) max[1] = maxs[idx+1]; 
      if (maxs[idx+2]>max[2]) max[2] = maxs[idx+2]; 
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
    sizes[i] = element_xsize[i] = lengths[idx]/(float)dims[0]; 
    element_ysize[i] = lengths[idx+1]/(float)dims[1]; 
    element_zsize[i] = lengths[idx+2]/(float)dims[2]; 
  }
  int count = nb; 
  sort_list(sizes, count); 
  printf(" there are %d different sizes\n", count); 
  num_levels = count; // the total levels in this data 


  level_minB = new float*[num_levels];  // now determin the min/max corners of each level 
  level_maxB = new float*[num_levels]; 
  level_xsize = new float[num_levels]; //  the physical lengths in each level 
  level_ysize = new float[num_levels]; 
  level_zsize = new float[num_levels]; 
  for (int i=0; i<num_levels; i++) {
    level_minB[i] = new float[3]; 
    level_maxB[i] = new float[3]; 
    level_minB[i][0] = level_minB[i][1] = level_minB[i][2] = 999999999; 
    level_maxB[i][0] = level_maxB[i][1] = level_maxB[i][2] = -999999999; 
  }
  // write the level info to each block 
  // also find the block length in x y z for all levels, as well 
  // as the min max spatial corners for all levels 

  for (int i=0; i<nb; i++) {
    for (int j=0; j<num_levels; j++) {
      if (element_xsize[i] == sizes[j]) {
	int level = num_levels-1-j; 
	its_level[i]= level; 
	level_xsize[level] = lengths[i*3]; 
	level_ysize[level] = lengths[i*3+1]; 
	level_zsize[level] = lengths[i*3+2]; 
	int idx = i*3; 

	// update the level min corner 
	if (mins[idx]<level_minB[level][0]) 
	  level_minB[level][0] = mins[idx]; 
	if (mins[idx+1]<level_minB[level][1]) 
	  level_minB[level][1] = mins[idx+1]; 
	if (mins[idx+2]<level_minB[level][2]) 
	  level_minB[level][2] = mins[idx+2]; 

	// update the level max corner 
	if (maxs[idx]>level_maxB[level][0]) 
	  level_maxB[level][0] = maxs[idx]; 
	if (maxs[idx+1]>level_maxB[level][1]) 
	  level_maxB[level][1] = maxs[idx+1]; 
	if (maxs[idx+2]>level_maxB[level][2]) 
	  level_maxB[level][2] = maxs[idx+2]; 
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////

main(int argc, char* argv[]) 
{

  float min[3], max[3]; 
  LoadAMR(argv[1], min, max); 

  len[0] = max[0]-min[0]; 
  len[1] = max[1]-min[1]; 
  len[2] = max[2]-min[2]; 
  center[0] = (min[0]+max[0])/2.0; 
  center[1] = (min[1]+max[1])/2.0; 
  center[2] = (min[2]+max[2])/2.0; 

  
  latticeAMR = new LatticeAMR(len[0], len[1], len[2], 1, num_levels); 
  for (int i=0; i<num_levels; i++)  {
    printf(" level size: [%f %f %f] min corner [%f %f %f] max corner [%f %f %f]\n", 
	   level_xsize[i], level_ysize[i], level_zsize[i], 
	   level_minB[i][0], level_minB[i][1], level_minB[i][2], 
	   level_maxB[i][0], level_maxB[i][1], level_maxB[i][2]); 

    latticeAMR->CreateLevel(i, level_xsize[i], level_ysize[i], level_zsize[i], 
			    dims[0], dims[1], dims[2], 
			    level_minB[i][0], level_maxB[i][0], level_minB[i][1], 
			    level_maxB[i][1], level_minB[i][2], level_maxB[i][2], 
			    0.0, 1.0); 

  }
  for (int i=0; i<nb; i++) {
    latticeAMR->CheckIn(its_level[i], centers[i*3], centers[i*3+1], centers[i*3+2], 0.0, 
			vectors[i]); 
  }

  latticeAMR->CompleteLevels(); // this call is very important. It finishes up all the 
                                // book-keeping there

  vb_list = latticeAMR->GetBoundsList(npart); 
  printf("nb = %d npart = %d \n", nb, npart);  // they should be equal 

  sl_list = new list<vtListSeedTrace*>[npart]; //one streamlines list for each subdomain 

  // now initialize the osuflow objects, one per block 
  osuflow_list = new OSUFlow*[npart]; 
  osuflow_seeds = new VECTOR3*[npart]; 
  osuflow_num_seeds = new int[npart]; 

  for (int i=0; i<npart; i++) {
    osuflow_list[i] = new OSUFlow(); 
    float **data = latticeAMR->GetDataPtr(i); // i is the rank 
    if (data == NULL) {
      printf(" panic\n"); 
      exit(1); 
    }
    float minB[3], maxB[3]; 
    minB[0] = vb_list[i].xmin; maxB[0] = vb_list[i].xmax;     
    minB[1] = vb_list[i].ymin; maxB[1] = vb_list[i].ymax;
    minB[2] = vb_list[i].zmin; maxB[2] = vb_list[i].zmax;
    // dims[0/1/2] are the x y z dat dimensions 
    // minB and maxB are the physical space min/max corners 
    printf(" %d =  dim %d %d %d ....\n", i, dims[0], dims[1], dims[2]); 
    osuflow_list[i]->CreateStaticFlowField(data[0], dims[0], dims[1], dims[2], 
					   minB, maxB); 
  }

  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  
  glutCreateWindow("Display streamlines"); 
  glutDisplayFunc(display); 

  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}
