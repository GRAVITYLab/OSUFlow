//------------------------------------------------------------------------------
//
// serial test draw
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

#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#ifdef LINUX
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

#ifdef MAC_OSX_10_4
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#endif

#include "OSUFlow.h"
#include "calc_subvolume.h"
#include "Lattice4D.h"

// defines related to drawing
#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

// drawing state
int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0;
float scale_size = 1;
int xform_mode = 0; 
bool toggle_bounds = true;

// drawing data
VECTOR4 *pt; // points in everyone's traces
int *npt; // everyone's number of points in their traces
int tot_ntrace; // total number of everyone's traces

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void PrintSeeds(int nblocks);
void ComputePathlines(int block_num);
void ComputeStreamlines(int block_num);
void GatherFieldlines();
void DrawFieldlines();
void Cleanup();
void draw_bounds(float *from, float *to);
void draw_cube(float r, float g, float b);
void display();
void timer(int val);
void mymouse(int button, int state, int x, int y);
void mymotion(int x, int y);
void mykey(unsigned char key, int x, int y);
void idle();
void Run();

// globals
static char filename[256]; // dataset file name
float size[3]; // spatial domain size
static int tsize; // temporal domain size
int *NumSeeds; // number of seeds
int *SizeSeeds; // size of seeds list (bytes)
VECTOR4 **Seeds; // list of seeds lists
VECTOR3 *seeds; // one temporary list of (3d) seeds
OSUFlow **osuflow; // one flow object for each block
list<vtListTimeSeedTrace*> *sl_list; // pathlines list
int nspart; // global total number of spatial blocks
int ntpart; // global total number of temporal blocks
int nblocks; // my number of blocks
Lattice4D* lat; // lattice
int tf; // max number of traces per block
int pf; // max number of points per trace
int max_rounds; // max number of rounds

// debug
#define MAX_RENDER_SEEDS 1000
VECTOR3 render_seeds[MAX_RENDER_SEEDS]; // seeds for rendering
int num_render_seeds = 0; // number of seeds for rendering
//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  char buf[256];

  GetArgs(argc, argv);
  Init();
  Run();

#ifdef GRAPHICS

  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  sprintf(buf, "Pathlines");
  glutCreateWindow(buf); 
  glutDisplayFunc(display); 
  glutIdleFunc(idle); 
  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 

#endif

  Cleanup();

}
//-----------------------------------------------------------------------
//
// Run
//
void Run() {

  int min_t, max_t; // subdomain temporal bounds
  float from[3], to[3]; // subdomain spatial bounds
  VECTOR3 min_s, max_s; // subdomain spatial bounds in VEC3 format
  int i, j, k;

  // init all blocks
  for (i = 0; i < nblocks; i++) {

    osuflow[i] = new OSUFlow;
    lat->GetVB(i, from, to, &min_t, &max_t);

    // init seeds for blocks at t = initial time
    if (min_t == 0) {

      osuflow[i]->SetRandomSeedPoints(from, to, tf); 
      seeds = osuflow[i]->GetSeeds(NumSeeds[i]); 

      while (SizeSeeds[i] < NumSeeds[i] * sizeof(VECTOR4)) {
	Seeds[i] = (VECTOR4 *)realloc(Seeds[i], SizeSeeds[i] * 2);
	assert(Seeds[i] != NULL);
	SizeSeeds[i] *= 2;
      }

      for (k = 0; k < NumSeeds[i]; k++)
	Seeds[i][k].Set(seeds[k][0], seeds[k][1], seeds[k][2], min_t);

    } // init seeds

  } // init all blocks

  // for all rounds
  for (j = 0; j < max_rounds; j++) {

#ifdef DEBUG
    PrintSeeds(nblocks);
#endif

    // for all blocks
    for (i = 0; i < nblocks; i++) {

      // load block for first round only
      if (j == 0) {
	lat->GetVB(i, from, to, &min_t, &max_t);
	min_s[0] = from[0]; min_s[1] = from[1]; min_s[2] = from[2];
	max_s[0] = to[0]; max_s[1] = to[1]; max_s[2] = to[2];
	osuflow[i]->LoadData(filename, false, min_s, max_s, min_t, max_t); 
	osuflow[i]->ScaleField(10.0); // improves visibility
      }

      if (tsize > 1)
	ComputePathlines(i);
      else
	ComputeStreamlines(i);

    } // for all blocks

    lat->SerExchangeNeighbors(Seeds, SizeSeeds, NumSeeds);

  } // for all rounds

  // gather fieldlines for rendering
  GatherFieldlines();

}
//-----------------------------------------------------------------------
//
// ComputeStreamlines
//
// block_num: local block number (0 to nblocks-1)
// not global partition number
//
void ComputeStreamlines(int block_num) {

  list<vtListSeedTrace*> list3; // 3D list of traces
  std::list<VECTOR3*>::iterator pt_iter3; // 3D iterator over pts in one trace
  std::list<vtListSeedTrace*>::iterator trace_iter3; // 3D iter. over traces
  VECTOR3 *Seeds3; // 3D seeds in current block
  VECTOR3 p3; // 3D current point
  VECTOR4 *p; // 4D current point
  vtListTimeSeedTrace *trace; // 4D single trace
  int i;

  if (NumSeeds[block_num]) {

    // make VECTOR3s (temporary)
    assert((Seeds3 = (VECTOR3 *)malloc(NumSeeds[block_num] * sizeof(VECTOR3)))
	   != NULL);
    for (i = 0; i < NumSeeds[block_num]; i++) {
      Seeds3[i][0]= Seeds[block_num][i][0];
      Seeds3[i][1]= Seeds[block_num][i][1];
      Seeds3[i][2]= Seeds[block_num][i][2];
    }
	  
    // perform the integration
    // todo: integrate in both directions
    osuflow[block_num]->SetIntegrationParams(1, 5); 
    osuflow[block_num]->GenStreamLines(Seeds3, FORWARD_DIR, 
				       NumSeeds[block_num], pf, list3); 

    // copy each 3D trace to a 4D trace and then to the streamline list
    // post end point of each trace to the send list
    for (trace_iter3 = list3.begin(); trace_iter3 != list3.end(); 
	 trace_iter3++) {

      if (!(*trace_iter3)->size())
	continue;

      trace = new vtListTimeSeedTrace;
      for (pt_iter3 = (*trace_iter3)->begin(); pt_iter3 != 
	     (*trace_iter3)->end(); pt_iter3++) {
	p3 = **pt_iter3;
	p = new VECTOR4;
	p->Set(p3[0], p3[1], p3[2], 0.0f);
	trace->push_back(p);

      }

      lat->PostPoint(block_num, *p); // last point only
      sl_list[block_num].push_back(trace); // for later rendering

    }

  }
  
}
//-----------------------------------------------------------------------
//
// ComputePathlines
//
// block_num: local block number (0 to nblocks-1)
// not global partition number
//
void ComputePathlines(int block_num) {
  list<vtListTimeSeedTrace*> list; // list of traces
  std::list<VECTOR4*>::iterator pt_iter; // iterator over pts in one trace
  std::list<vtListTimeSeedTrace*>::iterator trace_iter; // iter. over traces
  VECTOR4 p; // current point
  int i;

  if (NumSeeds[block_num]) {

    // perform the integration
    // todo: integrate in both directions
    osuflow[block_num]->SetIntegrationParams(1, 5); 
    osuflow[block_num]->GenPathLines(Seeds[block_num], list, FORWARD, 
				       NumSeeds[block_num], pf); 

    // copy each trace to the streamline list for later rendering
    // post end point of each trace to the send list
    for (trace_iter = list.begin(); trace_iter != list.end(); 
	 trace_iter++) {

      if (!(*trace_iter)->size())
	continue;

      // get the end point
      pt_iter = (*trace_iter)->end();
      pt_iter--;
      p = **pt_iter;

      lat->PostPoint(block_num, p); // last point only
      sl_list[block_num].push_back(*trace_iter); // for later rendering

    }

  }
  
}
//-----------------------------------------------------------------------
//
// GatherFieldlines
//
// gathers all fieldlines for rendering
//
void GatherFieldlines() {

  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  int i, j, k;
  int tot_npts = 0;

  // compute number of traces and points
  for (i = 0; i < nblocks; i++) {
    tot_ntrace += sl_list[i].size();
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++)
      tot_npts += (*trace_iter)->size();
  }

  // allocate rendering data
  assert((npt = new int[tot_ntrace]) != NULL);
  assert((pt = new VECTOR4[tot_npts]) != NULL); // points in everyones traces

  // compute number of points in each trace and collect the points
  j = 0;
  k = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
         pt_iter++)
	pt[k++] = **pt_iter;
      npt[j++] = (*trace_iter)->size();
    }
  }

}
//-----------------------------------------------------------------------
//
// GetArgs
//
// gets command line args
//
void GetArgs(int argc, char *argv[]) {

  VECTOR3 minLen, maxLen; // spatial data bounds
  int minTime, maxTime; // time data bounds

  assert(argc >= 8);

  strncpy(filename,argv[1],sizeof(filename));

  // hard code the minimum corner to be 0,0,0,0
  // need to allow for variable data origin in the future
  minLen[0] = minLen[1] = minLen[2] = 0.0;
  minTime = 0;

  maxLen[0] = atof(argv[2]) - 1.0f;
  maxLen[1] = atof(argv[3]) - 1.0f;
  maxLen[2] = atof(argv[4]) - 1.0f;
  maxTime   = (int)(atof(argv[5]) - 1.0f);

  size[0] = maxLen[0] - minLen[0] + 1.0f; // data sizes, space and time
  size[1] = maxLen[1] - minLen[1] + 1.0f;
  size[2] = maxLen[2] - minLen[2] + 1.0f;
  tsize   = (int)(maxTime - minTime + 1.0f);

  nspart = atoi(argv[6]); // total space partitions
  ntpart = atoi(argv[7]); // total time partitions

  tf = atoi(argv[8]); // traces per block
  pf = atoi(argv[9]); // points per trace
  max_rounds = atoi(argv[10]); // rounds

}
//-----------------------------------------------------------------------
//
// Init
//
// inits the app
//
void Init() {

  int nt; // max total number of traces
  int np; // max total number of points
  int ghost = 1; // number of ghost cells per spatial edge
  int i;

  // init lattice and osuflow
  lat = new Lattice4D((int)size[0], (int)size[1], (int)size[2], tsize, ghost, 
		      nspart, ntpart);
  lat->InitSeedLists(); 
  nblocks = nspart * ntpart;
  assert((osuflow = new OSUFlow*[nblocks]) != NULL);
  for (i = 0; i < nblocks; i++)
    osuflow[i] = NULL;

  // init seeds
  Seeds = (VECTOR4 **)malloc(nblocks * sizeof(VECTOR4 *));
  assert(Seeds != NULL);
  NumSeeds = (int *)malloc(nblocks * sizeof(int));
  assert(NumSeeds != NULL);
  SizeSeeds = (int *)malloc(nblocks * sizeof(int));
  assert(SizeSeeds != NULL);
  for (i = 0; i < nblocks; i++) {
    Seeds[i] = (VECTOR4 *)malloc(sizeof(VECTOR4));
    assert(Seeds[i] != NULL);
    SizeSeeds[i] = sizeof(VECTOR4);
  }

  // allocate streamline list for each block
  sl_list = new list<vtListTimeSeedTrace*>[nspart * ntpart];

  // print some of the args
#ifdef DEBUG
  fprintf(stderr,"Volume size: X %.3lf Y %.3lf Z %.3lf t %d\n",
	  size[0], size[1], size[2], tsize);
  fprintf(stderr, "Number of compute rounds: %d\n", max_rounds);
#endif

}
//-----------------------------------------------------------------------
//
// Cleanup
//
// frees memory and such
//
void Cleanup() {

  int i;

  delete [] pt;
  delete [] npt;
  delete [] sl_list;

  for (i = 0; i < nblocks; i++) {
    free(Seeds[i]);
    if (osuflow[i] != NULL)
      delete osuflow[i];
  }
  free(NumSeeds);
  free(Seeds);

  delete [] osuflow;
  delete lat;

}
//-----------------------------------------------------------------------
//
// PrintSeeds
//
// prints the current seeds
// (debug)
//
void PrintSeeds(int nblocks) {

  int i, j;

  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < NumSeeds[i]; j++) {
      fprintf(stderr,"Block %d Seed %d: %.3f\t%.3f\t%.3f\t%.3f\n", i, j, 
	      Seeds[i][j][0], Seeds[i][j][1], Seeds[i][j][2], Seeds[i][j][3]);
      render_seeds[num_render_seeds][0] = Seeds[i][j][0];
      render_seeds[num_render_seeds][1] = Seeds[i][j][1];
      render_seeds[num_render_seeds][2] = Seeds[i][j][2];
      render_seeds[num_render_seeds][3] = Seeds[i][j][3];
      if (num_render_seeds < MAX_RENDER_SEEDS - 2)
	num_render_seeds++;
    }
  }

  fprintf(stderr,"\n");


}
//--------------------------------------------------------------------------

#ifdef GRAPHICS

//
void draw_bounds(float *from, float *to) {

  float xmin = from[0];
  float ymin = from[1];
  float zmin = from[2];
  float xmax = to[0];
  float ymax = to[1];
  float zmax = to[2];

  glColor3f(1,0,0); 
  glBegin(GL_LINES); 

  glVertex3f(xmin, ymin, zmin); glVertex3f(xmax, ymin, zmin); 
  glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymax, zmin); 
  glVertex3f(xmax, ymax, zmin); glVertex3f(xmin, ymax, zmin); 
  glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, xmin, zmin); 

  glVertex3f(xmin, ymin, zmax); glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax); 
  glVertex3f(xmax, ymax, zmax); glVertex3f(xmin, ymax, zmax); 
  glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, xmin, zmax); 

  glVertex3f(xmin, ymin, zmin); glVertex3f(xmin, ymin, zmax); 
  glVertex3f(xmin, ymin, zmax); glVertex3f(xmin, ymax, zmax); 
  glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, ymax, zmin); 
  glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, ymin, zmin); 

  glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax); 
  glVertex3f(xmax, ymax, zmax); glVertex3f(xmax, ymax, zmin); 
  glVertex3f(xmax, ymax, zmin); glVertex3f(xmax, ymin, zmin); 

  glEnd(); 

}
//--------------------------------------------------------------------------
//
void DrawFieldlines() {
  
  float from[3], to[3]; // subdomain spatial bounds
  int min_t, max_t; // subdomain temporal bounds (not used)
  static int step = 0; // timestep number
  static int frame_num = 0; // number of frames at this timestep
  static int frames_per_step; // hold each time step for this many frames
  int i, j, k;

  // compute hold time
  frames_per_step = 2000 / tsize;

  glPushMatrix(); 

  glScalef(1.0f / (float)size[0], 1.0f / (float)size[0], 1.0f / (float)size[0]);
  glTranslatef(-size[0] / 2.0f, -size[1] / 2.0f, -size[2] / 2.0f); 
  glColor3f(1.0, 1.0, 0.0); 

  k = 0;

  // traces
  for (int i = 0; i < tot_ntrace; i++) {

    glBegin(GL_LINE_STRIP); 
    for (j = 0; j < npt[i]; j++) {
      if (pt[k][3] <= step)
	glVertex3f(pt[k][0], pt[k][1], pt[k][2]);
      k++;
    }
    glEnd(); 

  }

  // bounds
  if (toggle_bounds) {
    for (i = 0; i < nspart * ntpart; i++) {
      lat->GetGlobalVB(i, from, to, &min_t, &max_t);
      draw_bounds(from, to);
    }
  }

#ifdef DEBUG
  // seeds
  glPointSize(5);
  glBegin(GL_POINTS);
  for (i = 0; i < num_render_seeds; i++)
    glVertex3f(render_seeds[i][0],render_seeds[i][1],render_seeds[i][2]);
  glEnd();
#endif

  glPopMatrix(); 

  if (frame_num == frames_per_step) {
    frame_num = 0;
    if (step == tsize - 1)
      step = 0;
    else
      step++;
  }
  else
    frame_num++;

}
//-------------------------------------------------------------------------
//
void draw_cube(float r, float g, float b) {

  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube

}
//--------------------------------------------------------------------------
//
void display() {

  glEnable(GL_DEPTH_TEST); 
  glClearColor(0.0, 0.0, 0.0, 1.0); 
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

  DrawFieldlines(); 

  glPushMatrix(); 
  glScalef(1.0, size[1] / size[0], size[2] / size[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 
  
  glBegin(GL_LINES); 
  glColor3f(1,0,0); // red x-axis
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0);
  glColor3f(0,1,0); // green y-axis
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1); // blue z-axis
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 

  glutSwapBuffers(); 

}
//--------------------------------------------------------------------------
//
void timer(int val) {

  glutTimerFunc(10, timer, 0); 

}
//--------------------------------------------------------------------------
//
void idle() {

  glutPostRedisplay();

}
//--------------------------------------------------------------------------
//
void mymouse(int button, int state, int x, int y) {

  if (state == GLUT_DOWN) {

    press_x = x; press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE;
    else if (button == GLUT_RIGHT_BUTTON)
      xform_mode = XFORM_SCALE; 

  }

  else if (state == GLUT_UP)
    xform_mode = XFORM_NONE;

}
//--------------------------------------------------------------------------
//
void mymotion(int x, int y) {

  if (xform_mode==XFORM_ROTATE) {

    x_angle += (x - press_x)/5.0; 
    if (x_angle > 180)
      x_angle -= 360; 
    else if (x_angle <-180)
      x_angle += 360; 
    press_x = x; 
	   
    y_angle += (y - press_y)/5.0; 
    if (y_angle > 180)
      y_angle -= 360; 
    else if (y_angle <-180)
      y_angle += 360; 
    press_y = y; 

  }

  else if (xform_mode == XFORM_SCALE){

    float old_size = scale_size;
    scale_size *= (1+ (y - press_y)/60.0); 
    if (scale_size <0)
      scale_size = old_size; 
    press_y = y; 

  }

  glutPostRedisplay(); 

}
//--------------------------------------------------------------------------
//
void mykey(unsigned char key, int x, int y) {

  switch(key) {
  case 'q':
    exit(1);
    break; 
  case 'b':
    toggle_bounds = !toggle_bounds;
    break;
  default:
    break;
  }

}
//--------------------------------------------------------------------------

#endif
