//------------------------------------------------------------------------------
//
// mpi test draw
//
// Copyright (c) 2008 Han-Wei Shen and Tom Peterka
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
//------------------------------------------------------------------------------

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include <GL/gl.h>
#include <GL/glut.h>
#include <list>
#include <iterator>

#include "OSUFlow.h"
#include "calc_subvolume.h"
#include "Lattice.h"

#define MAX_RECV_ATTEMPTS 10
#define MAX_ITERATIONS 3

// defines related to drawing
#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

// globals related to drawing
int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0;
float scale_size = 1;
int xform_mode = 0; 
bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 

// function prototypes
void Config(int argc, char *argv[]);
void PrintSeeds(int nblocks);
void RecvToSeeds(int local_block, int global_block);
int EndTrace(list<vtListSeedTrace*> &list, VECTOR3 &p, int index);
void ComputeStreamlines();
void GatherStreamlines();
void draw_bounds(float xmin, float xmax, float ymin, float ymax, 
		 float zmin, float zmax);
void draw_streamlines();
void animate_streamlines();
void draw_cube(float r, float g, float b);
void display();
void timer(int val);
void mymouse(int button, int state, int x, int y);
void mymotion(int x, int y);
void mykey(unsigned char key, int x, int y);

// globals
static char filename[256]; // dataset file name
static VECTOR3 minLen, maxLen; // data bounds
VECTOR3 size; // domain size
static int bf; // number of blocks per rank
int *NumSeeds; // number of seeds
int *SizeSeeds; // size of seeds list (bytes)
VECTOR3 **Seeds; // seeds lists
OSUFlow **osuflow; // one flow object for each block
list<vtListSeedTrace*> *sl_list; // streamlines list
int npart; // global total number of blocks
volume_bounds_type *vb_list; // global subdomain volume bounds list
int *blocks; // local list of blocks
int nblocks; // local number of blocks
Lattice* lat; // lattice

//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  VECTOR3 minB, maxB; // subdomain bounds
  int nproc;  // mpi groupsize
  int rank; // mpi rank
  float from[3], to[3]; // seed points limits
  int i;
  int ghost; // number of ghost cells per side
  int init_seed_num = 10; // initial number of seeds
  char buf[256];

  // initialize mpi and the app
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Barrier(MPI_COMM_WORLD);
  Config(argc, argv);
  ghost = 1;

  // init the lattice
  lat = new Lattice(size[0], size[1], size[2], ghost, bf * nproc);
  vb_list = lat->GetBoundsList(); 
  lat->InitSeedLists(); 
  lat->RoundRobin_proc(nproc); 

  // create and populate the partitions
  if ((blocks = (int *)malloc(bf * sizeof(int))) == NULL)
    Error("Error: unable to allocate blocks\n");
  lat->GetPartitions(rank, blocks, nblocks);

  for (i = 0; i < nblocks; i++) {

    osuflow[i] = new OSUFlow();

    minB.Set(vb_list[blocks[i]].xmin, vb_list[blocks[i]].ymin, vb_list[blocks[i]].zmin);
    maxB.Set(vb_list[blocks[i]].xmax, vb_list[blocks[i]].ymax, vb_list[blocks[i]].zmax);

    fprintf(stderr,"Subdomain boundary: rank = %d i = %d global block = %d min = %.3lf %.3lf %.3lf max = %.3lf %.3lf %.3lf\n",rank,i,blocks[i],minB[0],minB[1],minB[2],maxB[0],maxB[1],maxB[2]);

    // read data
    osuflow[i]->ReadData(filename, true, minB, maxB, size); 
    if (rank == 0 && i == 0)
      fprintf(stderr,"read file %s\n", filename); 

    // init seeds
    from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
    to[0]   = maxB[0];   to[1]   = maxB[1];   to[2]   = maxB[2]; 
    osuflow[i]->SetRandomSeedPoints(from, to, init_seed_num); 
    Seeds[i] = osuflow[i]->GetSeeds(NumSeeds[i]); 
    SizeSeeds[i] = NumSeeds[i] * sizeof(VECTOR3);

    sl_list[i].clear();

  }

  // main loop for nondrawing procs
  if (rank > 0) {
    for (i = 0; i < MAX_ITERATIONS; i++) {
      ComputeStreamlines();
      GatherStreamlines();
    }
  }

  // main loop for drawing proc
  if (rank == 0) {

    glutInit(&argc, argv); 
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
    glutInitWindowSize(600,600); 
    sprintf(buf, "Streamlines");
    glutCreateWindow(buf); 
    glutDisplayFunc(display); 
    glutIdleFunc(ComputeStreamlines); 
    glutTimerFunc(10, timer, 0); 
    glutMouseFunc(mymouse); 
    glutMotionFunc(mymotion);
    glutKeyboardFunc(mykey); 
    glutMainLoop(); 

  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// ComputeStreamLines
//
void ComputeStreamlines() {

  std::list<vtListSeedTrace*>::iterator trace_iter; // iterator over seed traces
  vtListSeedTrace *trace; // single seed trace
  std::list<VECTOR3*>::iterator pt_iter; // iterator over points in one seed trace
  VECTOR3 p; // current point
  int neighbor; // neighbor's number (0-5)
  int ei, ej, ek; // neighbor's lattice position
  int i, j;

  // for all blocks, integrate points and send messages
  for(i = 0; i < nblocks; i++) {

//     // debug: print current seeds
//     fprintf(stderr, "Current seeds\n");
//     PrintSeeds(nblocks);

    if (NumSeeds[i]) {

      // declaration needs to be here to start with an empty list each time
      // need to learn how to clear a list w/o redeclaring it
      list<vtListSeedTrace*> list; // trace of seed points

      // perform the integration
      // todo: integrate in both directions
      osuflow[i]->SetIntegrationParams(1, 5); 
      osuflow[i]->GenStreamLines(Seeds[i], FORWARD_DIR, NumSeeds[i], 50, list); 

      // copy each trace to the streamline list for later rendering
      trace_iter = list.begin(); 
      for (; trace_iter != list.end(); trace_iter++) {
	trace = *trace_iter; 
	sl_list[i].push_back(trace); 
      }

      // redistribute boundary points to neighbors

      // for all boundary points
      for (trace_iter = list.begin(); trace_iter != list.end(); trace_iter++) {
	trace = *trace_iter;
	if (!trace->size())
	  continue;
	pt_iter = trace->end();
	pt_iter--;
	p = **pt_iter;

	// find which neighbor the point is in
	neighbor = lat->GetNeighbor(blocks[i], p[0], p[1], p[2], ei, ej, ek); 
	// post the point to the send list
	if (neighbor != -1)
	  lat->PostPoint(blocks[i], p, neighbor);

      } // for all boundary points

//       // debug
//       lat->PrintPost(blocks[i]);

      // send boundary list to neighbors
      lat->SendNeighbors(blocks[i], MPI_COMM_WORLD);

    } // if (NumSeeds[i])
      
  } // for all blocks

    // for all blocks, receive messages
  for (i = 0; i < nblocks; i++) {

    for (j = 0; j < MAX_RECV_ATTEMPTS; j++) {
      if (lat->ReceiveNeighbors(blocks[i], MPI_COMM_WORLD))
	break;
      usleep(100000);
    }

    if (j == MAX_RECV_ATTEMPTS)
      continue;

    // debug
    lat->PrintRecv(blocks[i]);

    // prepare for next iteration
    RecvToSeeds(i, blocks[i]);

  } // for all blocks

  // debug: print current seeds
  //   fprintf(stderr, "Final seeds\n");
  //   PrintSeeds(nblocks);

}
//-----------------------------------------------------------------------
//
// GatherStreamLines
//
// gathers all streamlines at the root for rendering
//
void GatherStreamLines() {

  MPI_rank rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {

  }

  if (rank > 0) {

  }

}
//-----------------------------------------------------------------------
//
// Config
//
// gets command line args
//
void Config(int argc, char *argv[]) {

  int rank;
  int nproc;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (argc < 6)
    Error("Error: Insufficient number of command line args\n");

  strncpy(filename,argv[1],sizeof(filename));

  // hard code the minimum corner to be 0,0,0
  // need to allow for variable data origin in the future
  minLen[0] = minLen[1] = minLen[2] = 0.0;

  maxLen[0] = atof(argv[2]) - 1.0f;
  maxLen[1] = atof(argv[3]) - 1.0f;
  maxLen[2] = atof(argv[4]) - 1.0f;

  size[0] = maxLen[0] - minLen[0] + 1.0f;
  size[1] = maxLen[1] - minLen[1] + 1.0f;
  size[2] = maxLen[2] - minLen[2] + 1.0f;

  bf = atoi(argv[5]);
  npart = bf * nproc;

  if (rank == 0)
    fprintf(stderr,"Volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
	   minLen[0], maxLen[0], minLen[1], maxLen[1], minLen[2], maxLen[2]); 

  // init seeds
  if ((Seeds = (VECTOR3 **)malloc(bf * sizeof(VECTOR3 *))) == NULL)
    Error("Error: Config() cannot allocate memory for Seeds\n");
  if ((NumSeeds = (int *)malloc(bf * sizeof(int))) == NULL)
    Error("Error: Config() cannot allocate memory for NumSeeds\n");
  if ((SizeSeeds = (int *)malloc(bf * sizeof(int))) == NULL)
    Error("Error: Config() cannot allocate memory for SizeSeeds\n");

  // init osuflow
  if ((osuflow = new OSUFlow*[bf]) == NULL)
    Error("Error: Config() cannot allocate memory for osuflow\n");

  // allocate streamline list for each block
  sl_list = new list<vtListSeedTrace*>[npart];

}
//-----------------------------------------------------------------------
//
// PrintSeeds
//
// prints the current seeds
// (debug)
//
void PrintSeeds(int nblocks) {

  int rank;
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < NumSeeds[i]; j++)
      fprintf(stderr,"Rank %d Block %d Seed %d: %.3f\t%.3f\t%.3f\n", rank, i, j, Seeds[i][j][0], Seeds[i][j][1], Seeds[i][j][2]);
  }

  fprintf(stderr,"\n");


}
//--------------------------------------------------------------------------
//
// RecvToSeeds
//
// copies received points to seeds
// local_block: block number within this process
// global_block: block number in entire domain
//
void RecvToSeeds(int local_block, int global_block) {

  int num;
  int i, j;

  num = lat->GetNumRecv(global_block);

  while (SizeSeeds[local_block] < num * sizeof(VECTOR3)) {

    Seeds[local_block] = (VECTOR3 *)realloc(Seeds[local_block], 
        SizeSeeds[local_block] * 2);
    if (Seeds[local_block] == NULL)
      Error("Error: RecvToSeeds() cannot reallocate memory\n");
    SizeSeeds[local_block] *= 2;

  }

  NumSeeds[local_block] = lat->CopyRecvToSeeds(global_block, 
       Seeds[local_block]);

}
//--------------------------------------------------------------------------
//
void draw_bounds(float xmin, float xmax, float ymin, float ymax, 
		 float zmin, float zmax) {

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
void draw_streamlines() {
  
  glPushMatrix(); 

  glScalef(1.0f / (float)size[0], 1.0f / (float)size[0], 1.0f / (float)size[0]); 
  glTranslatef(-size[0] / 2.0f, -size[1] / 2.0f, -size[2] / 2.0f); 

  std::list<vtListSeedTrace*>::iterator pIter; 

  for (int i = 0; i < npart; i++) {   // looping through all subdomains 

    pIter = sl_list[i].begin(); 

    for (; pIter != sl_list[i].end(); pIter++) {

      vtListSeedTrace *trace = *pIter; 
      std::list<VECTOR3*>::iterator pnIter; 
      pnIter = trace->begin(); 
      glColor3f(0.3,0.3,0.3); 
      glBegin(GL_LINE_STRIP); 
      for (; pnIter != trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
	//printf(" %f %f %f ", p[0], p[1], p[2]); 
	glVertex3f(p[0], p[1], p[2]); 
      }
      glEnd(); 

    }

    volume_bounds_type vb = vb_list[i];     
    draw_bounds(vb.xmin, vb.xmax, vb.ymin, vb.ymax, vb.zmin, vb.zmax); 

  }

  glPopMatrix(); 

}
//--------------------------------------------------------------------------
//
void animate_streamlines() {

}
//--------------------------------------------------------------------------
//
void draw_cube(float r, float g, float b) {

  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube 

}
//--------------------------------------------------------------------------
//
void display() {

  glEnable(GL_DEPTH_TEST); 
  glClearColor(1,1,1,1); 
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

//   ComputeStreamlines();

//   if (toggle_draw_streamlines == true)
    draw_streamlines(); 
//   else if (toggle_animate_streamlines == true)
//     animate_streamlines(); 

  glPushMatrix(); 
  glScalef(1.0, size[1] / size[0], size[2] / size[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 
  
  glBegin(GL_LINES); 
  glColor3f(1,0,0); 
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0);
  glColor3f(0,1,0);  
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1);  
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 


  glutSwapBuffers(); 

}
//--------------------------------------------------------------------------
//
void timer(int val) {

  if (toggle_animate_streamlines == true) {
    //    animate_streamlines(); 
    glutPostRedisplay(); 
  }
  glutTimerFunc(10, timer, 0); 

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
  else if (state == GLUT_UP) {
	  xform_mode = XFORM_NONE; 
  }
}
//--------------------------------------------------------------------------
//
void mymotion(int x, int y) {

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
//--------------------------------------------------------------------------
//
void mykey(unsigned char key, int x, int y) {

        switch(key) {
	case 'q': exit(1);
	  break; 
// 	case 's': compute_streamlines(); 
// 	  glutPostRedisplay(); 
// 	  break; 
// 	case 'd': 
// 	  toggle_draw_streamlines = !toggle_draw_streamlines; 
// 	  toggle_animate_streamlines = false; 
// 	  break; 
// 	case'a': 
// 	  toggle_animate_streamlines = !toggle_animate_streamlines; 
// 	  toggle_draw_streamlines = false; 
// 	  first_frame = 1; 
	}

}
//--------------------------------------------------------------------------

