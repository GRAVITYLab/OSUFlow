//------------------------------------------------------------------------------
//
// mpi test draw
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
//------------------------------------------------------------------------------

#include <mpi.h>
#include <omp.h>
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
bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 

// drawing data: only usable at the root process
VECTOR4 *pt; // points in everyone's traces, this round
int *npt; // everyone's number of points in their traces, this round
int tot_ntrace; // total number of everyone's traces, this round
VECTOR4 *Pt; // points in everyone's traces, all rounds
int *Npt; // everyone's number of points in their traces, all rounds
int Tot_ntrace; // total number of everyone's traces, all rounds

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void PrintSeeds(int nblocks);
int EndTrace(list<vtListSeedTrace*> &list, VECTOR3 &p, int index);
void ComputePathlines(int block_num);
void GatherPathlines();
int GatherNumPts(int* &ntrace);
void GatherPts(int *ntrace, int mynpt);
void DrawStreamlines();
void Cleanup();
void draw_bounds(float xmin, float xmax, float ymin, float ymax, 
		 float zmin, float zmax);
void draw_cube(float r, float g, float b);
void display();
void timer(int val);
void mymouse(int button, int state, int x, int y);
void mymotion(int x, int y);
void mykey(unsigned char key, int x, int y);
void idle();
void ComputeThread();
void IOThread();
void IOandCompute();
void ReceiveMessages();
int ComputeBlocksFit();
void MultiThreadEvictBlock(int round);
void SingleThreadEvictBlock();

// globals
static char filename[256]; // dataset file name
VECTOR3 size; // spatial domain size
static int tsize; // temporal domain size
int *NumSeeds; // number of seeds
int *SizeSeeds; // size of seeds list (bytes)
VECTOR4 **Seeds; // list of seeds lists
VECTOR3 *seeds; // one temporary list of (3d) seeds
OSUFlow **osuflow; // one flow object for each block
list<vtListTimeSeedTrace*> *sl_list; // pathlines list
int nspart; // global total number of spatial blocks
int ntpart; // global total number of temporal blocks
volume_bounds_type *vb_list; // global subdomain volume bounds list
int *blocks; // local list of blocks
int nblocks; // my number of blocks
Lattice4D* lat; // lattice
int tf; // max number of traces per block
int pf; // max number of points per trace
int max_rounds; // max number of rounds
int threads; // number of threads per process
int avail_mem; // memory space for dataset (MB)
int b_mem; // number of blocks to keep in memory
int max_bt; // max number of time steps in any block
int tr; // number of time partitions per round

// debug
#define MAX_RENDER_SEEDS 1000
VECTOR3 render_seeds[MAX_RENDER_SEEDS]; // seeds for rendering
int num_render_seeds = 0; // number of seeds for rendering
#define MAX_RENDER_PTS 20000
VECTOR3 render_pts[MAX_RENDER_PTS]; // seeds for rendering
int num_render_pts = 0; // number of seeds for rendering
//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  int nproc;  // mpi groupsize
  int rank; // mpi rank
  char buf[256];
  int mpi_thread_avail; // thread support in mpi implementation
  int i, j;

  // initialize mpi and the app
  GetArgs(argc, argv);
  if (threads == 1)
    MPI_Init(&argc, &argv);
  else {
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_thread_avail);
    if (mpi_thread_avail != MPI_THREAD_MULTIPLE)
      fprintf(stderr, "MPI thread support level = %d but it should be %d. Program will attempt to run, but results are undefined. Recommend running in single thread mode instead.\n",mpi_thread_avail, MPI_THREAD_MULTIPLE);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  Init();
  MPI_Barrier(MPI_COMM_WORLD);

  // multithread version

  if (threads == 2) {

    #pragma omp parallel num_threads(2)
    {
      // compute (master) thread
      if (!omp_get_thread_num())
	ComputeThread();

      // I/O (slave) thread
      if (omp_get_thread_num())
	IOThread();
    }

  }

  // single thread version

  if (threads == 1)
    IOandCompute();

#ifdef GRAPHICS

  // event loop for drawing
  if (rank == 0) {
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
  }

#endif

  Cleanup();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// ComputeThead
//
void ComputeThread() {

  int i, j;
  int done;
  int rank;
  int first; // first time
  int usec; // initial wait time (microseconds)
  int ngroups; // number of groups of blocks
  int g; // current group
  int sb, eb; // starting and ending block in the current group
  int bg; // number of blocks per group, except perhaps last group

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ngroups = ceil(ntpart / tr); // number of groups
  bg = floor(nblocks / ngroups); // number of blocks per group, except last

  // clear all blocks' compute status
  for (j = 0; j < nblocks; j++)
    lat->ClearComp(j);

  // for all groups
  for (g = 0; g < ngroups; g++) {

    if (rank == 0)
      fprintf(stderr, "** begin group %d **\n", g);

    // starting and ending blocks in the group
    sb = g * bg;
    eb = (g == ngroups - 1 ? nblocks : sb + bg);

    // debug
    fprintf(stderr, "sb = %d eb = %d\n",sb, eb);

    // for all rounds
    for (i = 0; i < max_rounds; i++) {

      if (rank == 0)
	fprintf(stderr, "begin round %d\n", i);

      // clear streamlines
      for (j = 0; j < nblocks; j++)
	sl_list[j].clear();

      done = 0;

      // until all blocks are computed
      usec = 1000;
      while (!done) {

	// for all blocks
	for (j = 0; j < nblocks; j++) {

	  // blocks that are not in this group do a null send to neighbors
	  if (i < sb || i >= eb) {
	    lat->SendNeighbors(blocks[i], MPI_COMM_WORLD);
#pragma omp flush
	    lat->SetComp(j, i);
#pragma omp flush
	  }

	  // blocks that are in this group get computed
	  if (!lat->GetComp(j, i) && lat->GetLoad(j)) {
	    ComputePathlines(j);
	    fprintf(stderr, "Computed block %d\n", j);
#pragma omp flush
	    lat->SetComp(j, i);
#pragma omp flush
	  }

	} // for all blocks

	// check if all done
	done = 1;
	for (j = 0; j < nblocks; j++) {
	  if (!lat->GetComp(j, i))
	    done = 0;
	}

	if (!done) {
	  usleep(usec);
	  usec *= 2;
	}

      } // until all blocks are computed

      ReceiveMessages();

    } // for all rounds

    if (rank == 0)
      fprintf(stderr, "Completed %d rounds\n", max_rounds);

  } // for all groups

  if (rank == 0)
    fprintf(stderr, "Completed %d groups\n", ngroups);

#ifdef GRAPHICS

  GatherPathlines();

#endif

}
//-----------------------------------------------------------------------
//
// IOThread
//
void IOThread() {

  VECTOR3 minB, maxB; // subdomain bounds
  float from[3], to[3]; // seed points limits
  int num_loaded; // number of blocks loaded into memory so far
  int rank;
  int first = 1; // first time
  int usec; // wait time in microseconds
  int i, j, k;
  int ngroups; // number of groups of blocks
  int g; // current group
  int sb, eb; // starting and ending block in the current group
  int bg; // number of blocks per group, except perhaps last group

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  num_loaded = 0;
  ngroups = ceil(ntpart / tr); // number of groups
  bg = floor(nblocks / ngroups); // number of blocks per group, except last

  // init all blocks
  for (i = 0; i < nblocks; i++) {

    lat->ClearLoad(i);
    osuflow[i] = new OSUFlow;
    minB.Set(vb_list[blocks[i]].xmin, vb_list[blocks[i]].ymin, 
	     vb_list[blocks[i]].zmin);
    maxB.Set(vb_list[blocks[i]].xmax, vb_list[blocks[i]].ymax, 
	     vb_list[blocks[i]].zmax);

    // init seeds for blocks at t = initial time
//     if (vb_list[blocks[i]].tmin == 0) {

      from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
      to[0]   = maxB[0];   to[1]   = maxB[1];   to[2]   = maxB[2]; 
      osuflow[i]->SetRandomSeedPoints(from, to, tf); 
      seeds = osuflow[i]->GetSeeds(NumSeeds[i]); 

      while (SizeSeeds[i] < NumSeeds[i] * sizeof(VECTOR4)) {
	Seeds[i] = (VECTOR4 *)realloc(Seeds[i], SizeSeeds[i] * 2);
	assert(Seeds[i] != NULL);
	SizeSeeds[i] *= 2;
      }

      for (k = 0; k < NumSeeds[i]; k++)
	Seeds[i][k].Set(seeds[k][0], seeds[k][1], seeds[k][2], 
			vb_list[blocks[i]].tmin);

//     } // init seeds

  } // init all blocks

  // for all groups
  for (g = 0; g < ngroups; g++) {

    if (rank == 0)
      fprintf(stderr, "** begin group %d **\n", g);

    // starting and ending blocks in the group
    sb = g * bg;
    eb = (g == ngroups - 1 ? nblocks : sb + bg);

    // debug
    fprintf(stderr, "sb = %d eb = %d\n",sb, eb);

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

      // for all blocks
      for (i = 0; i < nblocks; i++) {

	// block is either not in this group or is loaded already
	if (i < sb || i >= eb || lat->GetLoad(i))
	  continue;

	// make room for the next block
	if (num_loaded >= b_mem) {
	  MultiThreadEvictBlock(j);
	  num_loaded--;
	}

	// read the data
	fprintf(stderr, "block %d: min %.3lf %.3lf %.3lf %d max %.3lf %.3lf %.3lf %d\n", i, minB[0], minB[1], minB[2], vb_list[blocks[i]].tmin, maxB[0], maxB[1], maxB[2], vb_list[blocks[i]].tmax);
	osuflow[i]->LoadData(filename, false, minB, maxB, size, max_bt,
			     vb_list[blocks[i]].tmin, vb_list[blocks[i]].tmax, 
			     MPI_COMM_WORLD); 

	// update status
	lat->SetData(blocks[i], 1);
	fprintf(stderr, "Loaded block %d\n", i);
#pragma omp flush
	lat->SetLoad(i);
#pragma omp flush
	num_loaded++;
	osuflow[i]->ScaleField(10.0); // improves visibility

      } // for all blocks

      usleep(1000);

    } // for all rounds

  } // for all groups

}
//-----------------------------------------------------------------------
//
// MultiThreadEvictBlock
//
// evicts next block from memory (in block number order)
// multithread version
//
// round: round number
//
void MultiThreadEvictBlock(int round) {

  static int target = 0; // block to evict
  int first = 1; // first time
  int usec = 1000; // initial wait time (microseconds)

  // wait for the block to be computed before evicting
  while (!lat->GetComp(target, round)) {
    if (first)
      fprintf(stderr, "Waiting to evict block %d...\n", target);
    first = 0;
    usleep(usec);
    usec *= 2;
  }

  // evict it
  #pragma omp flush
  lat->ClearLoad(target);
  #pragma omp flush
  osuflow[target]->DeleteData();
  lat->SetData(blocks[target], 0);

  fprintf(stderr, "Evicted block %d\n", target);

  target = (target + 1) % nblocks;

}
//-----------------------------------------------------------------------
//
// IOandCompute
//
void IOandCompute() {

  VECTOR3 minB, maxB; // subdomain bounds
  float from[3], to[3]; // seed points limits
  int num_loaded; // number of blocks loaded into memory so far
  int rank;
  int i, j, k;
  int ngroups; // number of groups of blocks
  int g; // current group
  int sb, eb; // starting and ending block in the current group
  int bg; // number of blocks per group, except perhaps last group

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  num_loaded = 0; // number of blocks currently loaded
  ngroups = ceil(ntpart / tr); // number of groups
  bg = floor(nblocks / ngroups); // number of blocks per group, except last

  // init all blocks
  for (i = 0; i < nblocks; i++) {

    lat->ClearLoad(i);
    osuflow[i] = new OSUFlow;
    minB.Set(vb_list[blocks[i]].xmin, vb_list[blocks[i]].ymin, 
	     vb_list[blocks[i]].zmin);
    maxB.Set(vb_list[blocks[i]].xmax, vb_list[blocks[i]].ymax, 
	     vb_list[blocks[i]].zmax);

    // init seeds for blocks at t = initial time
//     if (vb_list[blocks[i]].tmin == 0) {

      from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
      to[0]   = maxB[0];   to[1]   = maxB[1];   to[2]   = maxB[2]; 
      osuflow[i]->SetRandomSeedPoints(from, to, tf); 
      seeds = osuflow[i]->GetSeeds(NumSeeds[i]); 

      while (SizeSeeds[i] < NumSeeds[i] * sizeof(VECTOR4)) {
	Seeds[i] = (VECTOR4 *)realloc(Seeds[i], SizeSeeds[i] * 2);
	assert(Seeds[i] != NULL);
	SizeSeeds[i] *= 2;
      }

      for (k = 0; k < NumSeeds[i]; k++)
	Seeds[i][k].Set(seeds[k][0], seeds[k][1], seeds[k][2], 
			vb_list[blocks[i]].tmin);

//     } // init seeds

  } // init all blocks

  // for all groups
  for (g = 0; g < ngroups; g++) {

    if (rank == 0)
      fprintf(stderr, "** begin group %d **\n", g);

    // starting and ending blocks in the group
    sb = g * bg;
    eb = (g == ngroups - 1 ? nblocks : sb + bg);

    // debug
    fprintf(stderr, "sb = %d eb = %d\n",sb, eb);

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

      if (rank == 0)
	fprintf(stderr, " * begin round %d *\n", j);

      // for all blocks
      for (i = 0; i < nblocks; i++) {

	// blocks that are not in this group do a null send to neighbors
	if (i < sb || i >= eb) {
	  lat->SendNeighbors(blocks[i], MPI_COMM_WORLD);
	  continue;
	}

	// if the block needs to be loaded
	if (!lat->GetLoad(i)) {

	  // make room for the next block
	  if (num_loaded >= b_mem) {
	    SingleThreadEvictBlock();
	    num_loaded--;
	  }

	  // read the data
	  fprintf(stderr, "rank %d: min %.3lf %.3lf %.3lf %d max %.3lf %.3lf %.3lf %d\n", blocks[i], minB[0], minB[1], minB[2], vb_list[blocks[i]].tmin, maxB[0], maxB[1], maxB[2], vb_list[blocks[i]].tmax);
	  osuflow[i]->LoadData(filename, false, minB, maxB, size, max_bt,
			       vb_list[blocks[i]].tmin, 
			       vb_list[blocks[i]].tmax, MPI_COMM_WORLD); 
	  lat->SetData(blocks[i], 1);
	  lat->SetLoad(i);
	  num_loaded++;
	  fprintf(stderr, "Loaded block %d\n", i);
	  osuflow[i]->ScaleField(10.0); // improves visibility

	} // if the block needs to be loaded

	// compute pathlines
	ComputePathlines(i);
	fprintf(stderr, "Computed block %d\n", i);

      } // for all blocks

      ReceiveMessages();

    } // for all rounds

    if (rank == 0)
      fprintf(stderr, "Completed %d rounds\n", max_rounds);

  } // for all groups

  if (rank == 0)
    fprintf(stderr, "Completed %d groups\n", ngroups);

  // gather pathlines at root for rendering
#ifdef GRAPHICS

  GatherPathlines();

#endif

}
//-----------------------------------------------------------------------
//
// SingleTheadEvictBlock
//
// evicts next block from memory (in block number order)
//
// single thread version
//
void SingleThreadEvictBlock() {

  static int target = 0; // block to evict

  lat->ClearLoad(target);
  osuflow[target]->DeleteData();
  lat->SetData(blocks[target], 0);
  fprintf(stderr, "Evicted block %d\n", target);
  target = (target + 1) % nblocks;
  
}
//-----------------------------------------------------------------------
//
// ComputePathlines
//
// computes path or streamlines
// 
// block_num: local block number (0 to nblocks-1)
// not global partition number
//
//
void ComputePathlines(int block_num) {

  list<vtListTimeSeedTrace*> list; // trace of seed points
  std::list<VECTOR4*>::iterator pt_iter; // iterator over pts in one trace
  std::list<vtListTimeSeedTrace*>::iterator trace_iter; // iter. over traces
  VECTOR4 p; // current point
  int neighbor; // neighbor's number (0-79)
  int i, j;

  if (NumSeeds[block_num]) {

    list.clear();

    // perform the integration
    // todo: integrate in both directions
    osuflow[block_num]->SetIntegrationParams(1, 5); 
    osuflow[block_num]->GenPathLines(Seeds[block_num], list, FORWARD, 
         NumSeeds[block_num], pf); 

    // debug
    fprintf(stderr, "NumSeeds[rank %d] = %d NumTraces[rank %d] = %d\n", 
	    blocks[block_num], NumSeeds[block_num], blocks[block_num], 
	    list.size());

    // copy each trace to the streamline list for later rendering
    for (trace_iter = list.begin(); trace_iter != list.end(); trace_iter++)
      sl_list[block_num].push_back(*trace_iter); 

    // redistribute end points of traces to neighbors

    // for all end points of traces
    for (trace_iter = list.begin(); trace_iter != list.end(); trace_iter++) {
      if (!(*trace_iter)->size())
	continue;
      pt_iter = (*trace_iter)->end();
      pt_iter--;
      p = **pt_iter;

      // find which neighbor the point is in
      neighbor = lat->GetNeighbor(blocks[block_num], p[0], p[1], p[2], p[3]);

      // post the point to the send list
      if (neighbor != -1)
	lat->PostPoint(blocks[block_num], p, neighbor);

    } // for all end points

    // send points to neighbors
    lat->SendNeighbors(blocks[block_num], MPI_COMM_WORLD);

  } // if (NumSeeds[block_num])
      
}
//-----------------------------------------------------------------------
//
// ReceiveMessages
//
// receives all pending messages
//
//
void ReceiveMessages() {

  int i;
  int npts;

  for (i = 0; i < nblocks; i++) {

    if ((npts = lat->ReceiveNeighbors(blocks[i], MPI_COMM_WORLD))) {

      NumSeeds[i] = npts;

      while (SizeSeeds[i] < NumSeeds[i] * sizeof(VECTOR4)) {
	Seeds[i] = (VECTOR4 *)realloc(Seeds[i], SizeSeeds[i] * 2);
	assert(Seeds[i] != NULL);
	SizeSeeds[i] *= 2;
      }

      lat->GetRecvPts(blocks[i], Seeds[i]);

    }

  }

}
//-----------------------------------------------------------------------
//
// GatherPathlines
//
// gathers all pathlines at the root for rendering
//
void GatherPathlines() {

  static int *ntrace = NULL; // number of traces for each proc
  int n; // total number of my points
  int i;

  // gather number of points in each trace at the root
  n = GatherNumPts(ntrace);
  
  // gather the actual points in each trace at the root
  GatherPts(ntrace, n);

}
//-----------------------------------------------------------------------
//
// GatherNumPts
//
// gathers number of points in each trace to the root
//
// ntrace: number of traces in each process (passed by reference)
// sb, eb: starting and ending blocks
//
// returns: total number of points in my process
//
int GatherNumPts(int* &ntrace) {

  int myntrace = 0; // my number of traces
  static int *ofst = NULL; // offsets into ntrace
  int *mynpt; // number of points in each of my traces
  int tot_mynpt = 0; // total number of my points
  int rank, nproc; // MPI usual
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (ntrace == NULL && (ntrace = new int[nproc]) == NULL)
    Error("Error: GatherNumPts cannot allocate ntrace\n");

  if (ofst == NULL && (ofst = new int[nproc]) == NULL)
    Error("Error: GatherNumPts cannot allocate ofst\n");

  // compute number of my traces
  for (i = 0; i < nblocks; i++)
    myntrace += sl_list[i].size();

  // gather number of traces at the root
  MPI_Gather(&myntrace, 1, MPI_INT, ntrace, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // compute number of points in each of my traces
  if ((mynpt = new int[myntrace]) == NULL)
    Error("Error: GatherNumPts cannot allocate mynpt\n");

  j = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      if (j >= myntrace) {
	fprintf(stderr,"Warning: GatherNumPts() j should be < myntrace but j = %d myntrace = %d. This should not happen\n",j,myntrace);
	break;
      }
      mynpt[j] = (*trace_iter)->size();
      tot_mynpt += mynpt[j++];
    }
  }

  // gather number of points in each trace at the root
  if (rank == 0) {

    tot_ntrace = 0;
    for (i = 0; i < nproc; i++) {
      ofst[i] = (i == 0) ? 0 : ofst[i - 1] + ntrace[i - 1];
      tot_ntrace += ntrace[i];
    }

    if (ofst[nproc - 1] + ntrace[nproc - 1] > nproc * nblocks * tf)
      Error("Error: GatherNumPts attempted to write beyond the end of npt\n");

  }

  MPI_Gatherv(mynpt, myntrace, MPI_INT, npt, ntrace, ofst, MPI_INT, 0, 
	      MPI_COMM_WORLD);

  delete[] mynpt;

  return tot_mynpt;

}
//-----------------------------------------------------------------------
//
// GatherPts
//
// gathers the points in each trace at the root
//
// ntrace: number of traces in each process
// mynpt: total number of points in my process
//
void GatherPts(int *ntrace, int mynpt) {

  static int *nflt = NULL; // number of floats in points from each proc
  static int *ofst = NULL; // offsets into pt
  VECTOR4 *mypt; // points in my traces
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  int rank, nproc; // MPI usual
  int i, j, k;
  static int next_pt = 0; // current next open slot in Pt (all points)
  static int next_trace = 0; // current next open slot in Npt (all traces)

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (nflt == NULL && (nflt = new int[nproc]) == NULL)
    Error("Error: GatherPts cannot allocate nflt\n");

  if (ofst == NULL && (ofst = new int[nproc]) == NULL)
    Error("Error: GatherPts cannot allocate ofst\n");

  if ((mypt = new VECTOR4[mynpt]) == NULL)
    Error("Error: GatherPts cannot allocate mypt\n");

  // collect my own points
  j = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
         pt_iter++)
	mypt[j++] = **pt_iter;
    }
  }

  // gather the points at the root
  if (rank == 0) {

    k = 0;
    for (i = 0; i < nproc; i++) {
      nflt[i] = 0;
      for (j = 0; j < ntrace[i]; j++)
	nflt[i] += (npt[k++] * 4);
      ofst[i] = (i == 0) ? 0 : ofst[i - 1] + nflt[i - 1];
    }

    if (ofst[nproc - 1] + nflt[nproc - 1] > 4 * nproc * nblocks * tf * pf)
      Error("Error: GatherPts attempted to write beyond the end of pt\n");

  }

  MPI_Gatherv(mypt, mynpt * 4, MPI_FLOAT, pt, nflt, ofst,
	      MPI_FLOAT, 0, MPI_COMM_WORLD);

  // concatenate points from current round with previous rounds
  k = 0;
  for (int i = 0; i < tot_ntrace; i++) {

    Npt[next_trace] = npt[i];
    next_trace++;
    for (j = 0; j < npt[i]; j++) {
      Pt[next_pt][0] = pt[k][0];
      Pt[next_pt][1] = pt[k][1];
      Pt[next_pt][2] = pt[k][2];
      Pt[next_pt][3] = pt[k][3];
      next_pt++;
      k++;
    }

  }
  Tot_ntrace += tot_ntrace;

  delete[] mypt;

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

  assert(argc >= 14);

  strncpy(filename,argv[1],sizeof(filename));

  // hard code the minimum corner to be 0,0,0,0
  // need to allow for variable data origin in the future
  minLen[0] = minLen[1] = minLen[2] = minTime = 0.0;

  maxLen[0] = atof(argv[2]) - 1.0f;
  maxLen[1] = atof(argv[3]) - 1.0f;
  maxLen[2] = atof(argv[4]) - 1.0f;
  maxTime   = atof(argv[5]) - 1.0f;

  size[0] = maxLen[0] - minLen[0] + 1.0f; // data sizes, space and time
  size[1] = maxLen[1] - minLen[1] + 1.0f;
  size[2] = maxLen[2] - minLen[2] + 1.0f;
  tsize   = maxTime - minTime + 1.0f;

  nspart = atoi(argv[6]); // total space partitions
  ntpart = atoi(argv[7]); // total time partitions
  tr = atoi(argv[8]); // number of time partitions per round

  tf = atoi(argv[9]); // traces per block
  pf = atoi(argv[10]); // points per trace
  max_rounds = atoi(argv[11]); // rounds
  threads = atoi(argv[12]) <= 1 ? 1 : 2; // threads per process
  avail_mem = atoi(argv[13]); // memory data size (MB)

}
//-----------------------------------------------------------------------
//
// Init
//
// inits the app
//
void Init() {

  int rank, nproc; // usual MPI
  int nt; // max total number of traces
  int np; // max total number of points
  int ghost = 1; // number of ghost cells per spatial edge
  int b_size; // data size in a typical block (Mbytes)
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  assert(nspart * ntpart >= nproc);
  assert(tr <= ntpart);

  // init lattice and osuflow
  lat = new Lattice4D(size[0], size[1], size[2], tsize, ghost, nspart, 
		      ntpart, 4);
  vb_list = lat->GetBoundsList(); 
  lat->InitSeedLists(); 
  lat->RoundRobin_proc(nproc); 
  nblocks = lat->GetNumPartitions(rank);
  blocks = (int *)malloc(nblocks * sizeof(int));
  assert(blocks != NULL);
  lat->GetPartitions(rank, blocks);
  osuflow = new OSUFlow*[nblocks];
  assert(osuflow != NULL);
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

  // allocate pts list and number of points list for rendering
  if (rank == 0) {
    nt = nspart * ntpart * tf; // max total traces, current round
    np = nt * pf; // max total points, current round
    npt = new int[nt]; // number of points in everyone's trace, current round
    assert(npt != NULL);
    pt = new VECTOR4[np]; // points in everyones traces, current round
    assert(pt != NULL);
    // number of points in everyone's traces, all rounds
    Npt = new int[nt * max_rounds];
    assert(Npt != NULL);
    // points ini everyone's traces, all rounds
    Pt = new VECTOR4[np * max_rounds];
    assert(Pt != NULL);
  }

  // max number of time steps in any block
  max_bt = ceil(tsize / ntpart) + 2 * ghost;

  // number of blocks to keep in memory
  b_size = size[0] * size[1] / 1048576.0f * size[2] / nspart * 
           tsize / ntpart * 3 * sizeof (float);
  b_mem = avail_mem / b_size;

  // print some of the args
  if (rank == 0) {
    fprintf(stderr,"Volume size: X %.3lf Y %.3lf Z %.3lf t %d\n",
	    size[0], size[1], size[2], tsize);
    fprintf(stderr, "Number of threads per process: %d\n", threads);
    fprintf(stderr, "Number of compute rounds: %d\n", max_rounds);
    fprintf(stderr, "Available dataset memory per process: %d MB\n", avail_mem);
    fprintf(stderr, "Number of blocks a process can fit in memory: %d\n", b_mem);
  }

}
//-----------------------------------------------------------------------
//
// Cleanup
//
// frees memory and such
//
void Cleanup() {

  int i;

  delete [] Pt;
  delete [] Npt;
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
  free(blocks);
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

  int rank;
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 0; i < nblocks; i++) {
    for (j = 0; j < NumSeeds[i]; j++) {
      fprintf(stderr,"Rank %d Block %d Seed %d: %.3f\t%.3f\t%.3f\t%.3f\n", rank, i, j, Seeds[i][j][0], Seeds[i][j][1], Seeds[i][j][2], Seeds[i][j][3]);
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
void DrawPathlines() {
  
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
  for (int i = 0; i < Tot_ntrace; i++) {

    glBegin(GL_LINE_STRIP); 
    for (j = 0; j < Npt[i]; j++) {
      if (Pt[k][3] <= step)
	glVertex3f(Pt[k][0], Pt[k][1], Pt[k][2]);
      k++;
    }
    glEnd(); 

  }

  for (i = 0; i < nspart * ntpart; i++) {
    volume_bounds_type vb = vb_list[i];     
    draw_bounds(vb.xmin, vb.xmax, vb.ymin, vb.ymax, vb.zmin, vb.zmax); 
  }

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

  DrawPathlines(); 

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
	case 'q': exit(1);
	  break; 
	}

}
//--------------------------------------------------------------------------

#endif
