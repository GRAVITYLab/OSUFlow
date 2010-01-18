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
//--------------------------------------------------------------------------

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

// drawing data: only usable at the root process
VECTOR4 *pt; // points in everyone's traces
int *npt; // everyone's number of points in their traces
int tot_ntrace; // total number of everyone's traces

// performance stats
int *block_stats; // block stats
double *time_stats; // time stats
int n_block_stats = 5; // number of block stats
int n_time_stats = 4; // number of time stats
int TotSeeds = 0; // total number of seeds for all blocks and all rounds
                  // in this process
int TotRounds = 0; // total number of rounds this process executed excluding
                   // idle rounds at the end
double TotIOBW = 0.0; // total IO bandwidth in MB/s
double TotTime = 0.0; // total time
double TotIOTime = 0.0; // total IO time
double TotCompTime = 0.0; // total computation time

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void PrintSeeds(int nblocks);
int EndTrace(list<vtListSeedTrace*> &list, VECTOR3 &p, int index);
void ComputeFieldlines(int block_num);
void ComputePathlines(int block_num);
void ComputeStreamlines(int block_num);
void GatherFieldlines();
int GatherNumPts(int* &ntrace, int all);
void GatherPts(int *ntrace, int mynpt);
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
void ComputeThread();
void IOThread();
void IOandCompute();
int ComputeBlocksFit();
void MultiThreadEvictBlock(int round);
void SingleThreadEvictBlock();
void PrintPerf();
void WriteFieldlines(int *ntrace, int mynpt, char *filename);
void Partition(char *filename);

// globals
static char filename[256]; // dataset file name
char part_file[256]; // partition file name
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
int threads; // number of threads per process
int avail_mem; // memory space for dataset (MB)
int b_mem; // number of blocks to keep in memory
int max_bt; // max number of time steps in any block
int tr; // number of time partitions per round
int data_mode; // data format

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

  // start the total time
  TotTime = MPI_Wtime();

  Init();
  MPI_Barrier(MPI_COMM_WORLD);

#ifndef MAC_OSX_10_4

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

#endif

  // single thread version

  if (threads == 1)
    IOandCompute();

  // stop the total time and print the performance stats
  MPI_Barrier(MPI_COMM_WORLD);
  TotTime = MPI_Wtime() - TotTime;
  PrintPerf();

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
// PrintPerf
//
void PrintPerf() {

  int nproc;  // mpi groupsize
  int rank; // mpi rank
  int *all_block_stats; // gathered block stats
  double *all_time_stats; // gathered time stats
  float TotCells; // total number of spatial data cells in billions
  float TotDataSize; // total size of a time step in GB
  int i;

  int tot_npart = 0; // number of partitions per proc
  int min_npart;
  int max_npart;
  int mean_npart;
  float var_npart = 0.0;
  float std_npart;

  int tot_nneigh = 0; // number of neighbors 
  int min_nneigh;
  int max_nneigh;
  int mean_nneigh;
  float var_nneigh = 0.0;
  float std_nneigh;

  int tot_nseed = 0; // number of seeds
  int min_nseed;
  int max_nseed;
  int mean_nseed;
  float var_nseed = 0.0;
  float std_nseed;

  int tot_nround = 0; // number of rounds before any null rounds at the end
  int min_nround;
  int max_nround;
  int mean_nround;
  float var_nround = 0.0;
  float std_nround;

  int tot_nptsend = 0; // total number of points sent
  int min_nptsend;
  int max_nptsend;
  int mean_nptsend;
  float var_nptsend = 0.0;
  float std_nptsend;

  double tot_iotime = 0.0; // I/O time
  double min_iotime;
  double max_iotime;
  double mean_iotime;
  double var_iotime = 0.0;
  double std_iotime;

  double tot_commtime = 0.0; // communication time
  double min_commtime;
  double max_commtime;
  double mean_commtime;
  double var_commtime = 0.0;
  double std_commtime;

  double tot_comptime = 0.0; // computation time
  double min_comptime;
  double max_comptime;
  double mean_comptime;
  double var_comptime = 0.0;
  double std_comptime;

  double tot_iobw = 0; // I/O bandwidth in MB/s
  double min_iobw;
  double max_iobw;
  double mean_iobw;
  double var_iobw = 0.0;
  double std_iobw;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // get stats of my process
  block_stats[0] = lat->GetMyNumPartitions(); // number of blocks
  block_stats[1] = lat->GetMyNumNeighbors(); // avg number of neighbors
  block_stats[2] = TotSeeds; // total number of particles advected
  block_stats[3] = TotRounds; // total number of rounds before any null rounds
                              // at the end
  block_stats[4] = lat->GetMyTotPtsSend(); // total points sent
  time_stats[0] = TotIOTime; // I/O time
  time_stats[1] = lat->GetMyCommTime(); // communication time
  time_stats[2] = TotCompTime; // computation time
  time_stats[3] = TotIOBW; // I/O bandwidth

  // alloc space and gather the stats
  assert((all_block_stats = (int *)malloc(n_block_stats * nproc * 
					  sizeof(int))) != NULL);
  assert((all_time_stats = (double *)malloc(n_time_stats * nproc * 
					    sizeof(double))) != NULL);
  MPI_Gather(block_stats, n_block_stats, MPI_INT, all_block_stats, 
	     n_block_stats, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(time_stats, n_time_stats, MPI_DOUBLE, all_time_stats, 
	     n_time_stats, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // print the stats
  if (rank == 0) {

    // totals, mins, maxs
    for (i = 0; i < nproc; i++) {

      // total
      tot_npart    += all_block_stats[n_block_stats * i];
      tot_nneigh   += all_block_stats[n_block_stats * i + 1];
      tot_nseed    += all_block_stats[n_block_stats * i + 2];
      tot_nround   += all_block_stats[n_block_stats * i + 3];
      tot_nptsend += all_block_stats[n_block_stats * i + 4];
      tot_iotime   += all_time_stats[n_time_stats * i];
      tot_commtime += all_time_stats[n_time_stats * i + 1];
      tot_comptime += all_time_stats[n_time_stats * i + 2];
      tot_iobw     += all_time_stats[n_time_stats * i + 3];

      // min, max
      if (i == 0) {
	min_npart = max_npart       = all_block_stats[n_block_stats * i];
	min_nneigh = max_nneigh     = all_block_stats[n_block_stats * i + 1];
	min_nseed = max_nseed       = all_block_stats[n_block_stats * i + 2];
	min_nround = max_nround     = all_block_stats[n_block_stats * i + 3];
	min_nptsend = max_nptsend   = all_block_stats[n_block_stats * i + 4];
	min_iotime = max_iotime     = all_time_stats[n_time_stats * i];
	min_commtime = max_commtime = all_time_stats[n_time_stats * i + 1];
	min_comptime = max_comptime = all_time_stats[n_time_stats * i + 2];
	min_iobw = max_iobw         = all_time_stats[n_time_stats * i + 3];
      }
      else {

	if (all_block_stats[n_block_stats * i] < min_npart)
	  min_npart = all_block_stats[n_block_stats * i];
	if (all_block_stats[n_block_stats * i] > max_npart)
	  max_npart = all_block_stats[n_block_stats * i];

	if (all_block_stats[n_block_stats * i + 1] < min_nneigh)
	  min_nneigh = all_block_stats[n_block_stats * i + 1];
	if (all_block_stats[n_block_stats * i + 1] > max_nneigh)
	  max_nneigh = all_block_stats[n_block_stats * i + 1];

	if (all_block_stats[n_block_stats * i + 2] < min_nseed)
	  min_nseed = all_block_stats[n_block_stats * i + 2];
	if (all_block_stats[n_block_stats * i + 2] > max_nseed)
	  max_nseed = all_block_stats[n_block_stats * i + 2];

	if (all_block_stats[n_block_stats * i + 3] < min_nround)
	  min_nround = all_block_stats[n_block_stats * i + 3];
	if (all_block_stats[n_block_stats * i + 3] > max_nround)
	  max_nround = all_block_stats[n_block_stats * i + 3];

	if (all_block_stats[n_block_stats * i + 4] < min_nptsend)
	  min_nptsend = all_block_stats[n_block_stats * i + 4];
	if (all_block_stats[n_block_stats * i + 4] > max_nptsend)
	  max_nptsend = all_block_stats[n_block_stats * i + 4];

	if (all_time_stats[n_time_stats * i] < min_iotime)
	  min_iotime = all_time_stats[n_time_stats * i];
	if (all_time_stats[n_time_stats * i] > max_iotime)
	  max_iotime = all_time_stats[n_time_stats * i];

	if (all_time_stats[n_time_stats * i + 1] < min_commtime)
	  min_commtime = all_time_stats[n_time_stats * i + 1];
	if (all_time_stats[n_time_stats * i + 1] > max_commtime)
	  max_commtime = all_time_stats[n_time_stats * i + 1];

	if (all_time_stats[n_time_stats * i + 2] < min_comptime)
	  min_comptime = all_time_stats[n_time_stats * i + 2];
	if (all_time_stats[n_time_stats * i + 2] > max_comptime)
	  max_comptime = all_time_stats[n_time_stats * i + 2];

	if (all_time_stats[n_time_stats * i + 3] < min_iobw)
	  min_iobw = all_time_stats[n_time_stats * i + 3];
	if (all_time_stats[n_time_stats * i + 3] > max_iobw)
	  max_iobw = all_time_stats[n_time_stats * i + 3];

      }
    }

    // means
    mean_npart = tot_npart / nproc;
    mean_nneigh = tot_nneigh / nproc;
    mean_nseed = tot_nseed / nproc;
    mean_nround = tot_nround / nproc;
    mean_nptsend = tot_nptsend / nproc;
    mean_iotime = tot_iotime / nproc;
    mean_commtime = tot_commtime / nproc;
    mean_comptime = tot_comptime / nproc;
    mean_iobw = tot_iobw / nproc;

    // variances
    for (i = 0; i < nproc; i++) {
      var_npart += (all_block_stats[n_block_stats * i] - mean_npart) *
	(all_block_stats[n_block_stats * i] - mean_npart);
      var_nneigh += (all_block_stats[n_block_stats * i + 1] - mean_nneigh) *
	(all_block_stats[n_block_stats * i + 1] - mean_nneigh);
      var_nseed += (all_block_stats[n_block_stats * i + 2] - mean_nseed) *
	(all_block_stats[n_block_stats * i + 2] - mean_nseed);
      var_nround += (all_block_stats[n_block_stats * i + 3] - mean_nround) *
	(all_block_stats[n_block_stats * i + 3] - mean_nround);
      var_nptsend += (all_block_stats[n_block_stats * i + 4] - mean_nptsend) *
	(all_block_stats[n_block_stats * i + 4] - mean_nptsend);
      var_iotime += (all_time_stats[n_time_stats * i] - mean_iotime) *
	(all_time_stats[n_time_stats * i] - mean_iotime);
      var_commtime += (all_time_stats[n_time_stats * i + 1] - mean_commtime) *
	(all_time_stats[n_time_stats * i + 1] - mean_commtime);
      var_comptime += (all_time_stats[n_time_stats * i + 2] - mean_comptime) *
	(all_time_stats[n_time_stats * i + 2] - mean_comptime);
      var_iobw += (all_time_stats[n_time_stats * i + 3] - mean_iobw) *
	(all_time_stats[n_time_stats * i + 3] - mean_iobw);
    }
    var_npart /= nproc;
    var_nneigh /= nproc;
    var_nseed /= nproc;
    var_nround /= nproc;
    var_nptsend /= nproc;
    var_iotime /= nproc;
    var_commtime /= nproc;
    var_comptime /= nproc;
    var_iobw /= nproc;

    // standard deviations
    std_npart = sqrt(var_npart);
    std_nneigh = sqrt(var_nneigh);
    std_nseed = sqrt(var_nseed);
    std_nround = sqrt(var_nround);
    std_nptsend = sqrt(var_nptsend);
    std_iotime = sqrt(var_iotime);
    std_commtime = sqrt(var_commtime);
    std_comptime = sqrt(var_comptime);
    std_iobw = sqrt(var_iobw);

    // misc: data size and aggregate bandwidth
    TotCells = size[0] * size[1] * size[2] / 1.0e6;
    TotDataSize = size[0] * size[1] * size[2] * 12 / 1048576;

    // print results
    fprintf(stderr, "----- Performance Summary -----\n");
    fprintf(stderr, "Number of procs = %d\n", nproc);
    fprintf(stderr, "Total time = %.2lf s\n", TotTime);
    fprintf(stderr, "Total data size = %.2f million cells = %.2f MB\n", TotCells, TotDataSize);
    fprintf(stderr, "Total particles = %.2f million\n", tot_nseed / 1.0e6);
    fprintf(stderr, "Aggregate I/O bandwidth = %.0lf MB/s\n", TotDataSize / mean_iotime);
    fprintf(stderr, "Blocks / proc %14s min = %-8d max = %-8d avg = %-8d var = %-8.0f std = %-8.0f\n", "", min_npart, max_npart, mean_npart, 
	    var_npart, std_npart);
    fprintf(stderr, "Neighbors / block %10s min = %-8d max = %-8d avg = %-8d var = %-8.0f std = %-8.0f\n", "", min_nneigh, max_nneigh, mean_nneigh, 
	    var_nneigh, std_nneigh);
    fprintf(stderr, "Total pts comp / proc %6s min = %-8d max = %-8d avg = %-8d var = %-8.0f std = %-8.0f\n", "", min_nseed, max_nseed, mean_nseed, 
	    var_nseed, std_nseed);
    fprintf(stderr, "Rounds / proc %14s min = %-8d max = %-8d avg = %-8d var = %-8.0f std = %-8.0f\n", "", min_nround, max_nround, mean_nround, 
	    var_nround, std_nround);
    fprintf(stderr, "Total pts sent / proc %6s min = %-8d max = %-8d avg = %-8d var = %-8.0f std = %-8.0f\n", "", min_nptsend, max_nptsend, mean_nptsend, 
	    var_nptsend, std_nptsend);
    fprintf(stderr, "I/O time / proc (s) %8s min = %-8.2lf max = %-8.2lf avg = %-8.2lf var = %-8.2lf std = %-8.2lf\n", "", min_iotime, max_iotime, mean_iotime, 
	    var_iotime, std_iotime);
    fprintf(stderr, "Comp time / proc (s) %7s min = %-8.2lf max = %-8.2lf avg = %-8.2lf var = %-8.2lf std = %-8.2lf\n", "", min_comptime, max_comptime, mean_comptime, 
	    var_commtime, std_commtime);
    fprintf(stderr, "Comm time / proc (s) %7s min = %-8.2lf max = %-8.2lf avg = %-8.2lf var = %-8.2lf std = %-8.2lf\n", "", min_commtime, max_commtime, mean_commtime, 
	    var_commtime, std_commtime);
    fprintf(stderr, "IO bw / proc (MB/s) %8s min = %-8.2lf max = %-8.2lf avg = %-8.2lf var = %-8.2lf std = %-8.2lf\n", "", min_iobw, max_iobw, mean_iobw, 
	    var_iobw, std_iobw);
    fprintf(stderr, "-------------------------------\n");

  } // rank = 0

}
//-----------------------------------------------------------------------
//
// ComputeThead
//
void ComputeThread() {

//   int i, j;
//   int done;
//   int rank;
//   int first; // first time
//   int usec; // initial wait time (microseconds)
//   int ngroups; // number of groups of blocks
//   int g; // current group
//   int sb, eb; // starting and ending block in the current group
//   int bg; // number of blocks per group, except perhaps last group
//   int last_round; // last non-null round

//   // init
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   ngroups = (int)(ceil(ntpart / tr)); // number of groups
//   bg = (int)(floor(nblocks / ngroups)); // number of blocks per group, except last

//   // clear all blocks' compute status
//   for (j = 0; j < nblocks; j++)
//     lat->ClearComp(j);

//   // for all groups
//   for (g = 0; g < ngroups; g++) {

// #ifdef DEBUG
//     if (rank == 0)
//       fprintf(stderr, "** begin group %d **\n", g);
// #endif

//     // starting and ending blocks in the group
//     sb = g * bg;
//     eb = (g == ngroups - 1 ? nblocks : sb + bg);

//     last_round = 0;

//     // for all rounds
//     for (i = 0; i < max_rounds; i++) {

// #ifdef DEBUG
//       if (rank == 0)
// 	fprintf(stderr, "begin round %d\n", i);
// #endif

//       done = 0;

//       // until all blocks are computed
//       usec = 1000;
//       while (!done) {

// 	// for all blocks
// 	for (j = 0; j < nblocks; j++) {

// // 	  // blocks that are not in this group do a null send to neighbors
// // 	  if (i < sb || i >= eb) {
// // 	    lat->SendNeighbors(i);
// // #pragma omp flush
// // 	    lat->SetComp(j, i);
// // #pragma omp flush
// // 	  }

// 	  // blocks that are in this group get computed
// 	  if (!lat->GetComp(j, i) && lat->GetLoad(j)) {
// 	    ComputeFieldlines(j);
// #ifdef DEBUG
// 	    fprintf(stderr, "Computed block %d\n", j);
// #endif
// #pragma omp flush
// 	    lat->SetComp(j, i);
// #pragma omp flush
// 	  }

// 	} // for all blocks

// 	// check if all done
// 	done = 1;
// 	for (j = 0; j < nblocks; j++) {
// 	  if (!lat->GetComp(j, i))
// 	    done = 0;
// 	}

// 	if (!done) {
// 	  usleep(usec);
// 	  usec *= 2;
// 	}

//       } // for all blocks

//     if (lat->ExchangeNeighbors(Seeds, SizeSeeds))
// 	last_round = i;

//     } // for all rounds

//     TotRounds = last_round + 1;

// #ifdef DEBUG
//     if (rank == 0)
//       fprintf(stderr, "Completed %d rounds\n", max_rounds);
// #endif

//   } // for all groups

// #ifdef DEBUG
//   if (rank == 0)
//     fprintf(stderr, "Completed %d groups\n", ngroups);
// #endif

// #ifdef GRAPHICS

//   GatherFieldlines();

// #endif

}
//-----------------------------------------------------------------------
//
// IOThread
//
void IOThread() {

//   int min_t, max_t; // subdomain time bounds
//   float from[3], to[3]; // seed points limits
//   int num_loaded; // number of blocks loaded into memory so far
//   int rank;
//   int first = 1; // first time
//   int usec; // wait time in microseconds
//   int i, j, k;
//   int ngroups; // number of groups of blocks
//   int g; // current group
//   int sb, eb; // starting and ending block in the current group
//   int bg; // number of blocks per group, except perhaps last group

//   // init
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   num_loaded = 0;
//   ngroups = (int)(ceil(ntpart / tr)); // number of groups
//   bg = (int)(floor(nblocks / ngroups)); // number of blocks per group, except last

//   // init all blocks
//   for (i = 0; i < nblocks; i++) {

//     lat->ClearLoad(i);
//     osuflow[i] = new OSUFlow;
//     lat->GetVB(i, from, to, &min_t, &max_t);

//     // init seeds for blocks at t = initial time
//     if (min_t == 0) {

//       osuflow[i]->SetRandomSeedPoints(from, to, tf); 
//       seeds = osuflow[i]->GetSeeds(NumSeeds[i]); 

//       while (SizeSeeds[i] < NumSeeds[i] * sizeof(VECTOR4)) {
// 	Seeds[i] = (VECTOR4 *)realloc(Seeds[i], SizeSeeds[i] * 2);
// 	assert(Seeds[i] != NULL);
// 	SizeSeeds[i] *= 2;
//       }

//       for (k = 0; k < NumSeeds[i]; k++)
// 	Seeds[i][k].Set(seeds[k][0], seeds[k][1], seeds[k][2], min_t);

//     } // init seeds

//   } // init all blocks

//   // for all groups
//   for (g = 0; g < ngroups; g++) {

// #ifdef DEBUG
//     if (rank == 0)
//       fprintf(stderr, "** begin group %d **\n", g);
// #endif

//     // starting and ending blocks in the group
//     sb = g * bg;
//     eb = (g == ngroups - 1 ? nblocks : sb + bg);

//     // for all rounds
//     for (j = 0; j < max_rounds; j++) {

//       // for all blocks
//       for (i = 0; i < nblocks; i++) {

// 	// block is either not in this group or is loaded already
// 	if (i < sb || i >= eb || lat->GetLoad(i))
// 	  continue;

// 	// make room for the next block
// 	if (num_loaded >= b_mem) {
// 	  MultiThreadEvictBlock(j);
// 	  num_loaded--;
// 	}

// 	// read the data
// 	lat->GetVB(i, from, to, &min_t, &max_t);
// 	osuflow[i]->LoadData(filename, from, to, size,
// 			     min_t, max_t, data_mode); 

// 	// update status
// #ifdef DEBUG
// 	fprintf(stderr, "Loaded block %d\n", i);
// #endif
// #pragma omp flush
// 	lat->SetLoad(i);
// #pragma omp flush
// 	num_loaded++;
// 	osuflow[i]->ScaleField(10.0); // improves visibility

//       } // for all blocks

//       usleep(1000);

//     } // for all rounds

//   } // for all groups

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
  fprintf(stderr, "Evicted block %d\n", target);

  target = (target + 1) % nblocks;

}
//-----------------------------------------------------------------------
//
// IOandCompute
//
void IOandCompute() {

  int min_t, max_t; // subdomain temporal bounds
  float from[3], to[3]; // seed points limits
  int num_loaded; // number of blocks loaded into memory so far
  int rank;
  int i, j, k;
  int ngroups; // number of groups of blocks
  int g; // current group
  int sb, eb; // starting and ending block in the current group
  int bg; // number of blocks per group, except perhaps last group
  int last_round; // last non-null round
  double t0; // temporary timer
  float TotDataRead = 0; // total size of data read in MB

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  num_loaded = 0; // number of blocks currently loaded
  ngroups = (int)(ceil(ntpart / tr)); // number of groups
  bg = (int)(floor(nblocks / ngroups)); // number of blocks per group, except last

  TotCompTime = MPI_Wtime();

  // init all blocks
  for (i = 0; i < nblocks; i++) {

    lat->ClearLoad(i);
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

  // for all groups
  for (g = 0; g < ngroups; g++) {

#ifdef DEBUG
    if (rank == 0)
      fprintf(stderr, "** begin group %d **\n", g);
#endif

    // starting and ending blocks in the group
    sb = g * bg;
    eb = (g == ngroups - 1 ? nblocks : sb + bg);

    last_round = 0;

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

#ifdef DEBUG
      if (rank == 0)
	fprintf(stderr, " * begin round %d *\n", j);
#endif

      // for all blocks
      for (i = 0; i < nblocks; i++) {

// 	// blocks that are not in this group do a null send to neighbors
// 	if (i < sb || i >= eb) {
// 	  lat->SendNeighbors(i);
// 	  continue;
// 	}

	// if the block needs to be loaded
	if (!lat->GetLoad(i)) {

// 	  // make room for the next block
// 	  if (num_loaded >= b_mem) {
// 	    SingleThreadEvictBlock();
// 	    num_loaded--;
// 	  }

	  // read the data
	  lat->GetVB(i, from, to, &min_t, &max_t);
	  t0 = MPI_Wtime();
	  osuflow[i]->LoadData(filename, from, to, size,
			       min_t, max_t, data_mode); 
	  TotIOTime += (MPI_Wtime() - t0);
	  TotDataRead += ((to[0] - from[0]) * (to[1] - from[1]) *
			  (to[2] - from[2]) * 12); // 3 comps * 4 bytes
	  lat->SetLoad(i);
	  num_loaded++;
#ifdef DEBUG
	  fprintf(stderr, "Loaded block %d\n", i);
#endif
	  osuflow[i]->ScaleField(10.0); // improves visibility


	} // if the block needs to be loaded

	// compute pathlines
	ComputeFieldlines(i);
#ifdef DEBUG
	fprintf(stderr, "Computed block %d\n", i);
#endif

      } // for all blocks

      if (lat->ExchangeNeighbors(Seeds, SizeSeeds, NumSeeds))
	last_round = j;

    } // for all rounds

    TotRounds = last_round + 1;

#ifdef DEBUG
    if (rank == 0)
      fprintf(stderr, "Completed %d rounds\n", max_rounds);
#endif

  } // for all groups

  TotCompTime = MPI_Wtime() - TotCompTime - TotIOTime - lat->GetMyCommTime();
  TotDataRead /= 1048576; // convert to MB
  TotIOBW = TotDataRead / TotIOTime;

#ifdef DEBUG
  if (rank == 0)
    fprintf(stderr, "Completed %d groups\n", ngroups);
#endif

  // gather fieldlines for rendering
  GatherFieldlines();

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
  fprintf(stderr, "Evicted block %d\n", target);

  target = (target + 1) % nblocks;
  
}
//-----------------------------------------------------------------------
//
// ComputeFieldlines
//
// block_num: local block number (0 to nblocks-1)
// not global partition number
//
//
void ComputeFieldlines(int block_num) {

  if (tsize > 1)
    ComputePathlines(block_num);
  else
    ComputeStreamlines(block_num);

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

    TotSeeds += NumSeeds[block_num];

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

    TotSeeds += NumSeeds[block_num];

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
// gathers all fieldlines at the root for rendering
//
void GatherFieldlines() {

  static int *ntrace = NULL; // number of traces for each proc
  int n; // total number of my points
  int i;

#ifdef GRAPHICS

  // gather number of points in each trace at the root
  n = GatherNumPts(ntrace, 0);
  
  // gather the actual points in each trace at the root
  GatherPts(ntrace, n);

  // write a file too
  WriteFieldlines(ntrace, n, (char *)"field_lines.out");

#else

  // gather number of points in each trace to everyone
  n = GatherNumPts(ntrace, 1);
  
  // write the traces collectively
  WriteFieldlines(ntrace, n, (char *)"field_lines.out");

#endif

}
//-----------------------------------------------------------------------
//
// GatherNumPts
//
// gathers number of points in each trace to the root
//
// ntrace: number of traces in each process (passed by reference)
// all: all = 0 gather to root, all = 1 gather to all
//
// returns: total number of points in my process
//
int GatherNumPts(int* &ntrace, int all) {

  int myntrace = 0; // my number of traces
  static int *ofst = NULL; // offsets into ntrace
  int *mynpt; // number of points in each of my traces
  int tot_mynpt = 0; // total number of my points
  int rank, nproc; // MPI usual
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // allocate memory
  if (ntrace == NULL)
    assert((ntrace = new int[nproc]) != NULL);
  if (ofst == NULL)
    assert((ofst = new int[nproc]) != NULL);

  // compute number of my traces
  for (i = 0; i < nblocks; i++)
    myntrace += sl_list[i].size();

  // gather number of traces
  MPI_Allgather(&myntrace, 1, MPI_INT, ntrace, 1, MPI_INT, MPI_COMM_WORLD);

  // compute number of points in each of my traces
  assert((mynpt = new int[myntrace]) != NULL);
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

  // gather number of points in each trace
  tot_ntrace = 0;
  for (i = 0; i < nproc; i++) {
    ofst[i] = (i == 0) ? 0 : ofst[i - 1] + ntrace[i - 1];
    tot_ntrace += ntrace[i];
  }
  assert((npt = new int[tot_ntrace]) != NULL);
  MPI_Allgatherv(mynpt, myntrace, MPI_INT, npt, ntrace, ofst, MPI_INT,
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
  int tot_npt = 0; // total number of points in all traces from everyone

  // init
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nflt == NULL)
    assert((nflt = new int[nproc]) != NULL);
  if (ofst == NULL)
    assert((ofst = new int[nproc]) != NULL);
  assert((mypt = new VECTOR4[mynpt]) != NULL);

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

  }
  for(i = 0; i < tot_ntrace; i++)
    tot_npt += npt[i];
  assert((pt = new VECTOR4[tot_npt]) != NULL);
  MPI_Gatherv(mypt, mynpt * 4, MPI_FLOAT, pt, nflt, ofst,
	      MPI_FLOAT, 0, MPI_COMM_WORLD);

  delete[] mypt;

}
//-----------------------------------------------------------------------
//
// WriteFieldlines
//
// writes field lines
//
// ntrace: number of traces in each process
// mynpt: total number of points in my process
// filename: output file name
//
void WriteFieldlines(int *ntrace, int mynpt, char *filename) {

  MPI_File fd;
  MPI_Status status;
  int myproc;
  int ofst; // offset into the file (bytes)
  int pts_ofst = 0; // number of points before mine
  float *mypt; // points in my traces
  VECTOR4 temp; // temporary point
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  float min[4], max[4]; // extents
  int delim = -1; // delimits numbers of points from the points in the file
  int i, j, n;
    
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  assert(MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | 
		       MPI_MODE_WRONLY, MPI_INFO_NULL, &fd) == MPI_SUCCESS);

  if (myproc == 0) {

    // write extents
    min[0] = min[1] = min[2] = min[3] = 0.0;
    max[0] = size[0]; max[1] = size[1]; max[2] = size[2]; max[3] = tsize - 1;
    MPI_File_write(fd, min, 4, MPI_FLOAT, &status);
    assert(status.count == 4 * sizeof(float)); // bytes
    MPI_File_write(fd, max, 4, MPI_FLOAT, &status);
    assert(status.count == 4 * sizeof(float)); // bytes
	   
    // write numbers of points in each trace
    assert(MPI_File_write(fd, npt, tot_ntrace, MPI_INT, 
			  &status) == MPI_SUCCESS);
    assert(status.count == tot_ntrace * sizeof(int));

    // write delimiter
    assert(MPI_File_write(fd, &delim, 1, MPI_INT, &status) == MPI_SUCCESS);

  }

  // set file pointer to start of my points
  ofst = 8 * sizeof(float) + (tot_ntrace + 1) * sizeof(int);
  n = 0;
  for (i = 0; i < myproc; i++) {
    for (j = 0; j < ntrace[i]; j++)
      pts_ofst += npt[n++];
  }
  ofst += pts_ofst * 4 * sizeof(float); // pts before mine

  // collect my points in a buffer
  assert((mypt = (float *)malloc(mynpt * 4 * sizeof(float))) != NULL);
  n = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
	   pt_iter++) {
	temp = **pt_iter;
	mypt[4 * n]     = temp[0];
	mypt[4 * n + 1] = temp[1];
	mypt[4 * n + 2] = temp[2];
	mypt[4 * n + 3] = temp[3];
	n++;
      }
    }
  }

  // write my points
  MPI_File_set_view(fd, ofst, MPI_FLOAT, MPI_FLOAT, (char *)"native", 
		    MPI_INFO_NULL);
  assert(MPI_File_write_all(fd, mypt, mynpt * 4, MPI_FLOAT, &status)
	 == MPI_SUCCESS);
  assert(status.count == mynpt * 4 * sizeof(float)); // in bytes

  free(mypt);
  MPI_File_close(&fd);

}
//-----------------------------------------------------------------------
//
// GetArgs
//
// gets command line args
//
void GetArgs(int argc, char *argv[]) {

  assert(argc >= 16);

  strncpy(filename,argv[1],sizeof(filename));
  size[0] = atof(argv[2]);
  size[1] = atof(argv[3]);
  size[2] = atof(argv[4]);
  tsize = atoi(argv[5]);
  nspart = atoi(argv[6]); // total space partitions
  ntpart = atoi(argv[7]); // total time partitions
  tr = atoi(argv[8]); // number of time partitions per round
  tf = atoi(argv[9]); // traces per block
  pf = atoi(argv[10]); // points per trace
  max_rounds = atoi(argv[11]); // rounds
  threads = atoi(argv[12]) <= 1 ? 1 : 2; // threads per process
  avail_mem = atoi(argv[13]); // memory data size (MB)
  strncpy(part_file, argv[14], sizeof(part_file));
  data_mode = atoi(argv[15]);

}
//-----------------------------------------------------------------------
//
// Init
//
// inits the app
//
void Init() {

  int myproc, nproc; // usual MPI
  int nt; // max total number of traces
  int np; // max total number of points
  int ghost = 1; // number of ghost cells per spatial edge
  int b_size; // data size in a typical block (Mbytes)
  int ngr; // number of groups of rounds
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  assert(nspart * ntpart >= nproc);
  assert(tr <= ntpart);

  // init lattice and osuflow
  if (part_file[0] == '!')
    lat = new Lattice4D((int)size[0], (int)size[1], (int)size[2], tsize, ghost, 
			nspart, ntpart, nproc, myproc);
  else
    lat = new Lattice4D(part_file, (int)size[0], (int)size[1], (int)size[2], 
			tsize, ghost, nspart, ntpart, nproc, myproc);

  lat->InitSeedLists(); 
  nblocks = lat->GetNumPartitions(myproc);
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

  // max number of time steps in any block
  max_bt = (int)(ceil(tsize / ntpart) + 2 * ghost);

  // number of blocks to keep in memory
//   b_size = (int)(size[0] * size[1] / 1048576.0f * size[2] / nspart * 
// 		 tsize / ntpart * 4 * sizeof (float));
//   b_mem = avail_mem / b_size;

  // performance stats
  assert((block_stats = (int *)malloc(n_block_stats * 
				      sizeof(int))) != NULL);
  assert((time_stats = (double *)malloc(n_time_stats * 
					sizeof(double))) != NULL);

  fflush(stderr);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  // print some of the args
#ifdef DEBUG
  if (myproc == 0) {
    fprintf(stderr,"Volume size: X %.3lf Y %.3lf Z %.3lf t %d\n",
	    size[0], size[1], size[2], tsize);
    fprintf(stderr, "Number of threads per process: %d\n", threads);
    fprintf(stderr, "Number of compute rounds: %d\n", max_rounds);
    fprintf(stderr, "Available dataset memory per process: %d MB\n", 
	    avail_mem);
    fprintf(stderr, "Number of blocks a process can fit in memory: %d\n", 
	    b_mem);
  }
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
