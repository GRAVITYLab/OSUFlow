//------------------------------------------------------------------------------
//
// mpi test main
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
#include "OSUFlow.h"

#include <list>
#include <iterator>

#include "calc_subvolume.h"
#include "Lattice.h"

#define MAX_RECV_ATTEMPTS 10
#define MAX_ITERATIONS 3

// function prototypes
void Config(int argc, char *argv[]);
void PrintSeeds(int nblocks);
void RecvToSeeds(int local_block, int global_block, Lattice *lat);
int EndTrace(list<vtListSeedTrace*> &list, VECTOR3 &p, int index);

// globals
static char filename[256]; // dataset file name
static VECTOR3 minLen, maxLen; // data bounds
VECTOR3 size; // domain size
static int bf; // number of blocks per rank
int *NumSeeds; // number of seeds
int *SizeSeeds; // size of seeds list (bytes)
VECTOR3 **Seeds; // seeds lists
OSUFlow **osuflow; // one flow object for each block

//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  VECTOR3 minB, maxB; // subdomain bounds
  volume_bounds_type *vb_list;  // volume bounds list
  int nproc;  // mpi groupsize
  int rank; // mpi rank
  float from[3], to[3]; // seed points limits
  std::list<vtListSeedTrace*>::iterator trace_iter; // iterator over seed traces
  vtListSeedTrace *trace; // single seed trace
  std::list<VECTOR3*>::iterator pt_iter; // iterator over points in one seed trace
  VECTOR3 p; // current point
  Lattice* lat; // lattice
  int lattice_xdim, lattice_ydim, lattice_zdim; // lattice bounds
  int neighbor; // neighbor's number (0-5)
  int ei, ej, ek; // neighbor's lattice position
  int i, j, k;
  int ghost; // number of ghost cells per side
  int init_seed_num = 10; // initial number of seeds
  int *blocks; // list of blocks
  int nblocks; // number of blocks
  int quit = 0;

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

  }


  // trace particles until they remain constant for some period of time
  for (j = 0; j < MAX_ITERATIONS; j++) {

    // for all blocks, integrate points and send messages
    for(i = 0; i < nblocks; i++) {

      // debug: print current seeds
//       fprintf(stderr, "Current seeds\n");
//       PrintSeeds(nblocks);

      if (NumSeeds[i]) {

	// declaration needs to be here to start with an empty list each time
	// need to learn how to clear a list w/o redeclaring it
	list<vtListSeedTrace*> list; // trace of seed points

	// perform the integration
	// todo: integrate in both directions
	osuflow[i]->SetIntegrationParams(1, 5); 
	osuflow[i]->GenStreamLines(Seeds[i], FORWARD_DIR, NumSeeds[i], 50, list); 

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

	// debug
// 	lat->PrintPost(blocks[i]);

	// send boundary list to neighbors
	lat->SendNeighbors(blocks[i], MPI_COMM_WORLD);

      } // if (NumSeeds[i])
      
    } // for all blocks

    // for all blocks, receive messages
    for (i = 0; i < nblocks; i++) {

      for (k = 0; k < MAX_RECV_ATTEMPTS; k++) {
	if (lat->ReceiveNeighbors(blocks[i], MPI_COMM_WORLD))
	  break;
	usleep(100000);
      }

      if (k == MAX_RECV_ATTEMPTS)
	continue;

      // debug
      lat->PrintRecv(blocks[i]);

      // prepare for next iteration
      RecvToSeeds(i, blocks[i], lat);

    } // for all blocks

  } // for all iterations

  // debug: print current seeds
//   fprintf(stderr, "Final seeds\n");
//   PrintSeeds(nblocks);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// Config
//
// gets command line args
//
void Config(int argc, char *argv[]) {

  int rank;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

}
//-----------------------------------------------------------------------
//
// EndTrace()
//
// returns the point at the end of the trace
// list: list of traces
// p: output point
// trace_index: which trace (0 - num_traces - 1)
//
// returns: 1 if the point is valid
//          0 if empty trace or no more traces
//
int EndTrace(list<vtListSeedTrace*> &list, VECTOR3 &p, int trace_index) {

  std::list<vtListSeedTrace*>::iterator trace_iter; // iterator through traces
  std::list<VECTOR3*>::iterator pt_iter; // iterator through points in trace

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  trace_iter = list.begin();

  for (int i = 0; i < trace_index; i++) {
    if (trace_iter == list.end())
      fprintf(stderr,"1: rank = %d\n", rank);
      return 0;
      trace_iter++;
  }

  vtListSeedTrace *trace = *trace_iter;

  if (trace->size() == 0) {
    fprintf(stderr,"2: rank = %d\n", rank);
    return 0;
  }

  pt_iter = trace->begin(); 
  for (; pt_iter != trace->end(); pt_iter++) {
    p = **pt_iter;
    fprintf(stderr,"EndTrace rank = %d p = %.3lf %.3lf %.3lf\n",rank,p[0],p[1],p[2]);
  }

  return 1;

}
//------------------------------------------------------------------------
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
// lat: lattice
//
void RecvToSeeds(int local_block, int global_block, Lattice *lat) {

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
