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

// defines
#define ERROR 0

// function prototypes
void Config(int argc, char *argv[]);
void Error(const char *fmt, ...);

// globals
static char filename[256]; // datset file name
static VECTOR3 minLen, maxLen; // data bounds

//---------------------------------------------------------------------------------

main(int argc, char *argv[]) {

  VECTOR3 minB, maxB; // volume bounds
  volume_bounds_type *vb_list;  // volume bounds list
  int nproc;  // mpi groupsize
  int rank; // mpi rank
  float from[3], to[3]; // seed points limits
  int nSeeds; // number of seeds
  VECTOR3* seeds; // seeds
  list<vtListSeedTrace*> list; // seed trace
  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  std::list<VECTOR3*>::iterator pnIter; 
  VECTOR3 p; // for debug purposes

  printf("Initializing...\n"); 

  // initialize mpi
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Barrier(MPI_COMM_WORLD);

  // get args, create classes
  Config(argc, argv);
  OSUFlow *osuflow = new OSUFlow(); 

  // subdivide the entire domain into nproc subdomains
  int* lattice;
  int lattice_xdim, lattice_ydim, lattice_zdim; 
  vb_list = calc_subvolume(maxLen[0] - minLen[0], maxLen[1] - minLen[1], 
     maxLen[2] - minLen[2], 2, nproc, &lattice, lattice_xdim, lattice_ydim, 
     lattice_zdim); 

//     printf("PE %d:  %d %d %d : %d %d %d\n", i, vb_list[i].xmin,  vb_list[i].ymin,  vb_list[i].zmin, 
// 	   vb_list[i].xmax,  vb_list[i].ymax,  vb_list[i].zmax); 

  // compute boundary
  minB[0] = vb_list[rank].xmin;  
  minB[1] = vb_list[rank].ymin;
  minB[2] = vb_list[rank].zmin; 
  maxB[0] = vb_list[rank].xmax;
  maxB[1] = vb_list[rank].ymax;
  maxB[2] = vb_list[rank].zmax; 

  // read data
  osuflow->LoadData(filename, true, minB, maxB); 
  if (rank == 0)
    printf("read file %s\n", filename); 

  // drop seeds
  from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
  to[0] = maxB[0];   to[1] = maxB[1];   to[2] = maxB[2]; 
  osuflow->SetRandomSeedPoints(from, to, 20); 
  seeds = osuflow->GetSeeds(nSeeds); 

  // debug
  if (rank == 0) {
    for (int j=0; j<nSeeds; j++) 
      printf(" seed no. %d : [%f %f %f]\n", j, seeds[j][0], seeds[j][1], seeds[j][2]); 
  }

  // integrate
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreamLines(list , FORWARD_DIR, 50, 0); 
  //     printf(" domain %d done integrations\n", i); 
  //     printf("list size = %d\n", list.size()); 

  // iterate through the result
  pIter = list.begin(); 
  for (; pIter!=list.end(); pIter++) {

    trace = *pIter; 
    pnIter = trace->begin(); 

    // debug
    for (; pnIter!=trace->end(); pnIter++) {
      p = **pnIter; 
      //      printf(" %f %f %f ", p[0], p[1], p[2]); 
    }

  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}
//---------------------------------------------------------------------------------
//
// Config
//
// gets command line args
//
void Config(int argc, char *argv[]) {

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 5)
    Error("Error: Insufficient number of command line args\n");

  strncpy(filename,argv[1],sizeof(filename));

  // hard code the minimum corner to be 0,0,
  // need to allow for variable data origin in the future
  minLen[0] = minLen[1] = minLen[2] = 0.0;

  maxLen[0] = atof(argv[2]);
  maxLen[1] = atof(argv[3]);
  maxLen[2] = atof(argv[4]);

  if (rank == 0)
    printf("Volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
	   minLen[0], maxLen[0], minLen[1], maxLen[1], minLen[2], maxLen[2]); 

}
//--------------------------------------------------------------------------------
//
// Error()
// mpi error handler
//
void Error(const char *fmt, ...){

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  sleep(5);
  MPI_Abort(MPI_COMM_WORLD,ERROR);

}
//--------------------------------------------------------------------------
