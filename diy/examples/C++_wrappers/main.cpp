//---------------------------------------------------------------------------
//
// diy C++ example
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2011 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------
#include <mpi.h>
#include "diy.h"

int main(int argc, char **argv) {

  // set some default parameters
  int dim = 3;
  int tot_blocks = 64;
  int64_t min[3] = {0, 0, 0};
  int64_t max[3] = {128, 128, 128};

  // start MPI
  MPI_Init(&argc, &argv);

  // start DIY
  DIY_begin(dim, min, max, tot_blocks, 0, MPI_COMM_WORLD, NULL, 0, NULL, NULL, 
	    NULL, NULL);

  // do some work

  // end DIY
  DIY_end();

  // end MPI
  MPI_Finalize();

  return 0;

}
