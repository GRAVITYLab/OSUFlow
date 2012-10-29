//---------------------------------------------------------------------------
//
// example of using DIY to perform blocking, assignment of blocks to
//  processes, and reading blocks from storage
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
#include "assignment.hpp"
#include "blocking.hpp"
#include "bil.h"

int main(int argc, char **argv) {

  int dim = 3;
  int tot_blocks = 8;
  int64_t data_size[3] = {10, 10, 10};
  int64_t min[3], max[3], size[3]; // block extents
  char filename[] = "test.dat"; // a pre-made 10x10x10 grid of ints
  int nblocks; // my local number of blocks
  int maxblocks; // maximum number of local blocks any process has
  int rank; // MPI process

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // create the blocking and default assignment
  // note in the blocking call that we are not adding extra ghost cells, but we
  // are sharing boundaries between blocks (share_face = 1)
  int64_t given[4] = {0, 0, 0, 0};
  Assignment *assignment = new Assignment(tot_blocks, nblocks, maxblocks,
					 MPI_COMM_WORLD);
  Blocking *blocking = new Blocking(dim, tot_blocks, data_size, 1, 0, 0, given, 
				    assignment, MPI_COMM_WORLD);

  // allocate pointers to data, in this example, the data type is int
  // the memset to 0 is needed to tell BIL to allocate the memory for us
  int *data[nblocks];
  memset(data, 0, sizeof(int*) * nblocks);

  // read blocks and print block bounds
  BIL_Init(MPI_COMM_WORLD);
  for (int i = 0; i < nblocks; i++) { // for all my blocks
    blocking->BlockStartsSizes(i, min, size);
    // convert from int64_t to int, and also reverse order for BIL
    // we are working to fix this discrepancy between BIL and blocking
    int bil_data_size[3] = { data_size[2], data_size[1], data_size[0] };
    int bil_min[3] = { min[2], min[1], min[0] };
    int bil_size[3] = { size[2], size[1], size[0] };
    // post a read for the block
    BIL_Add_block_raw(dim, bil_data_size, bil_min, bil_size, 
		      filename, MPI_INT, (void**)&(data[i]));
    // print the block bounds
    for (int j = 0; j < 3; j++)
      max[j] = min[j] + size[j] - 1;
    fprintf(stderr, "process rank = %d block local id = %d min = [%lld %lld %lld] "
	    "max = [%lld %lld %lld] size = [%lld %lld %lld]\n", 
	    rank, i, min[0], min[1], min[2], max[0], max[1], max[2],
	    size[0], size[1], size[2]);
  }

  // read all the blocks that were posted
  BIL_Read();

  // print the data values in the blocks
  for (int i = 0; i < nblocks; i++) {
    fprintf(stderr, "Data values for block %d:\n", i);
    for (int j = 0; j < blocking->TotalBlockSize(i); j++)
      fprintf(stderr, "%d ", data[i][j]);
    fprintf(stderr, "\n");
  }

  MPI_Finalize();

  return 0;

}
