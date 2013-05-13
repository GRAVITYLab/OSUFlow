//---------------------------------------------------------------------------
//
// example of using DIY to perform blocking and assignment of blocks to
//  processes
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
#include "io.hpp"
#include "merge.hpp"
#include "bil.h"

//
// Counts the number of data values in each bin
//
void Count(int *data, int *bins) {

}
//
// user-defined callback function for merging an array of items
//
// items: pointers to input items
// nitems: total number of input items
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: allocates resulting merged item
//
// returns: pointer to resulting merged item
//
char *ComputeMerge(char **items, int nitems) {

  // merging code

}
//
// user-defined callback function for creating a received item
//
// hdr: quantity information for allocating custom parts of the item
//  (not used in this example)
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: allocates the item
//
// returns: pointer to the item
//
char *CreateItem(int *hdr) {

  int *bins = new int[256]; // DIY will free this resource for you
  return (char *)bins;

}
//
// user-defined callback function for creating an MPI datatype for the
//   received item
//
// item: pointer to the item
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: creates & commits the MPI datatype
//
// returns: pointer to the datatype
//
MPI_Datatype *CreateType(char *item) {

  MPI_Datatype *dtype = new MPI_Datatype; // DIY will free this resource for you
  MPI_Type_contiguous(256, MPI_INT, dtype);
  MPI_Type_commit(dtype); // DIY will free this resource for you
  return dtype;

}
//
// main
//
int main(int argc, char **argv) {

  int dim = 3; // number of dimensions in the problem
  int tot_blocks = 8; // total number of blocks
  int64_t data_size[3] = {10, 10, 10}; // data size
  int64_t min[3], size[3]; // block extents
  int64_t given[3] = {0, 0, 0}; // constraints on blocking (none)
  int rounds = 1; // one round of merging
  int kvalues[1] = {2}; // k-way merging, eg 2-way merge in the only round
  char infile[] = "test.dat";
  char outfile[] = "test.out";

  int nblocks; // my local number of blocks
  int maxblocks; // maximum number of local blocks any process has

  MPI_Init(&argc, &argv);

  // create the blocking, default assignment i/o, and merge classes
  Assignment *assignment = new Assignment(tot_blocks, nblocks, maxblocks,
					 MPI_COMM_WORLD);
  Blocking *blocking = new Blocking(dim, tot_blocks, data_size, 1, 0, 0, given, 
				    assignment, MPI_COMM_WORLD);
  IO *io = new IO(dim, tot_blocks, maxblocks, MPI_COMM_WORLD);
  Merge *merge = new Merge(MPI_COMM_WORLD);

  // read data, assume integer, raw format
  int *data[nblocks];
  memset(data, 0, sizeof(int*) * nblocks); // memset tells BIL to allocate
                                           // data for us
  BIL_Init(MPI_COMM_WORLD);
  for (int i = 0; i < nblocks; i++) {
    blocking->BlockStartsSizes(i, min, size);
    // convert from int64_t to int, and also reverse order for BIL
    // we are working to fix this discrepancy between BIL and blocking
    int bil_data_size[3] = { data_size[2], data_size[1], data_size[0] };
    int bil_min[3] = { min[2], min[1], min[0] };
    int bil_size[3] = { size[2], size[1], size[0] };
    // post a read for the block
    fprintf(stderr, "block %d min = %d %d %d "
	    "size = %d %d %d\n", i, bil_min[2], bil_min[1], bil_min[0], 
	    bil_size[2], bil_size[1], bil_size[0]);
    BIL_Add_block_raw(dim, bil_data_size, bil_min, bil_size, 
		      infile, MPI_INT, (void**)&(data[i]));
  }
  BIL_Read();
  MPI_Barrier(MPI_COMM_WORLD); // everyone synchronizes after reading data

  // perform a local analysis, for example, compute a histogram
  int **bins; // histogram bins for each data block
  bins = new int*[nblocks];
  for (int b = 0; b < nblocks; b++) { // all my blocks
    bins[b] = new int[256]; // eg., 256 bins for each block
    Count(data[b], bins[b]);
  }

#if 0 // the rest of this ecample isn't quite finished yet

  // merge the analysis results
  // this could be a simplification of the analysis, eg. a probability density
  // or an image composition if the analysis is a local rendering
  float **pdf; // probability distribution function, in this example
  int nb_merged; // number of output merged blocks
  nb_merged = merge->MergeBlocks((char**)data, (int **)NULL, nblocks, 
				 (char** &)pdf[0], rounds, &kvalues[0], 
				 io, assignment, &ComputeMerge, &CreateItem, 
				 &CreateType);

  // write the results
  MPI_Datatype *dtype = new MPI_Datatype; // datatype for output
  MPI_Type_contiguous(256, MPI_FLOAT, dtype);
  MPI_Type_commit(dtype);
  MPI_Barrier(MPI_COMM_WORLD); // everyone synchronizes again
  io->WriteAnaInit(outfile, false);
  io->WriteAllAna((void **)pdf, nb_merged, nb_merged, dtype);
  io->WriteAnaFinalize();

  // cleanup
  for (int b = 0; b < nblocks; b++)
    delete[] bins[b];
  delete bins;
  MPI_Type_free(dtype);
  delete dtype;

#endif

  MPI_Finalize();

  return 0;

}
