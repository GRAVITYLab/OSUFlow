//---------------------------------------------------------------------------
//
// example of using DIY to perform a time-varying analysis consisting of 
//  blocking & process assignment, then during each time block,
//  parallel reading of data, local computation, and 
//  nearest neighbor communication of results,
//  After all time blocks are executed, parallel writing of analyis results
//
// Tom Peterka
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
#include "diy.h"
#include "assignment.hpp"
#include "blocking.hpp"
#include "io.hpp"
#include "neighborhoods.hpp"
#include "bil.h"

//
// delete local block lid (local ID)
//
void DeleteBlock(int lid) {

  // free memory here for this data block

}
//
// computation for my local block number lid (local ID)
//
void Compute(int lid, Blocking *blocking, Assignment *assignment,
	     Neighborhoods *nbhds) {

  // local computation here, producing an item that needs to be sent to 
  // a neighboring block

  // in this case, we are sending a 4D point
  float p[4] = { 0.0, 0.0, 0.0, 0.0 };

  // global id of my block
  int gid = assignment->RoundRobin_lid2gid(lid);

  // global id of neighboring block where p should go
  int neigh_gid = blocking->Pt2NeighGid(gid, p, 0, 0);

  // enqueue the item for sending to neighbor
  if (neigh_gid >= 0)
    nbhds->EnqueueItem(lid, (char *)p, 4 * sizeof(float), neigh_gid);

}
//
// makes MPI datatype for receiving one item
//
// cts: pointer to counts message
//
// side effects: allocates MPI datatype
//
// returns: pointer to MPI datatype
//
MPI_Datatype* RecvItemType(int *cts) {

  MPI_Datatype *dtype = new MPI_Datatype;
  MPI_Type_contiguous(4, MPI_FLOAT, dtype);

  return dtype;

}
//
// makes an MPI datatype for sending one item
//
// cts: pointer to counts message
// pts: pointer to points message
//
// side effects: allocates MPI datatype
//
// returns: pointer to MPI datatype
//
//
MPI_Datatype* SendItemType(int *cts, char** pts) {

  MPI_Datatype *dtype = new MPI_Datatype; // datatype for one point
  MPI_Type_contiguous(4, MPI_FLOAT, dtype);

  return dtype;

}
//
// main
//
int main(int argc, char **argv) {

  int dim = 4; // number of dimensions in the problem
  int space_blocks = 8; // number of spatial blocks
  int time_blocks = 2; // number of temporal blocks
  int tot_blocks = space_blocks * time_blocks; // total number of blocks
  int64_t data_size[4] = {10, 10, 10, 4}; // data size 10x10x10, 4 timesteps
  int64_t min[4], size[4]; // block extents
  int time_steps = data_size[3]; // number of timesteps
  int64_t given[4] = {0, 0, 0, 0}; // constraints on blocking (none so far)
  given[3] = time_blocks; // impose a constraint on the time blocking
  char *infiles[] = { (char *)"test0.dat", (char *)"test1.dat", 
		      (char *)"test2.dat", (char *)"test3.dat" };
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
  Neighborhoods *nbhds = 
    new Neighborhoods(blocking, assignment, MPI_COMM_WORLD);

  // initialize reading data, in this example, the data type is int
  int *data[nblocks];
  for (int i = 0; i < nblocks; i++)
    data[i] = new int[blocking->TotalBlockSize(i)];
  BIL_Init(MPI_COMM_WORLD);

  // analysis items
  vector<vector<char *> > items; // received items from neighbors
                                 // generic pointer to a byte, not a string

  // do the analysis
  for (int g = 0; g < time_blocks; g++) {  // for all time blocks

    // delete old blocks
    if (g > 0) {
      for (int i = 0; i < nblocks; i++) {
	if (blocking->InTimeBlock(g, i, time_steps, time_blocks))
	  DeleteBlock(i);
      }
    }

    // load blocks for this time block
    for (int i = 0; i < nblocks; i++) { // for all my blocks
      if (blocking->InTimeBlock(g, i, time_steps, time_blocks)) {
	blocking->BlockStartsSizes(i, min, size);
	// convert from int64_t to int, and also reverse order for BIL
	// we are working to fix this discrepancy between BIL and blocking
	int bil_data_size[3] = { data_size[2], data_size[1], data_size[0] };
	int bil_min[3] = { min[2], min[1], min[0] };
	int bil_size[3] = { size[2], size[1], size[0] };
	// post a BIL read for each time step in this block
	int *p = data[i];
	for (int j = min[3]; j < min[3] + size[3]; j++) {
	  fprintf(stderr, "group = %d block %d min = %d %d %d "
		  "size = %d %d %d timestep = %d\n", 
		  g, i, bil_min[2], bil_min[1], bil_min[0], bil_size[2], 
		  bil_size[1], bil_size[0], j);
	  BIL_Add_block_raw(3, bil_data_size, bil_min, bil_size, 
			    infiles[j], MPI_INT, (void **)&p);
	  p += size[0] * size[1] * size[2];
	}
      }
    }
    BIL_Read();
    fprintf(stderr, "Read group %d successfully\n", g);
    MPI_Barrier(MPI_COMM_WORLD); // everyone synchronizes after reading data

    // do the computation for the current blocks in memory
    for (int i = 0; i < nblocks; i++) {
      if (blocking->InTimeBlock(g, i, time_steps, time_blocks))
	Compute(i, blocking, assignment, nbhds);
    }

    // exchange neighbors
    nbhds->ExchangeNeighbors(items, 1.0, &RecvItemType, &SendItemType);

    // flush any remaining messages
    // if multiple rounds of compute / exchange neighbors, call FlushNeighbors
    // only once after those rounds complete
    nbhds->FlushNeighbors(items, &RecvItemType);

  } // for all time blocks

  // write the results
  MPI_Datatype *dtype = new MPI_Datatype; // datatype for output
  MPI_Type_contiguous(256, MPI_FLOAT, dtype);
  MPI_Type_commit(dtype);
  MPI_Barrier(MPI_COMM_WORLD); // everyone synchronizes again
  io->WriteAnaInit(outfile, false);
  // this example isn't quite done yet, not ready to write out
//   io->WriteAllAna((void **)&items[0], space_blocks, space_blocks, dtype);
  io->WriteAnaFinalize();

  // cleanup
  MPI_Finalize();

  return 0;

}
