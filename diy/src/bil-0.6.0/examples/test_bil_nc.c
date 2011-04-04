/***************************************************************************
 *   Copyright (C) 2010  Wes Kendall                                       *
 *   kendall@eecs.utk.edu                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Lesser General Public License for more details.                   *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// test_bil_nc.c
// Wes Kendall
// 05/29/2010
// program for testing the netcdf reader of BIL. it is assumed that the user
// has created raw test files with the create_nc_files program. BIL
// reads in blocks from these files with random block sizes. based on the
// block sizes and locations, BIL then asserts that the correct data was
// read in since the files were created with integers in increasing order
////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>
#include "bil.h"

// number of blocks per variable that each process reads
#define NUM_BLOCKS 2 
// number of variables in the netcdf files each process reads
#define NUM_VARS 2 

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    fprintf(stderr, "usage %s dim_size num_files <file_dir>\n", argv[0]);
    exit(1);
  }
  const int dim_size = atoi(argv[1]);
  const int num_files = atoi(argv[2]);
  const char *file_dir = (argc == 4) ? argv[3] : "./";

  // do some checking of arguments
  assert(strlen(file_dir) < 512);
  assert(dim_size <= 512 && dim_size > 0);
  assert(num_files <= 16 && num_files > 0);

  MPI_Init(NULL, NULL);

  // initialize BIL using MPI_COMM_WORLD
  BIL_Init(MPI_COMM_WORLD);

  int rank;
  int g, h, i, j, k, l;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int t = time(NULL);
  srand(rank + t);

  char file_names[num_files][1024];
  for (i = 0; i < num_files; i++) {
    memset(file_names[i], 0, 1024);
    sprintf(file_names[i], "%s/%d_%d_%d.%d.nc", file_dir,
        dim_size, dim_size, dim_size, i);
  }

  int block_starts[NUM_BLOCKS * NUM_VARS * num_files][3];
  int block_sizes[NUM_BLOCKS * NUM_VARS * num_files][3];
  int block = 0;
  int* block_data[NUM_BLOCKS * NUM_VARS * num_files];
  memset(block_data, 0, sizeof(int*) * NUM_BLOCKS * NUM_VARS * num_files);
  for (h = 0; h < num_files; h++) {
    for (g = 0; g < NUM_VARS; g++) {
      for (i = 0; i < NUM_BLOCKS; i++, block++) {
        // randomly generate block bounds
        block_starts[block][0] = ((float)rand() / RAND_MAX) * (dim_size - 1);
        block_starts[block][1] = ((float)rand() / RAND_MAX) * (dim_size - 1);
        block_starts[block][2] = ((float)rand() / RAND_MAX) * (dim_size - 1);
        
        block_sizes[block][0] = 
          ((float)rand() / RAND_MAX) * (dim_size - block_starts[block][0]) + 1;
        block_sizes[block][1] = 
          ((float)rand() / RAND_MAX) * (dim_size - block_starts[block][1]) + 1;
        block_sizes[block][2] = 
          ((float)rand() / RAND_MAX) * (dim_size - block_starts[block][2]) + 1;

        // add the variable blocks to BIL for reading in. BIL needs the
        // block bounds and some other information including the file and
        // variable name
        if (g == 0) {
          BIL_Add_block_nc(3, block_starts[block], 
              block_sizes[block], file_names[h], "var0",
              (void**)&(block_data[block]));
        } else if (g == 1) {
          BIL_Add_block_nc(3, block_starts[block], 
              block_sizes[block], file_names[h], "var1",
              (void**)&(block_data[block]));
        }
      }
    }
  }

  // read all the blocks. the data is stored as an array of pointers to
  // the blocks in the order they were added
  BIL_Read();

  // check all block values
  for (block = 0, h = 0; h < num_files; h++) {
    for (g = 0; g < NUM_VARS; g++) {
      for (i = 0; i < NUM_BLOCKS; i++, block++) {
        int off = 0;
        for (j = 0; j < block_sizes[block][0]; j++) {
          for (k = 0; k < block_sizes[block][1]; k++) {
            for (l = 0; l < block_sizes[block][2]; l++, off++) {
              int val = (h * dim_size * dim_size * dim_size)
                    + ((j + block_starts[block][0]) * dim_size * dim_size)
                    + ((k + block_starts[block][1]) * dim_size)
                    + block_starts[block][2] + l;
              if (g == 0) {
                if (block_data[block][off] != val) {
                  fprintf(stderr, "Error for rand seed %d, %d %d\n", 
                      t, block_data[block][off], val);
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              } else if (g == 1) {
                if (block_data[block][off] != val + 1) {
                  fprintf(stderr, "Error for rand seed %d, %d %d\n", 
                      t, block_data[block][off], val + 1);
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              }
            }
          }
        }
      }
    }
  }

  if (rank == 0) {
    printf("Test successfully completed\n");
  }

  // free the data
  for (h = 0; h < num_files * NUM_VARS * NUM_BLOCKS; h++) {
    free(block_data[h]);
  }

  // finalize BIL
  BIL_Finalize();
  MPI_Finalize();

  return 0;
}
