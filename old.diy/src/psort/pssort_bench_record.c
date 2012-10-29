#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>

#include "psort.h"

typedef struct {
  char val[100];
} Record;

int compare_record(const void* a, const void* b) {
  int cmp = memcmp(((Record*)a)->val, ((Record*)b)->val, 12);
  return cmp;
/*
  if (cmp == 0) {
    return 0;
  } else if (cmp > 0) {
    return 1;
  } else {
    return -1;
  }
*/
}

Record* create_rand_nums(int num_elements) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  srand(world_rank * time(NULL));
  Record* rand_nums = (Record*)malloc(sizeof(Record) * num_elements);
  assert(rand_nums != NULL);
  int i;
  int val[3];
  for (i = 0; i < num_elements; i++) {
    int r = rand();
    val[0] = r;
    r = rand();
    val[1] = r;
    r = rand();
    val[2] = r;
    memcpy(rand_nums[i].val, val, 12);
  }
  return rand_nums;
}

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  Pssort_Object pssort_object;
  pssort_init(&pssort_object, MPI_COMM_WORLD);
  int i;

  /* Weak scaling start */
  int num_sample_sorters = 1;
  int num_samples = 8;
  for (num_samples = 32; num_samples <= 128; num_samples *= 2) {
  for (num_sample_sorters = 1; num_sample_sorters <= world_size;
       num_sample_sorters *= 2) {
    pssort_set_num_samples(&pssort_object, num_samples);
    pssort_set_sample_sorter_factor(&pssort_object,
                                    world_size / num_sample_sorters);
    int num_elements = 41943;
    int initial_num_elements = num_elements;
    Record* rand_nums = create_rand_nums(num_elements);
    MPI_Barrier(MPI_COMM_WORLD);
    double t = -MPI_Wtime();
    pssort_a(&pssort_object, (void**)&rand_nums, &num_elements, sizeof(Record),
             compare_record);
    MPI_Barrier(MPI_COMM_WORLD);
    t += MPI_Wtime();
    if (world_rank == 0) {
      fprintf(stderr, "weak:%d:integers_per_proc:%d:total_GB:"
              "%lf:num_sample_sorters:%d:num_samples:%d:time:%lf\n",
               world_size, initial_num_elements,
              ((int64_t)initial_num_elements * sizeof(Record) * world_size)
              / (1024.0 * 1024.0 * 1024.0), num_sample_sorters, num_samples, t);
                 
    }  
    MPI_Barrier(MPI_COMM_WORLD);

    free(rand_nums);
  }
  }
  /* Weak scaling done */

  pssort_finalize(&pssort_object);
  MPI_Finalize();
  return 0;
}
