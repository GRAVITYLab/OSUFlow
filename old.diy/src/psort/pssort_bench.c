#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#include "psort.h"

int* create_rand_nums(int num_elements) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  srand(world_rank * time(NULL));
  int* rand_nums = (int*)malloc(sizeof(int) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = rand();
  }
  return rand_nums;
}

int compare_int(const void* a, const void* b) {
  if (*(int*)a < *(int*)b) {
    return -1;
  } else if (*(int*)a > *(int*)b) {
    return 1;
  } else {
    return 0;
  }
}

int checksum(int* data, int num_elements) {
  int sum = 0;
  int i;
  for (i = 0; i < num_elements; i++) {
   sum += data[i];
  }

  int global_sum = 0;
  MPI_Reduce(&sum, &global_sum, 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
  return global_sum;
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
    int num_elements = 1048576;
    int initial_num_elements = num_elements;
    int* rand_nums = create_rand_nums(num_elements);
    int init_checksum = checksum(rand_nums, num_elements);

    MPI_Barrier(MPI_COMM_WORLD);
    double t = -MPI_Wtime();
    pssort_a(&pssort_object, (void**)&rand_nums, &num_elements, sizeof(int),
             compare_int);
    MPI_Barrier(MPI_COMM_WORLD);
    t += MPI_Wtime();
    if (world_rank == 0) {
      fprintf(stderr, "weak:%d:integers_per_proc:%d:total_GB:"
              "%lf:num_sample_sorters:%d:num_samples:%d:time:%lf\n",
               world_size, initial_num_elements,
              ((int64_t)initial_num_elements * sizeof(int) * world_size)
              / (1024.0 * 1024.0 * 1024.0), num_sample_sorters, num_samples, t);
                 
    }  
    MPI_Barrier(MPI_COMM_WORLD);

    int final_checksum = checksum(rand_nums, num_elements);

    if (num_elements > 0) {
      for (i = 1; i < num_elements; i++) {
        assert(rand_nums[i] >= rand_nums[i - 1]);
      }
    }
    MPI_Status recv_status;
    int num_recv_items;
    if (world_size > 1) {
      int my_max_num = (num_elements == 0) ? -1 : rand_nums[num_elements - 1];
      int my_min_num = (num_elements == 0) ? 0 : rand_nums[0];
      int neighbor_max_num;
      if (world_rank == 0) {
        MPI_Send(&my_max_num, (num_elements == 0) ? 0 : 1, MPI_INT,
                 world_rank + 1, 0, MPI_COMM_WORLD);
      } else if (world_rank == world_size - 1) {
        MPI_Recv(&neighbor_max_num, 1, MPI_INT, world_rank - 1, 0,
                 MPI_COMM_WORLD, &recv_status);
      } else {
        MPI_Sendrecv(&my_max_num, (num_elements == 0) ? 0 : 1, MPI_INT,
            world_rank + 1, 0, &neighbor_max_num, 1, MPI_INT,
            world_rank - 1, 0, MPI_COMM_WORLD, &recv_status);
      }
      MPI_Get_count(&recv_status, MPI_INT, &num_recv_items);
      if (world_rank != 0 && num_recv_items != 0 && num_elements != 0) {
        assert(my_min_num >= neighbor_max_num);
      }
    }

    if (world_rank == 0) {
      assert(init_checksum == final_checksum);
    }  
    free(rand_nums);
  }
  }
  /* Weak scaling done */

  pssort_finalize(&pssort_object);
  MPI_Finalize();
  return 0;
}
