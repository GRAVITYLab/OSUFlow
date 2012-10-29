#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#include "psort.h"

int* create_rand_nums(int num_elements) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  srand(world_rank);
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

int checksum(int32_t* data, int32_t num_elements) {
  int32_t sum = 0;
  int32_t i;
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

  int32_t num_elements = 100;
  int32_t* rand_nums = create_rand_nums(num_elements);

  int32_t init_checksum = checksum(rand_nums, num_elements);

  // After sorting is complete, "rand_nums" will contain a new buffer of
  // numbers and "num_elements" will contain the new amount of nums in the
  // buffer.
  pssort((void**)&rand_nums, &num_elements, sizeof(int), compare_int);
//  pmsort((void**)&rand_nums, &num_elements, sizeof(int), compare_int);
  
  int final_checksum = checksum(rand_nums, num_elements);

  // Check for correctness.
  if (num_elements > 0) {
    int i;
    for (i = 1; i < num_elements; i++) {
      assert(rand_nums[i] >= rand_nums[i - 1]);
    }
  }

  // Check neighboring processes for correctness.
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
    printf("Sorting complete and correct\n");
  }  

  free(rand_nums);
  MPI_Finalize();
  return 0;
}
