#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

#include "psort.h"
#include "util.h"

/*
  Creates a schedule of mergers for processes during distributed merge sort.
*/
int* create_merge_schedule(int world_size, int world_rank,
                           int* num_merges_ret) {
  int i;

  /* Find the closest power of two to determine the size of the swap schedule */
  i = 0;
  while ((1 << i) < world_size) {
    i++;
  }
  *num_merges_ret = (1 << i) - 1;

  int* swap_schedule = (int*)cmalloc(sizeof(int) * (*num_merges_ret));
  for (i = 0; i < (*num_merges_ret); i++) {
    swap_schedule[i] = (world_rank & ~(i + 1)) | (~world_rank & (i + 1));
  }

  return swap_schedule;
}

/*
  Merges two lists backwards. Given the first list of size m and the second
  list of size n, it will merge the last n elements.
*/
void merge_items_backward(void* data, void* recv_data, void* merge_data,
                          int num_items, int num_recv_items, int item_size,
                          int (*compare)(const void*, const void*)) {
  int i, m, n;
  int data_size = num_items * item_size;
  int recv_size = num_recv_items * item_size;

  for (i = data_size - item_size, m = recv_size - item_size,
       n = data_size - item_size; i >= 0 && m >= 0 && n >= 0; i -= item_size) {
    if (compare(recv_data + m, data + n) > 0) {
      memcpy(merge_data + i, recv_data + m, item_size);
      m -= item_size;
    } else {
      memcpy(merge_data + i, data + n, item_size);
      n -= item_size;
    }
  }

  /* Copy remainder */
  if (m < 0) {
    memcpy(merge_data, data + n - i, i + item_size);
  } else if (n < 0) {
    memcpy(merge_data, recv_data + m - i, i + item_size);
  }
}

/*
  Merges two lists forwards. Given the first list of size m and the second
  list of size n, it will merge the first m elements.
*/
void merge_items_forward(void* data, void* recv_data, void* merge_data,
                         int num_items, int recv_num_items, int item_size,
                         int (*compare)(const void*, const void*)) {
  int i, m, n;
  int data_size = num_items * item_size;
  int recv_size = recv_num_items * item_size;

  /* Merge items forward, keeping num_items elements */
  for (i = 0, m = 0, n = 0; m < recv_size && n < data_size && i < data_size;
       i += item_size) {
    if (compare(recv_data + m, data + n) < 0) {
      memcpy(merge_data + i, recv_data + m, item_size);
      m += item_size;
    } else {
      memcpy(merge_data + i, data + n, item_size);
      n += item_size;
    }
  }

  /* Copy remainder */
  if (m == recv_size) {
    memcpy(merge_data + i, data + n, data_size - i);
  } else if (n == data_size) {
    memcpy(merge_data + i, recv_data + m, data_size - i);
  }

}

/*
  Initializes state for parallel merge sorting.
*/
void pmsort_init(Pmsort_Object* pmsort_object, MPI_Comm world_comm) {
  pmsort_object->world_comm = world_comm;
}

/*
  Currently does nothing since it does not allocate memory.
*/
void pmsort_finalize(Pmsort_Object* pmsort_object) {
  return;
}

/*
  Wrapper around pmsort_a.
*/
int pmsort(void** data_ret, int* num_items_ret, int item_size,
           int32_t (*compare)(const void*, const void*)) {
  Pmsort_Object pmsort_object;
  pmsort_init(&pmsort_object, MPI_COMM_WORLD);
  pmsort_a(&pmsort_object, data_ret, num_items_ret, item_size, compare);
  pmsort_finalize(&pmsort_object);
}

/*
  Performs a distributed merge sort amongst processes. The data and num_items
  that are returned to the user will be different than the data buffers they
  supplied because of the nature of how the algorithm works. This
  implementation of parallel merge sort is primarily used by the parallel
  sample sort algorithm.
*/
int pmsort_a(Pmsort_Object* pmsort_object, void** data_ret, int* num_items_ret,
             int item_size, int32_t (*compare)(const void*, const void*)) {
  MPI_Comm world_comm = pmsort_object->world_comm;
  assert(world_comm != MPI_COMM_NULL);
  int32_t world_size;
  MPI_Comm_size(world_comm, &world_size);
  int32_t world_rank;
  MPI_Comm_rank(world_comm, &world_rank);

  /* Sort data to be merged */
  qsort(*data_ret, *num_items_ret, item_size, compare);

  /* Handle base case */
  if (world_size == 1) {
    return 1;
  }

  /* For convenience */
  MPI_Datatype item_type;
  MPI_Type_contiguous(item_size, MPI_BYTE, &item_type);
  MPI_Type_commit(&item_type);
  void* data = *data_ret;
  int num_items = *num_items_ret;

  /* Build the merging schedule */
  int num_merges;
  int* merge_schedule = create_merge_schedule(world_size, world_rank,
                                              &num_merges);

  /* Find the max number of items that will be received per any stage */
  int max_num_items = 0;
  MPI_Allreduce(&num_items, &max_num_items, 1, MPI_INT, MPI_MAX, world_comm);

  /* Allocate merger receive buffer and merge buffer */
  void* recv_data = cmalloc(item_size * max_num_items);
  void* merge_data = cmalloc(item_size * num_items);
  
  int i;
  void* merger_min_max = cmalloc(item_size);
  for (i = 0; i < num_merges; i++) {
    /* For non-power-of-two process counts, the merge_schedule may contain
       invalid ranks */
    if (merge_schedule[i] < world_size) {
      int send_num_items = num_items;
      int send_offset = 0;
      MPI_Status recv_status;
      int num_recv_items;

      /* Swap mins and maxs to minimize communication */
      if (world_rank < merge_schedule[i]) {
        MPI_Sendrecv(data + ((num_items - 1) * item_size),
                     (num_items == 0) ? 0 : 1, item_type,
                     merge_schedule[i], i, merger_min_max, 1, item_type,
                     merge_schedule[i], i, world_comm, &recv_status);
        MPI_Get_count(&recv_status, item_type, &num_recv_items);
        /* If sender or receiver has no elements, go to the next merger */
        if (num_items == 0 || num_recv_items == 0) {
          continue;
        }
        /* Find the first element to send to the other merger */
        send_offset = (bsearchx(merger_min_max, data, num_items,
                                item_size, compare) + 1) / 2;
        send_num_items = num_items - send_offset;
      } else {
        MPI_Sendrecv(data, (num_items == 0) ? 0 : 1, item_type,
                     merge_schedule[i], i, merger_min_max, 1, item_type,
                     merge_schedule[i], i, world_comm, &recv_status);
        MPI_Get_count(&recv_status, item_type, &num_recv_items);
        /* If sender or receiver has no elements, go to the next merger */
        if (num_items == 0 || num_recv_items == 0) {
          continue;
        }
        send_num_items = (bsearchx(merger_min_max, data,
                                   num_items, item_size, compare) + 1) / 2;
      }

      MPI_Sendrecv(data + (send_offset * item_size), send_num_items, 
                   item_type, merge_schedule[i], i, recv_data, max_num_items,
                   item_type, merge_schedule[i], i, world_comm, &recv_status);
      MPI_Get_count(&recv_status, item_type, &num_recv_items);

      /* Merge the items on each merger, either keeping the upper or lower
         part of elements */
      if (world_rank < merge_schedule[i]) {
        merge_items_forward(data, recv_data, merge_data, num_items, 
                            num_recv_items, item_size, compare);
      } else {
        merge_items_backward(data, recv_data, merge_data, num_items, 
                             num_recv_items, item_size, compare);
      }

      /* Update your data */
      void* temp = merge_data;
      merge_data = data;
      data = temp;
    }
  }
  cfree(merger_min_max);

  /* Clean up */
  cfree(merge_schedule);
  cfree(merge_data);
  cfree(recv_data);
  MPI_Type_free(&item_type);
  
  *num_items_ret = num_items;
  *data_ret = data;
}
