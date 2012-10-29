#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include "psort.h"
#include "util.h"

/* This is used to determine the number of sample sorters */
#define SAMPLE_SORTER_FACTOR 32

/* This determines the amount of random samples to choose per proc */
#define NUM_SAMPLES 32

/*
  Utility function to split an amount into near-equally-sized chunks.
*/
void get_chunk_sizes(int total_items, int split_size, int* chunk_sizes) {
  int i;
  int chunk_size = ceil_pos((float)total_items / split_size);
  int num_big_chunks = split_size - ((split_size * chunk_size) - total_items);
  int num_small_chunks = split_size - num_big_chunks;
  chunk_size--;
  for (i = 0; i < num_small_chunks; i++) {
    chunk_sizes[i] = chunk_size;
  }
  chunk_size++;
  for (i = num_small_chunks; i < split_size; i++) {
    chunk_sizes[i] = chunk_size;
  }
}

/*
  Obtains the splitter elements for the bins. To obtain the splitter elements,
  the samples are first sorted in parallel and then the splitter elements are
  chosen at equal steps through the sample size per proc
*/
void get_splitters(void* data, int num_items, const int num_samples,
                   int item_size, MPI_Comm world_comm,
                   MPI_Comm sample_sort_comm, MPI_Comm sample_gather_comm,
                   int (*compare)(const void *, const void *),
                   void* splitters) {
  MPI_Datatype item_type;
  MPI_Type_contiguous(item_size, MPI_BYTE, &item_type);
  MPI_Type_commit(&item_type);
  int world_size;
  MPI_Comm_size(world_comm, &world_size);
  int world_rank;
  MPI_Comm_rank(world_comm, &world_rank);

  int *chunk_sizes = cmalloc(sizeof(int) * num_samples);
  get_chunk_sizes(num_items, num_samples, chunk_sizes);

  /* TODO obtain samples randomly */
  void* samples = cmalloc(item_size * num_samples);
  int64_t sample_offset = 0;
  int i;
  if (num_items != 0) {
    for (i = 0; i < num_samples; i++) {
      memcpy(samples + (i * item_size),
             data + sample_offset, item_size);
      sample_offset += (chunk_sizes[i] * item_size);
    }
  } else {
    memset(samples, 0, item_size * num_samples);
  }

  int sample_gather_size;
  MPI_Comm_size(sample_gather_comm, &sample_gather_size);
  int sample_gather_rank;
  MPI_Comm_rank(sample_gather_comm, &sample_gather_rank);
  /* Gather the samples to sort. First gather the amount of samples per
     proc. */
  void* gathered_samples = NULL;
  int num_gathered_samples = num_samples * sample_gather_size;
  if (sample_gather_rank == 0) {
    gathered_samples = cmalloc(num_samples * sample_gather_size * item_size);
  }
  MPI_Gather(samples, num_samples, item_type, gathered_samples,
             num_samples, item_type, 0, sample_gather_comm);
  /* Now sort the gathered samples if you are a sorter and obtain the
     splitters */
  if (sample_gather_rank == 0) {
    int num_sorters;
    MPI_Comm_size(sample_sort_comm, &num_sorters);
    /* Parallel distributed merge sort. Warning, this could potentially produce
       incorrect results if the data from each sorter is not equal. This case
       occurs when the sample gatherers do no have equal sizes. Since this
       often results in only minor noise when choosing splitter keys, it is
       ignored. Also, the splitter keys are resorted just to be sure this
       doesn't produce incorrect results in parallel sample sort */
    Pmsort_Object pmsort_object;
    pmsort_init(&pmsort_object, sample_sort_comm);
    pmsort_a(&pmsort_object, &gathered_samples, &num_gathered_samples,
             item_size, compare);
    pmsort_finalize(&pmsort_object);

   /* Choose the last sorted sample from each process. Although we will gather
      world_size amount of splitter keys, only the first world_size - 1 will be
      used as splitters. This is to allow proc 0 to obtain all the values < the
      first splitter and vice versa for the last proc */
    for (i = 0; i < sample_gather_size; i++) {
      memcpy(gathered_samples + (i * item_size), gathered_samples
             + (((num_samples * i) + num_samples - 1) * item_size), item_size);
    }

    /* If each sample gather group is equal in size, we dont have to worry
       about the case where pmsort might produce incorrect results */
    if (world_size % SAMPLE_SORTER_FACTOR == 0 ||
        world_size < SAMPLE_SORTER_FACTOR) {
      MPI_Gather(gathered_samples, sample_gather_size, item_type,
                 splitters, sample_gather_size, item_type,
                 0, sample_sort_comm);
    } else {
      /* TODO handle this case */
      assert(0);
    }
    cfree(gathered_samples);
  }

  MPI_Bcast(splitters, world_size, item_type, 0, world_comm);

  // clean up
  MPI_Type_free(&item_type);
  cfree(chunk_sizes);
  cfree(samples);
}

void bin_data(MPI_Comm world_comm, void** data_ret, int* num_items_ret,
              int item_size, int (*compare)(const void *, const void *),
              void* splitters) {
  void* data = *data_ret;
  int num_items = *num_items_ret;
  MPI_Datatype item_type;
  MPI_Type_contiguous(item_size, MPI_BYTE, &item_type);
  MPI_Type_commit(&item_type);
  int world_size;
  MPI_Comm_size(world_comm, &world_size);
  int world_rank;
  MPI_Comm_rank(world_comm, &world_rank);
  int i;

  /* Allocate the binned data and figure out the bins of each item */
  int* send_counts = cmalloc(sizeof(int) * world_size);
  memset(send_counts, 0, sizeof(int) * world_size);
  void* binned_data = cmalloc(item_size * num_items);
  int* which_bins = cmalloc(sizeof(int) * num_items);
  for (i = 0; i < num_items; i++) {
    int which_bin = (bsearchx(data + (i * item_size), splitters,
                              world_size - 1, item_size, compare) + 1) / 2;
    which_bins[i] = which_bin;
    send_counts[which_bin]++;
  }

  /* Make the offset array for MPI_Alltoallv */
  int* send_offsets = cmalloc(sizeof(int) * world_size);
  send_offsets[0] = 0;
  for (i = 1; i < world_size; i++) {
    send_offsets[i] = send_offsets[i - 1] + send_counts[i - 1];
  }

  /* Keep track of the bin items left for creating the binned data */
  int* binned_items = cmalloc(sizeof(int) * world_size);
  memset(binned_items, 0, sizeof(int) * world_size);

  /* Create the binned data */
  for (i = 0; i < num_items; i++) {
    memcpy(binned_data + ((binned_items[which_bins[i]]
           + send_offsets[which_bins[i]]) * item_size), 
           data + (i * item_size), item_size);
    binned_items[which_bins[i]]++;
  }

  cfree(data);
  cfree(which_bins);
  cfree(binned_items);

  /* Find out the sizes that other processes will contribute to your bin */
  int* recv_counts = cmalloc(sizeof(int) * world_size);
  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, world_comm);

  /* Find the size of your bin */
  int num_bin_items = 0;
  for (i = 0; i < world_size; i++) {
    num_bin_items += recv_counts[i];
  }

  /* Allocate your new data */
  void* bin = cmalloc(num_bin_items * item_size);

  /* Make the offset array for MPI_Alltoallv */
  int* recv_offsets = cmalloc(sizeof(int) * world_size);
  recv_offsets[0] = 0;
  for (i = 1; i < world_size; i++) {
    recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
  }

  /* Bin the data */
  MPI_Alltoallv(binned_data, send_counts, send_offsets, item_type,
                bin, recv_counts, recv_offsets, item_type, world_comm);

  /* Clean up */
  cfree(send_counts);
  cfree(send_offsets);
  cfree(recv_counts);
  cfree(recv_offsets);
  cfree(binned_data);
  MPI_Type_free(&item_type);

  *data_ret = bin;
  *num_items_ret = num_bin_items;
}

/*
  Takes the amount of sample sorts to use and creates the necessary
  communicators to take advantage of using a smaller amount of
  sample sorters. It returns a communicator for gathering the samples down
  to a smaller amount of processes, and then another communicator for
  sorting those samples.
*/
void init_sample_comms(Pssort_Object* pssort_object) {
  MPI_Comm world_comm = pssort_object->world_comm;
  int world_size;
  MPI_Comm_size(world_comm, &world_size);
  int world_rank;
  MPI_Comm_rank(world_comm, &world_rank);
  int num_sample_sorters =
    ceil_pos((float)world_size / pssort_object->sample_sorter_factor);

  if (world_rank < num_sample_sorters) {
    MPI_Comm_split(world_comm, 0, world_rank, 
                   &(pssort_object->sample_sort_comm));
  } else {
    MPI_Comm_split(world_comm, MPI_UNDEFINED, 0,
                   &(pssort_object->sample_sort_comm));
  }

  if (world_rank < num_sample_sorters) {
    MPI_Comm_split(world_comm, world_rank, 0,
                   &(pssort_object->sample_gather_comm));
  } else {
    int* sample_gather_sizes = cmalloc(sizeof(int) * num_sample_sorters);
    get_chunk_sizes(world_size - num_sample_sorters, num_sample_sorters, 
                    sample_gather_sizes);
    int i;
    int amount = 0;
    for (i = 0; i < num_sample_sorters; i++)
    {
      if (world_rank < sample_gather_sizes[i] + num_sample_sorters + amount) {
        MPI_Comm_split(world_comm, i, world_rank,
                       &(pssort_object->sample_gather_comm));
        break;
      }
      amount += sample_gather_sizes[i];
    }
    free(sample_gather_sizes);
  }
}

/*
  Frees communicators for gathering and sorting samples.
*/
void finalize_sample_comms(Pssort_Object* pssort_object) {
  if (pssort_object->sample_gather_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&(pssort_object->sample_gather_comm));
  }
  if (pssort_object->sample_sort_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&(pssort_object->sample_sort_comm));
  }
}

/*
  Initializes the advanced mode of using pssort.
*/
void pssort_init(Pssort_Object* pssort_object, MPI_Comm world_comm) {
  pssort_object->world_comm = world_comm;
  pssort_object->num_samples = NUM_SAMPLES;
  pssort_object->sample_sorter_factor = SAMPLE_SORTER_FACTOR;
  init_sample_comms(pssort_object);
}

/*
  Sets the number of samples each proc uses during sorting.
*/
void pssort_set_num_samples(Pssort_Object* pssort_object, int num_samples) {
  pssort_object->num_samples = num_samples;
}

/*
  Sets the factor of sample sorters to use.
*/
void pssort_set_sample_sorter_factor(Pssort_Object* pssort_object, 
                                     int sample_sorter_factor) {
  pssort_object->sample_sorter_factor = sample_sorter_factor;
  finalize_sample_comms(pssort_object);
  init_sample_comms(pssort_object);
}

/*
  Frees sample communicators.
*/
void pssort_finalize(Pssort_Object* pssort_object) {
  if (pssort_object == NULL) {
    return;
  }
  finalize_sample_comms(pssort_object);
}

/*
  Wrapper around pssort_a.
*/
int pssort(void** data_ret, int* num_items_ret, int item_size,
           int (*compare)(const void *, const void *)) {
  Pssort_Object pssort_object;
  pssort_init(&pssort_object, MPI_COMM_WORLD);
  pssort_a(&pssort_object, data_ret, num_items_ret, item_size, compare);
  pssort_finalize(&pssort_object);
}

/*
  Performs the parallel sample sort algorithm. The data given by the user
  will be set to a different buffer that will most likely have a different
  number of items than the original data. It returns 1 always, for now.
*/
int pssort_a(Pssort_Object* pssort_object, void** data_ret,
             int* num_items_ret, int item_size,
             int (*compare)(const void *, const void *)) {
  int world_size;
  MPI_Comm_size(pssort_object->world_comm, &world_size);
  /* Handle base case */
  if (world_size == 1) {
    qsort(*data_ret, *num_items_ret, item_size, compare);
    return 1;
  }

#ifdef __PSSORT_TIMING
  double splitter_sort_t = -MPI_Wtime();
#endif

  /* Find the splitter items for bins */
  void* splitters = cmalloc(item_size * world_size);
  get_splitters(*data_ret, *num_items_ret, pssort_object->num_samples, 
                item_size, pssort_object->world_comm, 
                pssort_object->sample_sort_comm,
                pssort_object->sample_gather_comm, compare, splitters);

#ifdef __PSSORT_TIMING
  splitter_sort_t += MPI_Wtime();
  double bin_t = -MPI_Wtime();
#endif

  /* Bin the data and sort the bins */
  bin_data(pssort_object->world_comm, data_ret, num_items_ret, item_size,
           compare, splitters);

#ifdef __PSSORT_TIMING
  bin_t += MPI_Wtime();
  double sort_t = -MPI_Wtime();
#endif

  qsort(*data_ret, *num_items_ret, item_size, compare);

#ifdef __PSSORT_TIMING
  sort_t += MPI_Wtime();
#endif

#ifdef __PSSORT_TIMING
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) {
    fprintf(stderr, "splitter sort %lf bin %lf sort %lf total %lf\n",
           splitter_sort_t, bin_t, sort_t, splitter_sort_t + bin_t + sort_t);
  }
#endif

  /* Clean up */
  cfree(splitters);
  return 1;
}
