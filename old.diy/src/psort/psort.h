#ifndef PSORT_H
#define PSORT_H 1

#ifdef __cplusplus
extern "C"
{
#endif

/* Basic functions */
int pmsort(void** data_ret, int* num_items_ret, int item_size,
           int (*compare)(const void*, const void*));

int pssort(void** data_ret, int* num_items_ret, int item_size,
           int (*compare)(const void*, const void*));

/* Advanced functions */
typedef struct {
  int num_samples; /* Num samples per proc, default = 32 */
  int sample_sorter_factor; /* Factor of procs to sort samples, default = 32 */
  MPI_Comm world_comm; /* Communicator for sorting. */
  MPI_Comm sample_gather_comm; /* Communicator for gathering samples */
  MPI_Comm sample_sort_comm; /* Communicator for sorting samples */
} Pssort_Object;

void pssort_init(Pssort_Object* pssort_object, MPI_Comm world_comm);

void pssort_set_num_samples(Pssort_Object* pssort_object, int num_samples);

void pssort_set_sample_sorter_factor(Pssort_Object* pssort_object,
                                     int sample_sorter_factor);

int pssort_a(Pssort_Object* pssort_object, void** data_ret,
             int* num_elements_ret, int element_size,
             int (*compare)(const void*, const void*));

void pssort_finalize(Pssort_Object* pssort_object);

typedef struct {
  MPI_Comm world_comm; /* Communicator for sorting. */
} Pmsort_Object;

void pmsort_init(Pmsort_Object* pmsort_object, MPI_Comm world_comm);

int pmsort_a(Pmsort_Object* pmsort_object, void** data_ret,
             int* num_elements_ret, int element_size,
             int (*compare)(const void*, const void*));

void pmsort_finalize(Pmsort_Object* pmsort_object);

#ifdef __cplusplus
}
#endif

#endif
