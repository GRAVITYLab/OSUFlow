//---------------------------------------------------------------------------
//
// diy wrappers, callable from C, C++, and Fortran
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
#include "diy.h"
#include "blocking.hpp"
#include "io.hpp"

// generic pointers to class objects
static void *io; // io object
static void *merge; // merge object
static void *blocking; // blocking object

// diy global data
static struct bb_t extents; // grid extents
static int64_t sizes[MAX_DIM]; // grid size (extents.max - extents.min + 1)
static int64_t tb; // total number of blocks in the domain
static int64_t nb; // my number of blocks
static bb_t *bb_list; // block bounds
static int dim; // number of dimensions
static int dtype; // datatype: char, unsigned char, etc.
//--------------------------------------------------------------------------
// C and C++ wrapper
void DIY_begin(int dim, int64_t *min, int64_t *max, int tot_blocks, 
	       int dtype, MPI_Comm comm, MPI_Datatype *atype, 
	       int64_t num_blocks, int64_t *block_gids, int64_t *neigh_gids, 
	       int64_t *neigh_procs, struct bb_t *bb) {

  for (int i = 0; i < dim; i++) {
    assert(max[i] > min[i]); // sanity
    extents.min[i] = min[i];
    extents.max[i] = max[i];
    sizes[i] = extents.max[i] - extents.min[i] + 1;
  }
  ::dim = dim;
  tb = tot_blocks;
  ::dtype = dtype;
//   blocking = new Blocking(dim, tot_blocks, comm);
  io = new IO(dim, tot_blocks, 1, comm);

}
// fortran wrapper
extern "C"
void diy_begin_(int dim, int64_t *extents, int *tot_blocks, int *dtype, 
		MPI_Comm comm, MPI_Datatype *atype, int64_t *block_gids, 
		int64_t *neigh_gids, int64_t *neigh_procs, 
		struct bb_t *bb) {


}
//--------------------------------------------------------------------------
void DIY_end() {

}

extern "C"
void diy_end_() {

}
//--------------------------------------------------------------------------
void DIY_InitWriteResults(char *filename, int write_footer) {
  ((IO *)io)->WriteAnaInit(filename, write_footer);
}

extern "C"
void diy_init_write_results_(char *filename, int write_footer) {
  ((IO *)io)->WriteAnaInit(filename, write_footer);
}
//--------------------------------------------------------------------------
void DIY_FinalizeWriteResults() {
  ((IO *)io)->WriteAnaFinalize();
}

extern "C"
void diy_finalize_write_results_() {
  ((IO *)io)->WriteAnaFinalize();
}
//--------------------------------------------------------------------------
void DIY_InitReadData(char *filename) {
  ((IO *)io)->ReadDataInit(filename);
}

extern "C"
void diy_init_read_data_(char *filename) {
  ((IO *)io)->ReadDataInit(filename);
}
//--------------------------------------------------------------------------
void DIY_FinalizeReadData() {
  ((IO *)io)->ReadDataFinalize();
}

extern "C"
void diy_finalize_read_data_() {
  ((IO *)io)->ReadDataFinalize();
}
//--------------------------------------------------------------------------
void DIY_ReadData(void ***data) {

  assert(dtype >= 0 && dtype < 9); // sanity
  switch(dtype) {
  case 0:
    ((IO *)io)->ReadAllData<char>((char***)data, sizes, nb, 
				  (Blocking *)blocking);
  break;
  case 1:
    ((IO *)io)->ReadAllData<unsigned char>((unsigned char***)data, sizes, nb, 
					   (Blocking *)blocking);
  break;
  case 2:
    ((IO *)io)->ReadAllData<int16_t>((int16_t***)data, sizes, nb,
				     (Blocking *)blocking);
  break;
  case 3:
    ((IO *)io)->ReadAllData<uint16_t>((uint16_t***)data, sizes, nb, 
				      (Blocking *)blocking);
  break;
  case 4:
    ((IO *)io)->ReadAllData<int32_t>((int32_t***)data, sizes, nb, 
				     (Blocking *)blocking);
  break;
  case 5:
    ((IO *)io)->ReadAllData<uint32_t>((uint32_t***)data, sizes, nb, 
				      (Blocking *)blocking);
  break;
  case 6:
    ((IO *)io)->ReadAllData<float>((float***)data, sizes, nb, 
				   (Blocking *)blocking);
  break;
  case 7:
    ((IO *)io)->ReadAllData<double>((double***)data, sizes, nb, 
				    (Blocking *)blocking);
  break;
  case 8:
    ((IO *)io)->ReadAllData<long double>((long double***)data, sizes, nb, 
					 (Blocking *)blocking);
  break;
  }

}

extern "C"
void diy_read_data_(void ***data) {

}
//--------------------------------------------------------------------------

