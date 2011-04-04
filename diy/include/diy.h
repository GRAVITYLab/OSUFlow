/*---------------------------------------------------------------------------
 *
 * diy C, C++ interface
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * (C) 2011 by Argonne National Laboratory.
 * See COPYRIGHT in top-level directory.
 *
--------------------------------------------------------------------------*/

#ifndef _DIY
#define _DIY

#include <stdint.h>
#include <assert.h>
#include "mpi.h"

/* #define MAX_ELEMENTS 256 /\* maximum number of elements in an MPI_Datatype *\/ */

#define MAX_HDR_ELEMENTS 256 /* max. number of elements in the header */
#define HDR_SIZE (MAX_HDR_ELEMENTS * sizeof(int)) /* size of header (bytes) */

/* order of block sizes for footer */
#define PROC 0 /* ordered by process */
#define BLOCK 1 /*ordered by  block global id */

#define MAX_DIM 10 /* maximum number of dimensions */

#define MAX_K 64 /* maximum k value of any round */
#define MAX_R 64 /* maximum number of rounds */

/* datatypes */
#define DIY_INT8     0 /* 1 byte signed integer */
#define DIY_UINT8    1 /* 1 byte unsigned integer */
#define DIY_INT16    2 /* 2 bytes signed integer */
#define DIY_UINT16   3 /* 2 bytes unsigned integer */
#define DIY_INT32    4 /* 4 bytes signed integer */
#define DIY_UINT32   5 /* 4 bytes unsigned integer */
#define DIY_FLOAT32  6 /* 4 bytes floating point */
#define DIY_FLOAT64  7 /* 8 bytes floating point */
#define DIY_FLOAT128 8 /* 8 bytes floating point */

/* block bounds */
struct bb_t { 
  int64_t min[MAX_DIM];
  int64_t max[MAX_DIM];
};

/* one global block */
struct gb_t {
  int gid; /* global block id */
  int proc; /* process to which block is assigned */
};

/* one local block */
struct lb_t {
  int gid; /* global block id */
  int loaded; /* whether block has has been loaded in memory */
};

/*-------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif
void DIY_begin(int dim, int64_t *min, int64_t *max, int tot_blocks, 
	       int dtype, MPI_Comm comm, MPI_Datatype *atype, 
	       int64_t num_blocks, int64_t *block_gids, int64_t *neigh_gids, 
	       int64_t *neigh_procs, struct bb_t *bb);

#ifdef __cplusplus
extern "C"
#endif
void DIY_end();

#ifdef __cplusplus
extern "C"
#endif
void DIY_InitWriteResults(const char *filename, int write_footer);

#ifdef __cplusplus
extern "C"
#endif
void DIY_FinalizeWriteResults();

#ifdef __cplusplus
extern "C"
#endif
void DIY_InitReadData(const char *filename);

#ifdef __cplusplus
extern "C"
#endif
void DIY_FinalizeReadData();

#ifdef __cplusplus
extern "C"
#endif
void DIY_ReadData(void ***data);

/*-------------------------------------------------------------------------*/

#endif
