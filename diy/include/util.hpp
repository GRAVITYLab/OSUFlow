//---------------------------------------------------------------------------
//
// utilities
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
//----------------------------------------------------------------------------
#include <vector>
#include <mpi.h>
#include <math.h>
#include <diy.h>

#ifdef ZLIB
#include "zlib.h"
#define CHUNK 262144 // chunk size for zlib
#endif

using namespace std;

#ifndef __UTIL
#define __UTIL

void swap(char *n, int nitems, int item_size);
void swap8(char *n);
void swap4(char *n);
void swap2(char *n);
void CreateDtype(MPI_Aint addr, vector<map_block_t> *map, 
		 MPI_Datatype *type);
void CompressBlock(void* addr, MPI_Datatype dtype, MPI_Comm comm,
		   vector<unsigned char> *comp_buf, int *comp_size);
void DecompressBlockToDatatype(unsigned char *in_buf, int in_size,
			       void* addr, MPI_Datatype dtype, 
			       MPI_Comm comm);
void DecompressBlockToBuffer(unsigned  char* in_buf, int in_size, 
			     vector<unsigned char> *decomp_buf, 
			     int *decomp_size);

#endif
