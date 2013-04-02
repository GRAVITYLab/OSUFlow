//---------------------------------------------------------------------------
//
// swap class
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

#ifndef _SWAP
#define _SWAP

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include "mpi.h"
#include "comm.hpp"
#include "assignment.hpp"
#include "diy.h"

class Swap {

public:
  Swap(int start_b, MPI_Comm comm);
  ~Swap(){}
  void SwapBlocks(int did, char **its, int **hdrs, int nr, int *kv, 
		  int num_elems, int *starts, int *sizes,
		  Comm *cc, Assignment *assign,
		  void (*reduce_func)(char **, int *, int, int), 
		  char* (*recv_create_func)(int *, int),
		  void (*recv_detroy_func)(void *),
		  void* (*send_type_func)(void*, MPI_Datatype*, int, int),
		  void (*recv_type_func)(void*, MPI_Datatype*, int));

private:

  inline void GetGrpPos(int cur_r, const int *kv, int gid, 
			int &g, int &p);
  inline void GetPartners(const int *kv, int cur_r, int gid, 
			  int *partners, int &grp, int &pos);
  void ReduceBlocks(int did, char** its, int cur_r, int *kv,
		    int *sz_part, Comm *cc, Assignment *assign,
		    void (*reduce_func)(char **, int *, int, int), 
		    char* (*recv_create_func)(int *, int),
		    void (*recv_detroy_func)(void *),
		    void (*recv_type_func)(void*, MPI_Datatype*, int));

  int start_b; // starting block global id, number of blocks in prior domains
  MPI_Comm comm; // communicator

};

#endif
