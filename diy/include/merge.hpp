//---------------------------------------------------------------------------
//
// merge class
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

#ifndef _MERGE
#define _MERGE

#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include "mpi.h"
#include "comm.hpp"
#include "assignment.hpp"
#include "diy.h"

class Merge {

public:

  Merge(int start_b, MPI_Comm comm);
  ~Merge(){}
  int MergeBlocks(int did, char **its, int **hdrs, 
		  int nr, int *kv, Comm *cc, Assignment *assign,
		  void (*merge_func)(char **, int *, int, int *),
		  char * (*create_func)(int *),
		  void (*destroy_func)(void *),
		  void (*type_func)(void*, MPI_Datatype*, int *));
  int AsyncMergeBlocks(int did, char **its, int **hdrs, float wf,
		       int nr, int *kv, Comm *cc, Assignment *assign,
		       void (*merge_func)(char **, int *, int, int *),
		       char * (*item_func)(int *),
		       void (*destroy_func)(void *),
		       void (*type_func)(void*, MPI_Datatype*, int *));

private:

  inline bool GetPartners(const int *kv, int cur_r, 
			  int gid, int *partners);

  int start_b; // starting block global id, number of blocks in prior domains
  MPI_Comm comm; // communicator

};

#endif
