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
#include "io.hpp"
#include "assignment.hpp"
#include "diy.h"

class Merge {

public:
  Merge(MPI_Comm comm);
  ~Merge(){}
  int MergeBlocks(char **mh_in, int **hdrs, int nb_in, char** &mh_out,
		  int nr, int *kv, IO *io, Assignment *assign,
		  char * (*merge_func)(char **, int),
		  char * (*item_func)(int *),
		  MPI_Datatype* (*type_func)(char *));

private:

  inline void GetGrpPos(int r, const int *kv, int gid, 
			int &g, int &p);
  inline bool GetPartners(const int *kv, int cur_r, 
			  int gid, int *partners, int grp);

  MPI_Comm comm; // communicator

};

#endif
