//---------------------------------------------------------------------------
//
// blocking class
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

#ifndef _BLOCKING
#define _BLOCKING

#include <vector>
#include <mpi.h>
#include <math.h>
#include "diy.h"
#include "assignment.hpp"

using namespace std;

class Blocking {

 public:

  Blocking(int dim, int tot_b, int64_t *data_size, bool share_face,
	   int ghost, int ghost_dir, int64_t *given, Assignment *assignment,
	   MPI_Comm comm);
  Blocking(int dim, int tot_b, int64_t *data_size, bool share_face,
	   int ghost, int* ghost_dir, int* ghost_dim, int64_t *given, 
	   Assignment *assignment, MPI_Comm comm);
  ~Blocking();

  void Init(int dim, int tot_b, int64_t *data_size, bool share_face,
	    int ghost, int* ghost_dir, int* ghost_dim, int64_t *given,
	    Assignment *assignment, MPI_Comm comm);
  void ComputeBlocking(bool share_face, int ghost, int* ghost_dir, 
		       int* ghost_dim, int64_t *given);
  int64_t *GetBlockSizes() { return block_size; }
  void GetBlockBounds(int lid, int64_t* from, int64_t* to);
  void GetRealBlockBounds(int lid, int64_t* from, int64_t* to);
  void GetTimeBounds(int g, int64_t* tmin, int64_t* tmax);
  void GetRealTimeBounds(int g, int64_t* tmin, int64_t* tmax);
  int64_t *GetLatSizes() { return lat_size; }
  int64_t BlockStartsSizes(int lid, int64_t *starts, int64_t *sizes);
  int64_t BlockSizes(int lid, int64_t *sizes);
  void BlockStarts(int lid, int64_t *starts);
  int64_t TotalBlockSize(int lid);
  void GetNeighbors(int lid, vector<struct gb_t>& neighbors);
  void Gid2Indices(int gid, int& i, int& j);
  void Gid2Indices(int gid, int& i, int& j, int& k);
  void Gid2Indices(int gid, int& i, int& j, int& k, int &l);
  int Indices2Gid(int i, int j);
  int Indices2Gid(int i, int j, int k);
  int Indices2Gid(int i, int j, int k, int l);
  int Pt2NeighGid(int gid, float *pt, int ghost, int ghost_dir);
  bool IsIn(float *pt, int *bi, int ghost, int ghost_dir);
  bool InTimeBlock(int g, int lid, int tsize, int tb);

  Assignment *assign; // assignment class

private:

  void FactorDims(int64_t *given);
  void ApplyGhost(int ghost, int* ghost_dir, int* ghost_dim);

  int dim; // number of dimensions
  int tot_b; // total number of blocks in the domain
  int nb; // my local number of blocks
  int rank; // my MPI process rank
  int groupsize; // MPI groupsize
  int64_t data_size[MAX_DIM]; // number of vertices in each dimension
  int64_t lat_size[MAX_DIM]; // number of blocks in each dimension
  int64_t block_size[MAX_DIM]; //  size of blocks in each dimension
  MPI_Comm comm; // MPI communicator

  bb_t *bb_list;  // block bounds list, includes ghost cells, only local blocks
  bb_t *rbb_list; // real block bounds list, does not include ghost cells.
		  // currently stores all blocks.

  tb_t *tb_list;  // time bounds list, includes ghost cells
  tb_t *rtb_list; // real time bounds list, does not include ghost cells


};

#endif
