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
	   int *ghost, int64_t *given, Assignment *assignment, 
	   MPI_Comm comm);
  Blocking(int dim, int *gids, bb_t *bounds,
	   Assignment *assignment, MPI_Comm comm);
  ~Blocking();

  void ComputeBlocking(int64_t *given);
  int64_t *GetBlockSizes() { return block_size; }
  int64_t *GetLatSizes() { return lat_size; }
  int GetDim() { return dim; }
  int64_t BlockStartsSizes(int lid, int64_t *starts, int64_t *sizes);
  int64_t BlockSizes(int lid, int64_t *sizes);
  void BlockStarts(int lid, int64_t *starts);
  int64_t TotalBlockSize(int lid);
  void BlockBounds(int lid, bb_t *bounds);
  void NoGhostBlockBounds(int lid, bb_t *bounds);
  void GetNeighbors(int lid, vector<struct gb_t>& neighbors);
  void Gid2Indices(int gid, int& i, int& j);
  void Gid2Indices(int gid, int& i, int& j, int& k);
  void Gid2Indices(int gid, int& i, int& j, int& k, int &l);
  int Indices2Gid(int i, int j);
  int Indices2Gid(int i, int j, int k);
  int Indices2Gid(int i, int j, int k, int l);
  bool InTimeBlock(int g, int lid, int tsize, int tb);
  void NumLatBlocks(int64_t *lat_nblocks);
  int Lid2Gid(int lid);
  int Gid2Lid(int gid);

private:

  void FactorDims(int64_t *given);
  void ApplyGhost();

  int dim; // number of dimensions
  int tot_b; // total number of blocks in the domain
  gb_t *blocks; // my local blocks
  bb_t *rbb; // no-ghost bounds of my local blocks
  int nb; // my local number of blocks
  int rank; // my MPI process rank
  int groupsize; // MPI groupsize
  float data_min[DIY_MAX_DIM]; // global data minimum
  float data_max[DIY_MAX_DIM]; // global data maximum
  int64_t data_size[DIY_MAX_DIM]; // number of vertices in each dimension
  int64_t lat_size[DIY_MAX_DIM]; // number of blocks in each dimension
  int64_t block_size[DIY_MAX_DIM]; //  size of blocks in each dimension
  MPI_Comm comm; // MPI communicator
  Assignment *assign; // assignment class
  bool share_face; // whether neighboring blocks share a common face or not
  int ghost[2 * DIY_MAX_DIM]; // ghost layer per side

};

#endif
