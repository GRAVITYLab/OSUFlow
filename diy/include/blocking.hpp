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

// full tree node
struct kd_node_t {
  bb_t bounds; // bounds of region correponding to this tree node
  int proc; // process to which this node is block is or will be assigned
  int l_child; // array index of left child, -1 = leaf
  int r_child; // array index of right child, -1 = leaf
  int parent; // array index of parent, -1 = root
};

class Blocking {

 public:

  Blocking(int start_b, int did, int dim, int tot_b, int *data_size, 
	   bool share_face, int *ghost, int *given, Assignment *assignment, 
	   MPI_Comm comm);
  Blocking(int start_b, int did, int dim, int tot_b, int *gids, bb_t *bounds,
	   Assignment *assignment, MPI_Comm comm);
  ~Blocking();

  void ComputeBlocking(int *given);
  int *GetBlockSizes() { return block_size; }
  int *GetLatSizes() { return lat_size; }
  int GetDim() { return dim; }
  int64_t BlockStartsSizes(int lid, int *starts, int *sizes);
  int64_t BlockSizes(int lid, int *sizes);
  void BlockStarts(int lid, int *starts);
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
  void NumLatBlocks(int *lat_nblocks);
  int Lid2Gid(int lid);
  int Gid2Lid(int gid); // only for gids on this process
  void BuildTree(float *pts, int loc_num_pts, int glo_num_pts,
		 int num_levels, int num_bins);
  int SearchTree(float *pt, int start_node = 0);
  void GetLeaf(int index, leaf_t *leaf);

private:

  void FactorDims(int *given);
  void ApplyGhost();
  int BinarySearch(int start, int num_vals, int *vals, int target);
  void AddChildren(int parent, int split_dir, float split_frac);

  int did; // domain id
  int dim; // number of dimensions
  int start_b; // starting gid of this domain (num blocks in prior domains)
  int tot_b; // total number of blocks in the domain
  gb_t *blocks; // my local blocks
  bb_t *rbb; // no-ghost bounds of my local blocks
  int nb; // my local number of blocks
  int rank; // my MPI process rank
  int groupsize; // MPI groupsize
  float data_min[DIY_MAX_DIM]; // global data minimum
  float data_max[DIY_MAX_DIM]; // global data maximum
  int data_size[DIY_MAX_DIM]; // number of vertices in each dimension
  int lat_size[DIY_MAX_DIM]; // number of blocks in each dimension
  int block_size[DIY_MAX_DIM]; //  size of blocks in each dimension
  MPI_Comm comm; // MPI communicator
  Assignment *assign; // assignment class
  bool share_face; // whether neighboring blocks share a common face or not
  int ghost[2 * DIY_MAX_DIM]; // ghost layer per side
  vector <kd_node_t> kd_tree; // KD tree implemented as a vector

};

// callback functions are not class members
static void KdTree_MergeHistogram(char **items, int *gids, int num_items, 
				  int *hdr);
static char *KdTree_CreateHistogram(int *hdr);
static void KdTree_DestroyHistogram(void *item);
static void KdTree_CreateHistogramType(void *item, DIY_Datatype *dtype, 
					int *hdr);

#endif
