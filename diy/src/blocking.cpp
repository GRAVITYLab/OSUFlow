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

#include "blocking.hpp"

//----------------------------------------------------------------------------
//
// constructor for creating blocking from scratch
//
// start_b: starting block global id (number of blocks in prior domains)
// did: domain id
// dim: number of dimensions
// tot_b: total number of blocks
// data_size: data size in up to 4 dimensions
// share_face: whether neighboring blocks share a common face or are
//  separated by a gap of one unit
// ghost: ghost layer for each dimension and side (min, max)
//  each entry can be 0 (no ghost) or > 0 (this many ghost cells per side)
//   {x min side ghost, x max side ghost, y min side ghost, y max side ghost...}
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//   eg., {0, 0, 0, t} would result in t blocks in the 4th dimension
// assignment: pointer to asignment class
// comm: MPI communicator
//
Blocking::Blocking(int start_b, int did, int dim, int tot_b, 
		   int *data_size, bool share_face,
		   int *ghost, int *given, 
		   Assignment *assignment, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->start_b = start_b;
  this->did = did;
  this->dim = dim;
  this->tot_b = tot_b;
  this->comm = comm;
  this->share_face = share_face;
  assign = assignment;
  nb = assign->NumBlks();

  for (int i = 0; i < dim; i++) {
    data_min[i] = 0.0f; // for now, until data min and max are added to API
    data_max[i] = data_size[i] - 1.0f; // for now
    (this->data_size)[i] = data_size[i];
    (this->ghost)[2 * i] = ghost[2 * i];
    (this->ghost)[2 * i + 1] = ghost[2 * i + 1];
  }

  // save local blocks list
  blocks = new gb_t[nb];
  rbb = new bb_t[nb];
  ComputeBlocking(given);

#pragma omp parallel for
  for (int i = 0; i < nb; i++) {
    blocks[i].gid = assign->AssignGid(i);
    blocks[i].proc = -1; // no need to store my own process id
  }

}
//----------------------------------------------------------------------------
//
// constructor for adopting an existing blocking
//
// start_b: starting block global id (number of blocks in prior domains)
// did: domain id
// dim: number of dimensions
// tot_b: total number of blocks
// gids: global ids of my local blocks (unique across all domains)
// bounds: block bounds (extents) of local blocks
// assignment: pointer to asignment class
// comm: MPI communicator
//
Blocking::Blocking(int start_b, int did, int dim, int tot_b, int *gids, 
		   bb_t *bounds, Assignment *assignment, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->start_b = start_b;
  this->did = did;
  this->dim = dim;
  this->tot_b = tot_b;
  this->comm = comm;
  assign = assignment;
  nb = assign->NumBlks();

  for (int i = 0; i < dim; i++)
    data_size[i] = -1; // unknown

  blocks = new gb_t[nb];

#pragma omp parallel for
  for (int i = 0; i < nb; i++) {
    blocks[i].gid = gids[i]; // user's gids need to be unique across all domains
    blocks[i].proc = -1; // no need to store my own process id
    for (int j = 0; j < dim; j++) {
      blocks[i].bb.min[j] = bounds[i].min[j];
      blocks[i].bb.max[j] = bounds[i].max[j];
    }
  }

}
//----------------------------------------------------------------------------
//
// destructor
//
Blocking::~Blocking() {

  delete[] blocks;
  if (data_size[0] != -1) // new blocking was computed
    delete[] rbb;

}
//----------------------------------------------------------------------------
//
// number of blocks in each dimension
//
// nblocks: number of blocks in each dimension
//
void Blocking::NumLatBlocks(int *lat_nblocks) {

  for (int i = 0; i < dim; i++)
    lat_nblocks[i] = lat_size[i];

}
//----------------------------------------------------------------------------
//
// block starts and sizes
// for blocks consisting of discrete, regular grid points
//
// lid: local block id
// starts: pointer to allocated array of starting block extents (output), index
//  of starting grid point (not cell) in each direction
// sizes: pointer to allocated array of block sizes (output), number of grid
//  points (not cells) in each direction
//
// returns: total size of block (product of sizes in each dimension)
//
int64_t Blocking::BlockStartsSizes(int lid, int *starts, int *sizes) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++) {
    starts[i] = blocks[lid].bb.min[i];
    sizes[i] = blocks[lid].bb.max[i] - starts[i] + 1;
    tot_size *= sizes[i];
  }

  return tot_size;

}
//--------------------------------------------------------------------------
//
// block sizes
// for blocks consisting of discrete, regular grid points
//
// lid: local block id
// sizes: pointer to allocated array of block sizes (output), number of grid
//  points (not cells) in each direction
//
// returns: total size of block (product of sizes in each dimension)
//
int64_t Blocking::BlockSizes(int lid, int *sizes) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++) {
    sizes[i] = blocks[lid].bb.max[i] - blocks[lid].bb.min[i] + 1;
    tot_size *= sizes[i];
  }

  return tot_size;

}
//--------------------------------------------------------------------------
//
// block starts
// for blocks consisting of discrete, regular grid points
//
// lid: local block id
// starts: pointer to allocated array of starting block extents (output), index 
//  of starting grid point (not cell) in each direction
//
void Blocking::BlockStarts(int lid, int *starts) {

  for (int i = 0; i < dim; i++)
    starts[i] = blocks[lid].bb.min[i];

}
//--------------------------------------------------------------------------
//
// total block size
// for blocks consisting of discrete, regular grid points
//
// lid: local block id
//
// returns: total size of block, product of number of grid pints (not cells) 
//  in each dimension
//
int64_t Blocking::TotalBlockSize(int lid) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++)
    tot_size *= (int64_t)(blocks[lid].bb.max[i] - blocks[lid].bb.min[i] + 1.0f);

  return tot_size;

}
//--------------------------------------------------------------------------
//
// block bounds, including ghost
// for blocks consisting of continuous regions
//
//
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
//
void Blocking::BlockBounds(int lid, bb_t *bounds) {

  for (int i = 0; i < dim; i++) {
    bounds->min[i] = blocks[lid].bb.min[i];
    bounds->max[i] = blocks[lid].bb.max[i];
  }

}
//--------------------------------------------------------------------------
//
// block bounds, excluding ghost
// for blocks consisting of continuous regions
//
//
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
void Blocking::NoGhostBlockBounds(int lid, bb_t *bounds) {

  // existing decomposition w/ unknown ghost, just return original bounds
  if (data_size[0] == -1) {
    for (int i = 0; i < dim; i++) {
      bounds->min[i] = blocks[lid].bb.min[i];
      bounds->max[i] = blocks[lid].bb.max[i];
    }
  }

  // otherwise return bounds w/o ghost
  else {
    for (int i = 0; i < dim; i++) {
      bounds->min[i] = rbb[lid].min[i];
      bounds->max[i] = rbb[lid].max[i];
    }
  }

}
//--------------------------------------------------------------------------
//
// calculates a default blocking of the data
// by factoring the total number of verices in alternating x, y, z, t directions
// supports 2D, 3D, and 4D depending on number of dimensions when
// constructing the object
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//   eg., {0, 0, 0, t} would result in t blocks in the 4th dimension
//
void Blocking::ComputeBlocking(int *given) {

  // factor data dimensions
  FactorDims(given);

  // volume bounds in row-major index order (x, y, z)
  // x changes fastest, z slowest
  int lid = 0;
  int gid = start_b;
  int nl = (dim < 4 ? 1 : lat_size[3]);
  int nk = (dim < 3 ? 1 : lat_size[2]);
  int nj = (dim < 2 ? 1 : lat_size[1]);

  // new version
  int gap = (share_face ? 0 : 1);
  int cur[4]; // current location in the lattice
  int d[4] = {0, 0, 0, 0}; // actual current block size
  cur[3] = 0;
  for (int l = 0; l < nl; l++) {
    cur[2] = 0;
    for (int k = 0; k < nk; k++) {
      cur[1] = 0;
      for (int j = 0; j < nj; j++) {
	cur[0] = 0;
	for (int i = 0; i < (int)lat_size[0]; i++) {

	  // x
	  if (lat_size[0] * block_size[0] > data_size[0] &&
	      (i + 1) * block_size[0] + 
	      (lat_size[0] - i - 1) * (block_size[0] - 1) >=
	      data_size[0])
	    d[0] = block_size[0] - 1;
	  else
	    d[0] = block_size[0];

	  // only store my local blocks
	  if (assign->Gid2Proc(gid) == rank) {
	    blocks[lid].bb.min[0] = cur[0];
	    if (i == (int)lat_size[0] - 1)
	      blocks[lid].bb.max[0] = data_size[0] - 1;
	    else
	      blocks[lid].bb.max[0] = cur[0] + d[0] - gap;
	  }

	  // y
	  if (lat_size[1] * block_size[1] > data_size[1] &&
	      (j + 1) * block_size[1] + 
	      (lat_size[1] - j - 1) * (block_size[1] - 1) >=
	      data_size[1])
	    d[1] = block_size[1] - 1;
	  else
	    d[1] = block_size[1];

	  if (assign->Gid2Proc(gid) == rank) {
	    blocks[lid].bb.min[1] = cur[1];
	    if (j == (int)lat_size[1] - 1)
	      blocks[lid].bb.max[1] = data_size[1] - 1;
	    else
	      blocks[lid].bb.max[1] = cur[1] + d[1] - gap;
	  }

	  // z
	  if (dim > 2) {
	    if (lat_size[2] * block_size[2] > data_size[2] &&
		(k + 1) * block_size[2] + 
		(lat_size[2] - k - 1) * (block_size[2] - 1) >=
		data_size[2])
	      d[2] = block_size[2] - 1;
	    else
	      d[2] = block_size[2];

	    if (assign->Gid2Proc(gid) == rank) {
	      blocks[lid].bb.min[2] = cur[2];
	      if (k == (int)lat_size[2] - 1)
		blocks[lid].bb.max[2] = data_size[2] - 1;
	      else
		blocks[lid].bb.max[2] = cur[2] + d[2] - gap;
	    }
	  }

	  // t
	  if (dim > 3) {
	    if (lat_size[3] * block_size[3] > data_size[3] &&
		(l + 1) * block_size[3] + 
		(lat_size[3] - l - 1) * (block_size[3] - 1) >=
		data_size[3])
	      d[3] = block_size[3] - 1;
	    else
	      d[3] = block_size[3];

	    if (assign->Gid2Proc(gid) == rank) {
	      blocks[lid].bb.min[3] = cur[3];
	      if (l == (int)lat_size[3] - 1)
		blocks[lid].bb.max[3] = data_size[3] - 1;
	      else
		blocks[lid].bb.max[3] = cur[3] + d[3] - gap;
	      if (data_size[3] == 1) // special case for 4D but static
		blocks[lid].bb.max[3] = 0;
	    }
	  }

	  if (assign->Gid2Proc(gid) == rank)
	    lid++;
	  gid++;
	  cur[0] += d[0];
	}
	cur[1] += d[1];
      }
      cur[2] += d[2];
    }
    cur[3] += d[3];
  }

  // ghost cells (also stores no-ghost bounds)
  ApplyGhost();

  // debug: print the bb list
//   for (int i = 0; i < nb; i++) {
//     gid = assign->Lid2Gid(i);
//     fprintf(stderr, "lid = %d gid = %d  min = [%.1f %.1f %.1f %.1f] "
// 	    "max = [%.1f %.1f %.1f %.1f]\n", 
// 	    i, gid, blocks[i].bb.min[0], blocks[i].bb.min[1], blocks[i].bb.min[2],
// 	    blocks[i].bb.min[3], blocks[i].bb.max[0], blocks[i].bb.max[1], 
// 	    blocks[i].bb.max[2], blocks[i].bb.max[3]);
//   }

}
//---------------------------------------------------------------------------
//
// factors the total number of verices in alternating x, y, z, t directions
//
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//
void Blocking::FactorDims(int *given) { 

  int rem = tot_b; // unfactored remaining portion of tot_b
  int block_dim[DIY_MAX_DIM]; // current block size
  int max; // longest remaining direction (0, 1, 2)
  int i, j;

  // init
  for (i = 0; i < dim; i++) {
    if (given[i] == 0) {
      lat_size[i] = 1;
      block_dim[i] = data_size[i];
    }
    else {
      lat_size[i] = given[i];
      if (rem % given[i])

#ifdef MAC
	fprintf(stderr,"Unable to block the volume with given[%d] = %d "
		"dimension. Please provide different 'given' constraints and "
		"rerun.\n", i, given[i]);
#else
	fprintf(stderr,"Unable to block the volume with given[%d] = %d "
		"dimension. Please provide different 'given' constraints and "
		"rerun.\n", i, given[i]);
#endif

      assert(rem % given[i] == 0);
      rem /= given[i];
    }
  }

  // compute factorization of data dimensions into lattice dimensions
  while (1) {

    // find longest division direction
    max = 0;
    for(i = 1; i < dim; i++) {
      if (given[i] == 0 && block_dim[i] > block_dim[max])
	max = i;
    }

    // smallest factor remaining gets assigned to this direction
    for (j = 2; j <= rem; j++) {
      if (rem % j == 0) {
	lat_size[max] *= j;
	block_dim[max] /= j;
	rem /= j;
	break;
      }
    }

    if (rem == 1)
      break;

    if (j > rem)
      fprintf(stderr,"Unable to block the volume into %d blocks. "
	      "Please select a different number of blocks and rerun.\n", tot_b);
    assert(j <= rem);

  }

  // sanity check
  int prod_blocks = 1;
  for (i = 0; i < dim; i++)
    prod_blocks *= lat_size[i];
  assert(prod_blocks == tot_b);

  // block sizes
  for(i = 0; i < dim; i++)
    block_size[i] = (int)(roundf((float)data_size[i] / (float)lat_size[i]));

  // debug
//   fprintf(stderr, "block sizes = [ ");
//   for (i = 0; i < dim; i++)
//     fprintf(stderr, "%d ", block_size[i]);
//   fprintf(stderr, "]\n");

}
//---------------------------------------------------------------------------
//
// copies original bounds and adds the ghost layer to the block bounds
//
void Blocking::ApplyGhost() {

  // copy original bounds
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j < dim; j++) {
	rbb[i].min[j] = blocks[i].bb.min[j];
	rbb[i].max[j] = blocks[i].bb.max[j];
    }
  }

  // add ghost
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j < dim; j++) {
      if (blocks[i].bb.min[j] - ghost[2 * j] >= data_min[j])
	blocks[i].bb.min[j] -= ghost[2 * j];
      if (blocks[i].bb.max[j] + ghost[2 * j + 1] <= data_max[j])
	blocks[i].bb.max[j] += ghost[2 * j + 1];
    }
  }

}
//---------------------------------------------------------------------------
//
// the following functions only work for a round robin assignment
// we still need to generalize these to process order assingment, and
// perhaps some of the utility functions should be part of the assignment
// class instead of the blocking class?
//
//---------------------------------------------------------------------------
//
// gets all neighbor blocks of a local block
// neighbors include the block itself
// neighbor info includes the global id as well as the process rank
//
// lid: my local block number
// neighbors: global block ids and process ids of neighbor blocks
//
void Blocking::GetNeighbors(int lid, vector<struct gb_t>& neighbors) {

  int mi, mj, mk, ml; // my lattice coords
  int i, j, k, l; // offset lattice coords to get from me to neighbor
  int gid; // my neighbor blocks's global id
  gb_t neigh; // one neighbor

  switch(dim) {
  case 2:
    Gid2Indices(Lid2Gid(lid), mi, mj);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	gid = Indices2Gid(mi + i, mj + j);
	if (gid >= 0) {
	  neigh.gid = gid;
	  neigh.proc = assign->Gid2Proc(gid);
	  neighbors.push_back(neigh);
	}
      }
    }
    break;
  case 3:
    Gid2Indices(Lid2Gid(lid), mi, mj, mk);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	for (k = -1; k <= 1; k++) {
	  gid = Indices2Gid(mi + i, mj + j, mk + k);
	  if (gid >= 0) {
	    neigh.gid = gid;
	    neigh.proc = assign->Gid2Proc(gid);
	    neighbors.push_back(neigh);
	  }
	}
      }
    }
    break;
  case 4:
    Gid2Indices(Lid2Gid(lid), mi, mj, mk, ml);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	for (k = -1; k <= 1; k++) {
	  for (l = -1; l <= 1; l++) {
	    gid = Indices2Gid(mi + i, mj + j, mk + k, ml + l);
	    if (gid >= 0) {
	      neigh.gid = gid;
	      neigh.proc = assign->Gid2Proc(gid);
	      neighbors.push_back(neigh);
	    }
	  }
	}
      }
    }
    break;
  }
  
}
//---------------------------------------------------------------------------
//
// gets lattice indices given a global block id (2D version)
//
// gid: global block id
// i, j: lattice indices (output)
//
void Blocking::Gid2Indices(int gid, int& i, int& j) {

  j = (gid - start_b) / lat_size[0];
  i = (gid - start_b) % lat_size[0];

}
//---------------------------------------------------------------------------
//
// gets lattice indices given a global block id (3D version)
//
// gid: global block id
// i, j, k: lattice indices (output)
//
void Blocking::Gid2Indices(int gid, int& i, int& j, int& k) {

  k = (gid - start_b) / (lat_size[0] * lat_size[1]) ; 
  j = ((gid - start_b) % (lat_size[0] * lat_size[1])) / lat_size[0]; 
  i = (gid - start_b) % lat_size[0]; 

}
//---------------------------------------------------------------------------
//
// gets lattice indices given a global block id (4D version)
//
// gid: global block id
// i, j, k, l: lattice indices (output)
//
void Blocking::Gid2Indices(int gid, int& i, int& j, int& k, int &l) {

  l = (gid - start_b) / (lat_size[0] * lat_size[1] * lat_size[2]); 
  int r = (gid - start_b) % (lat_size[0] * lat_size[1] * lat_size[2]); 
  k = r / (lat_size[0] * lat_size[1]) ; 
  j = (r % (lat_size[0] * lat_size[1])) / lat_size[0]; 
  i = r % lat_size[0]; 

}
//---------------------------------------------------------------------------
//
// return the global id of a block given its lattice coords (2D version)
// returns -1 if out of the domain
//
int Blocking::Indices2Gid(int i, int j) {

  if (i < 0 || i >= (int)lat_size[0] || j < 0 || j >= (int)lat_size[1]) 
    return(-1); 

  return(start_b + j * lat_size[0] + i);

}
//----------------------------------------------------------------------------
//
// return the global id of a block given its lattice coords (3D version)
// returns -1 if out of the domain
//
int Blocking::Indices2Gid(int i, int j, int k) {

  if (i < 0 || i >= (int)lat_size[0] || j < 0 || j >= (int)lat_size[1] || 
      k < 0 || k >= (int)lat_size[2]) 
    return(-1); 

  return(start_b + k * lat_size[0] * lat_size[1] + j * lat_size[0] + i);

}
//----------------------------------------------------------------------------
//
// return the global id of a block given its lattice coords (4D version)
// returns -1 if out of the domain
//
int Blocking::Indices2Gid(int i, int j, int k, int l) {

  if (i < 0 || i >= (int)lat_size[0] || j < 0 || j >= (int)lat_size[1] || 
      k < 0 || k >= (int)lat_size[2] || l < 0 || l >= (int)lat_size[3]) 
    return(-1); 

  return(start_b + l * lat_size[0] * lat_size[1] * lat_size[2] + 
	 k * lat_size[0] * lat_size[1] + j * lat_size[0] + i);

}
//----------------------------------------------------------------------------
//
// tests whether block b is in time block g
//
// g: current time block
// lid: local block id
// tsize: total number of timesteps
// tb: total number of global time blocks
//
bool Blocking::InTimeBlock(int g, int lid, int tsize, int tb) {

  int starts[4]; // block starts

  BlockStarts(lid, starts);
  if (tsize == 1 || tb == 1 || starts[3] == g * tsize / tb)
    return true;

  return false;

}
//-----------------------------------------------------------------------
//
// local block id to global block id
//
int Blocking::Lid2Gid(int lid) {

    return(blocks[lid].gid);

}
//----------------------------------------------------------------------------
//
// global block id to local block id
// only for global blocks on this process, see assign class for all gids
//
// returns -1 if gid not found
//
int Blocking::Gid2Lid(int gid) {

  int i;

  for (i = 0; i < nb; i++) {
    if (blocks[i].gid == gid)
      break;
  }

  if (i < nb)
    return i;
  else
    return -1;

}
//----------------------------------------------------------------------------
//
// build kd tree (prototype)
//
// pts: point locations to be indexed in kd-tree (for now)
// loc_num_pats: local number of points
// glo_num_pats: global number of points
// num_levels: number of tree levels, counting root
// num_bins: number of histogram bins at all levels
//
void Blocking::BuildTree(float *pts, int loc_num_pts, int glo_num_pts,
			 int num_levels, int num_bins) {

  int num_hists; // number of historgrams in this level
  int median = glo_num_pts / 2; // desired median
  int parent = 0; // curent parent tree node
  int tot_num_bins; // total number of bins in all histograms for a block

  // headers for each block
  int **hdrs = new int*[nb];
  for (int i = 0; i < nb; i++)
    hdrs[i] = new int[1];

  // todo: only the right number of bins for global range
  // yet to implement local ranges, especially for later levels

  // allocate histograms, hists[i] allocated and freed during each level
  int **hists;
  hists = new int*[nb];

  // initialize kd-tree with level 0
  kd_node_t node;
  for (int i = 0; i < dim; i++) {
    node.bounds.min[i] = data_min[i];
    node.bounds.max[i] = data_max[i];
  }
  node.proc = 0; // default
  node.l_child = -1; // empty, to be filled in later
  node.r_child = -1; // ditto
  node.parent = -1; // will remain empty for root
  kd_tree.push_back(node);

  for (int level = 1; level < num_levels; level++) { // tree levels

    int dir = (level - 1) % dim;
    num_hists = ((level - 1) ? num_hists * 2 : 1);
    // toto: reduce histogram size with level (after one full cycle of all dirs)
    // turned off for now
//     if (level - 1)
//       num_bins = (num_bins >= 2 * min_num_bins ? num_bins / 2 : min_num_bins);
    tot_num_bins = num_hists * num_bins;

    for (int b = 0; b < nb; b++)
      hists[b] = new int[tot_num_bins];

    // init histograms
    for (int b = 0; b < nb; b++) {
      for (int i = 0; i < tot_num_bins; i++)
	hists[b][i] = 0;
    }

    for (int b = 0; b < nb; b++) { // local blocks

      // scan and bin objects
      // todo: is data_max and data_min set correctly? I doubt it for particles
      float bin_width = (data_max[dir] - data_min[dir]) / num_bins;
      for (int i = 0; i < loc_num_pts; i++) {

	// search for particle in tree
	// todo: don't have to start at root each time
	int pt_node = SearchTree(&pts[3 * i], 0);
	// pos of node in level
	int temp = (int)(pow(2, level - 1));
	int node_level_pos = pt_node + 1 - temp;

	// debug
// 	fprintf(stderr, "pt: %.1f %.1f %.1f in node %d: node_level_pos %d\n",
// 		pts[3 * i], pts[3 * i + 1], pts[3 * i + 2],
// 		pt_node, node_level_pos);

	int bin = pts[3 * i + dir] / bin_width;
	if (bin >= num_bins)
	  bin = num_bins - 1;
	bin += num_bins * node_level_pos; // move to correct histogram
	hists[b][bin]++;

      } // objects

      hdrs[b][0] = tot_num_bins;

    } // local blocks

    // debug: print the histograms
//     for (int b = 0; b < nb; b++) {
//       for (int i = 0; i < num_bins; i++)
// 	fprintf(stderr, "hists[%d][%d] = %d\n", b, i, hists[b][i]);
//     }

    // merge the histograms
    // todo: change merge and swap API to take a target k and figure out
    // rounds and kvalues itself
    int rounds = log2f((float)tot_b); // todo: assumes power of 2 blocks
    int kvalues[rounds];
    for (int i = 0; i < rounds; i++)
      kvalues[i] = 2;
    int nb_merged; // number of output merged blocks

    DIY_Merge_blocks(did, (char**)hists, hdrs, rounds, kvalues, 
		     &KdTree_MergeHistogram, &KdTree_CreateHistogram, 
		     &KdTree_DestroyHistogram, 
		     &KdTree_CreateHistogramType, &nb_merged);

    // find median split points in histograms
    int split_index[num_hists]; // split indices in cumulative mass function
    if (rank == groupsize - 1) {

      assert(nb_merged == 1); // sanity

      // debug: print the merged histogram
//       for (int i = 0; i < num_bins; i++)
// 	fprintf(stderr, "hist[%d] = %d\n", i, hists[0][i]);

      // convert histogram to cumulative mass function; prefix sum
      for (int i = 0; i < num_hists; i++) {
	int ofst = i * num_bins; // start of this histogram
	for (int j = 1; j < num_bins; j++)
	  hists[0][ofst + j] += hists[0][ofst + j - 1];
      }

      // debug: print the CMF
//       for (int i = 0; i < num_bins; i++)
// 	fprintf(stderr, "cmf[%d] = %d\n", i, hists[0][i]);

      // find split index of CMF
      for (int i = 0; i < num_hists; i++) {
	int ofst = i * num_bins; // start of this histogram
	split_index[i] = BinarySearch(ofst, num_bins, hists[0], median);
	// debug
	fprintf(stderr, "level = %d split_index[%d] = %d median = %d\n", 
		level, i, split_index[i], median);
      }

    }

    // broadcast the split points
    MPI_Bcast(split_index, num_hists, MPI_INT, groupsize - 1, comm);

    // add new level to kd-tree
    for (int i = 0; i < num_hists; i++) {
      int ofst = i * num_bins; // start of this histogram
      AddChildren(parent + i, dir, ((float)split_index[i] - ofst)/ num_bins);
    }

    parent += num_hists;
    median /= 2;

    // cleanup
    for (int b = 0; b < nb; b++)
      delete[] hists[b];

  // debug: print the kd tree
//   if (rank == 0) { // duplicated on all ranks, print only once
//     for (int i = 0; i < kd_tree.size(); i++)
//       fprintf(stderr, "kd tree node %d: proc %d min[%.1f %.1f %.1f] "
// 	      "max[%.1f %.1f %.1f] l_child %d r_child %d parent %d\n",
// 	      i, kd_tree[i].proc, kd_tree[i].bounds.min[0], 
// 	      kd_tree[i].bounds.min[1], kd_tree[i].bounds.min[2], 
// 	      kd_tree[i].bounds.max[0], kd_tree[i].bounds.max[1], 
// 	      kd_tree[i].bounds.max[2], kd_tree[i].l_child, kd_tree[i].r_child,
// 	      kd_tree[i].parent);
//   }

  } // tree levels

  // debug: print the kd tree
  if (rank == 0) { // duplicated on all ranks, print only once
    for (int i = 0; i < (int)kd_tree.size(); i++)
      fprintf(stderr, "kd tree node %d: proc %d min[%.1f %.1f %.1f] "
	      "max[%.1f %.1f %.1f] l_child %d r_child %d parent %d\n",
	      i, kd_tree[i].proc, kd_tree[i].bounds.min[0], 
	      kd_tree[i].bounds.min[1], kd_tree[i].bounds.min[2], 
	      kd_tree[i].bounds.max[0], kd_tree[i].bounds.max[1], 
	      kd_tree[i].bounds.max[2], kd_tree[i].l_child, kd_tree[i].r_child,
	      kd_tree[i].parent);
  }

  // cleanup
  delete[] hists;

}
//----------------------------------------------------------------------------
//
// search the tree looking for a point
//
// pt: target point
// start_node: index of starting node of search (usually 0 (root), but
//  the caller may have more information and can shorten the search by
//  providing a node closer to the leaves)
//
// returns: leaf node index containing the target 
//  -1 if not found, indicates either an erroneous tree 
//  or a target point out of bounds of the entire domain at the root node
//
int Blocking::SearchTree(float *pt, int start_node) {

  int node = start_node;
  int i;

  while (1) {

    if (kd_tree[node].l_child == -1) // leaf; done
      return node;

    // check left child
    for (i = 0; i < dim; i++) {
      if (kd_tree[kd_tree[node].l_child].bounds.min[i] > pt[i] ||
	  kd_tree[kd_tree[node].l_child].bounds.max[i] < pt[i])
	break;
    }

    if (i == dim) {
      node = kd_tree[node].l_child;
      continue;
    }

    // check right child
    for (i = 0; i < dim; i++) {
      if (kd_tree[kd_tree[node].r_child].bounds.min[i] > pt[i] ||
	  kd_tree[kd_tree[node].r_child].bounds.max[i] < pt[i])
	break;
    }

    if (i == dim) {
      node = kd_tree[node].r_child;
      continue;
    }

    fprintf(stderr, "Error: TreeSearch() could not find target point\n");
    return -1;

  }

}
//----------------------------------------------------------------------------
//
// retrieves a tree leaf node
//
// index; leaf node index
// leaf: (output) leaf data
//
void Blocking::GetLeaf(int index, leaf_t *leaf) {

  leaf->gid = index + start_b;
  leaf->proc = kd_tree[index].proc;
  for (int i = 0; i < dim; i++) {
    leaf->bounds.min[i] = kd_tree[index].bounds.min[i];
    leaf->bounds.max[i] = kd_tree[index].bounds.max[i];
  }

}
//----------------------------------------------------------------------------
//
// finds index of target value of a sorted array using binary search
//
// start: starting index (eg., 0)
// num_vals: number of values starting at starting index
// vals: array of values
// target: target value
//
// returns: nearest index to target
//
int Blocking::BinarySearch(int start, int num_vals, int *vals, int target) {

    int lo = start;
    int hi = start + num_vals - 1;
    int mid = (lo + hi) / 2;
    while (hi - lo > 1 && mid > 0) {
      if (vals[lo] >= target) {
	mid = lo;
	break;
      }
      if (target >= vals[hi]) {
	mid = hi;
	break;
      }
      if (vals[mid] < target)
	lo = mid;
      else if (target < vals[mid])
	hi = mid;
      else
	break;
      mid = (lo + hi) / 2;
    }

    return mid;

}
//----------------------------------------------------------------------------
//
// Add new children to a parent in the tree
//
// parent: index of parent node
// split_dir: split direction (0 to dim - 1)
// split_frac: fraction of parent bounds where to split children (0.0 - 1.0)
//
void Blocking::AddChildren(int parent, int split_dir, float split_frac) {

  // debug
//   fprintf(stderr, "ready to add node parent %d split dir %d split_frac %.1f\n",
// 	  parent, split_dir, split_frac);

  kd_node_t node; // one kd tree node

  // for all children
  for (int child = 0; child < 2; child++) {

    node.parent = parent;
    node.proc = (int)kd_tree.size() % groupsize; // round robin for now
    node.l_child = -1; // currently a leaf node
    node.r_child = -1;

    for (int i = 0; i < dim; i++) {

      if (i == split_dir) { // the split direction
	// map split index to point in the bounds
	float split_point = kd_tree[node.parent].bounds.min[i] + split_frac *
	  (kd_tree[node.parent].bounds.max[i] - 
	   kd_tree[node.parent].bounds.min[i]);
	if (child == 0) { // left child
	  node.bounds.min[i] = kd_tree[node.parent].bounds.min[i];
	  node.bounds.max[i] = split_point;
	}
	else { // right child
	  node.bounds.min[i] = split_point;
	  node.bounds.max[i] = kd_tree[node.parent].bounds.max[i];
	}
      }
      else { // other directions unaffected
	node.bounds.min[i] = kd_tree[node.parent].bounds.min[i];
	node.bounds.max[i] = kd_tree[node.parent].bounds.max[i];
      }

    }

    // add the node and point parent to it
    kd_tree.push_back(node);
    if (child == 0)
      kd_tree[node.parent].l_child = (int)kd_tree.size() - 1;
    else
      kd_tree[node.parent].r_child = (int)kd_tree.size() - 1;

    // debug
//     fprintf(stderr, "adding node: proc %d min[%.1f %.1f %.1f] "
// 	    "max[%.1f %.1f %.1f] l_child %d r_child %d parent %d\n",
// 	    node.proc, node.bounds.min[0], 
// 	    node.bounds.min[1], node.bounds.min[2], 
// 	    node.bounds.max[0], node.bounds.max[1], 
// 	    node.bounds.max[2], node.l_child, node.r_child,
// 	    node.parent);

  } // children

}
//----------------------------------------------------------------------------
//
// callback function to compute a global histogram
//  by merging individual histograms
//
// items: pointers to input / output items, result in items[0]
// char * is used as a generic pointers to bytes, not necessarily to strings
// gids: gloabl ids of items to be reduced (unused)
// num_items: total number of input items
// hdr: quantity information
//
static void KdTree_MergeHistogram(char **items, int *gids, int num_items,
				  int *hdr) {

  gids = gids; // quiet compiler warning

  // todo: need to offset histograms for range if/when histograms are not the
  // global range, as they are now

  // add histograms
  for (int i = 1; i < num_items; i++) {
    for (int j = 1; j < hdr[0]; j++)
      ((int **)items)[0][j] += ((int **)items)[i][j];
  }

}
//----------------------------------------------------------------------------
//
// callback function to create a received item
//
// hdr: quantity information
//
// char * is used as a generic pointers to bytes, not necessarily to strings
//
// side effects: allocates the item
//
// returns: pointer to the item
//
static char *KdTree_CreateHistogram(int *hdr) {

  int *bins = new int[hdr[0]];
  return (char *)bins;

}
//----------------------------------------------------------------------------
//
// callback function to destroy a received item
//
// item: item to be destroyed
//
static void KdTree_DestroyHistogram(void *item) {

  delete[] (int *)item;

}
//----------------------------------------------------------------------------
//
// callback function to create a DIY datatype for received item being merged
//
// item: pointer to the item (nused)
// dtype: pointer to the datatype
// hdr: quantity information
//
// side effects: commits the datatype but DIY will cleanup datatype for you
//
static void KdTree_CreateHistogramType(void *item, DIY_Datatype *dtype,
					int *hdr) {

  item = item; // quiet compiler warning

  DIY_Create_vector_datatype(hdr[0], 1, DIY_INT, dtype);

}
//----------------------------------------------------------------------------
