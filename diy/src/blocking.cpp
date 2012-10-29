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
Blocking::Blocking(int dim, int tot_b, int64_t *data_size, bool share_face,
		   int *ghost, int64_t *given, 
		   Assignment *assignment, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->dim = dim;
  this->tot_b = tot_b;
  this->comm = comm;
  this->share_face = share_face;
  assign = assignment;
  nb = assign->NumBlks();

  for (int i = 0; i < dim; i++) {
    data_min[i] = 0.0; // for now, until data min and max are added to API
    data_max[i] = data_size[i] - 1.0; // for now
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
    blocks[i].gid = assign->Lid2Gid(i);
    blocks[i].proc = -1; // no need to store my own process id
  }

}
//----------------------------------------------------------------------------
//
// constructor for adopting an existing blocking
//
// dim: number of dimensions
// gids: global ids of my local blocks
// bounds: block bounds (extents) of local blocks
// assignment: pointer to asignment class
// comm: MPI communicator
//
Blocking::Blocking(int dim, int *gids, bb_t *bounds,
		   Assignment *assignment, MPI_Comm comm) {

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->dim = dim;
  tot_b = -1; // unknown
  this->comm = comm;
  assign = assignment;
  nb = assign->NumBlks();

  for (int i = 0; i < dim; i++)
    data_size[i] = -1; // unknown

  blocks = new gb_t[nb];

#pragma omp parallel for
  for (int i = 0; i < nb; i++) {
    blocks[i].gid = gids[i];
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
void Blocking::NumLatBlocks(int64_t *lat_nblocks) {

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
int64_t Blocking::BlockStartsSizes(int lid, int64_t *starts, int64_t *sizes) {

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
int64_t Blocking::BlockSizes(int lid, int64_t *sizes) {

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
void Blocking::BlockStarts(int lid, int64_t *starts) {

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
    tot_size *= (blocks[lid].bb.max[i] - blocks[lid].bb.min[i] + 1);

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

  // DEPRECATED
//   for (int i = 0; i < dim; i++) {
//     if (blocks[lid].bb.min[i] > data_min[i])
//       bounds->min[i] = blocks[lid].bb.min[i] + ghost[2 * i];
//     else
//       bounds->min[i] = blocks[lid].bb.min[i];
//     if (blocks[lid].bb.max[i] < data_max[i])
//       bounds->max[i] = blocks[lid].bb.max[i] - ghost[2 * i + 1];
//     else
//       bounds->max[i] = blocks[lid].bb.max[i];
//   }

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
void Blocking::ComputeBlocking(int64_t *given) {

  // factor data dimensions
  FactorDims(given);

  // volume bounds in row-major index order (x, y, z)
  // x changes fastest, z slowest
  int lid = 0;
  int gid = 0;
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
void Blocking::FactorDims(int64_t *given) { 

  int rem = tot_b; // unfactored remaining portion of tot_b
  int64_t block_dim[DIY_MAX_DIM]; // current block size
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
	fprintf(stderr,"Unable to block the volume with given[%d] = %Ld "
		"dimension. Please provide different 'given' constraints and "
		"rerun.\n", i, given[i]);
#else
	fprintf(stderr,"Unable to block the volume with given[%d] = %ld "
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
    block_size[i] = (int64_t)(roundf((float)data_size[i] / (float)lat_size[i]));

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
    Gid2Indices(assign->Lid2Gid(lid), mi, mj);
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
    Gid2Indices(assign->Lid2Gid(lid), mi, mj, mk);
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
    Gid2Indices(assign->Lid2Gid(lid), mi, mj, mk, ml);
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

  j = gid / lat_size[0];
  i = gid % lat_size[0];

}
//---------------------------------------------------------------------------
//
// gets lattice indices given a global block id (3D version)
//
// gid: global block id
// i, j, k: lattice indices (output)
//
void Blocking::Gid2Indices(int gid, int& i, int& j, int& k) {

  k = gid / (lat_size[0] * lat_size[1]) ; 
  j = (gid % (lat_size[0] * lat_size[1])) / lat_size[0]; 
  i = gid % lat_size[0]; 

}
//---------------------------------------------------------------------------
//
// gets lattice indices given a global block id (4D version)
//
// gid: global block id
// i, j, k, l: lattice indices (output)
//
void Blocking::Gid2Indices(int gid, int& i, int& j, int& k, int &l) {

  l = gid / (lat_size[0] * lat_size[1] * lat_size[2]); 
  int r = gid % (lat_size[0] * lat_size[1] * lat_size[2]); 
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

  return(j * lat_size[0] + i);

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

  return(k * lat_size[0] * lat_size[1] + j * lat_size[0] + i);

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

  return(l * lat_size[0] * lat_size[1] * lat_size[2] + k * lat_size[0] * 
	 lat_size[1] + j * lat_size[0] + i);

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

  int64_t starts[4]; // block starts

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
