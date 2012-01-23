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
// constructor
//
// dim: number of dimensions
// tot_b: total number of blocks
// data_size: data size in up to 4 dimensions
// share_face: whether neighboring blocks share a common face or are
//  separated by a gap of one unit
// ghost: ghost layer per side
// ghost_dir: where to apply ghost cells, for each dimension
// -1 = apply ghost layer to minimum sides of block only,
//  1 = apply ghost layer to maximum sides of block only,
//  0 = apply ghost layer equally to all sides of block
// ghost_dim: 
//  which dimensions to apply a ghost layer. a 1 means to apply a ghost layer
//  to that dimension, while a 0 means do not apply a ghost layer in that
//  dimension
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//   eg., {0, 0, 0, t} would result in t blocks in the 4th dimension
// assignment: pointer to asignment class
// comm: MPI communicator
//
Blocking::Blocking(int dim, int tot_b, int64_t *data_size, bool share_face,
		   int ghost, int* ghost_dir, int* ghost_dim, int64_t *given, 
		   Assignment *assignment, MPI_Comm comm) {

  Init(dim, tot_b, data_size, share_face, ghost, ghost_dir, ghost_dim, given,
       assignment, comm);
}

Blocking::Blocking(int dim, int tot_b, int64_t *data_size, bool share_face,
		   int ghost, int ghost_dir, int64_t *given, 
		   Assignment *assignment, MPI_Comm comm) {

  // assume ghost cells are applied to all dimensions.
  // assume use the same ghost_dir for all dimensions.
  int ghost_dim[4];
  int ghost_dir2[4];
  for(int i=0; i<dim; i++) {
    ghost_dim[i] = 1;
    ghost_dir2[i] = ghost_dir;
  }

  Init(dim, tot_b, data_size, share_face, ghost, ghost_dir2, ghost_dim, given,
       assignment, comm);

}

// initialize the object. common constructor code.
void Blocking::Init(int dim, int tot_b, int64_t *data_size, bool share_face,
		   int ghost, int* ghost_dir, int* ghost_dim, int64_t *given, 
		   Assignment *assignment, MPI_Comm comm) {
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  this->dim = dim;
  this->tot_b = tot_b;
  this->comm = comm;
  this->assign = assignment;
  nb = assign->NumBlks();
  this->time_starts = NULL;

  for (int i = 0; i < dim; i++)
    (this->data_size)[i] = data_size[i];

  bb_list = new bb_t[nb];
  rbb_list = new bb_t[tot_b];
  ComputeBlocking(share_face, ghost, ghost_dir, ghost_dim, given);
}

//----------------------------------------------------------------------------
//
// destructor
//
Blocking::~Blocking() {

  delete[] bb_list;
  delete[] rbb_list;
  if(time_starts != NULL)
  {
    delete[] time_starts;
  }

}
//----------------------------------------------------------------------------
//
// block starts and sizes
//
// lid: local block id
// starts: pointer to allocated array of starting block extents (output)
// sizes: pointer to allocated array of block sizes (output)
//
// returns: total size of block (product of sizes in each dimension)
//
int64_t Blocking::BlockStartsSizes(int lid, int64_t *starts, int64_t *sizes) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++) {
    starts[i] = bb_list[lid].min[i];
    sizes[i] = bb_list[lid].max[i] - starts[i] + 1;
    tot_size *= sizes[i];
  }

  return tot_size;

}
//----------------------------------------------------------------------------
//
// get the block bounds given the local block id
//
// lid: local block id
// from: the lower bounds of the block (output)
// to: the upper bounds of the block (output)
//
void Blocking::GetBlockBounds(int lid, int64_t* from, int64_t* to) {

  for(int i=0; i<dim; i++)
  {
    from[i] = bb_list[lid].min[i];
    to[i] = bb_list[lid].max[i];
  }
}
//----------------------------------------------------------------------------
//
// get the real block bounds given the local block id
// the real block bounds are the ones not counting ghost cells
//
// lid: local block id
// from: the lower bounds of the block (output)
// to: the upper bounds of the block (output)
//
void Blocking::GetRealBlockBounds(int lid, int64_t* from, int64_t* to) {

  int gid = assign->RoundRobin_lid2gid(lid);
  for(int i=0; i<dim; i++)
  {
    from[i] = rbb_list[gid].min[i];
    to[i] = rbb_list[gid].max[i];
  }
}
//--------------------------------------------------------------------------
//
// block sizes
//
// lid: local block id
// sizes: pointer to allocated array of block sizes (output)
//
// returns: total size of block (product of sizes in each dimension)
//
int64_t Blocking::BlockSizes(int lid, int64_t *sizes) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++) {
    sizes[i] = bb_list[lid].max[i] - bb_list[lid].min[i] + 1;
    tot_size *= sizes[i];
  }

  return tot_size;

}
//--------------------------------------------------------------------------
//
// block starts
//
// lid: local block id
// starts: pointer to allocated array of block starts (output)
//
void Blocking::BlockStarts(int lid, int64_t *starts) {

  for (int i = 0; i < dim; i++)
    starts[i] = bb_list[lid].min[i];

}
//--------------------------------------------------------------------------
//
// total block size
//
// lid: local block id
//
// returns: total size of block (product of sizes in each dimension)
//
int64_t Blocking::TotalBlockSize(int lid) {

  int64_t tot_size = 1;

  for (int i = 0; i < dim; i++)
    tot_size *= (bb_list[lid].max[i] - bb_list[lid].min[i] + 1);

  return tot_size;

}
//--------------------------------------------------------------------------
//
// calculates a default blocking of the data
// by factoring the total number of verices in alternating x, y, z, t directions
// supports 2D, 3D, and 4D depending on number of dimensions when
// constructing the object
// share_face: whether neighboring blocks share a common face or are
//  separated by a gap of one unit
// ghost: ghost layer per side
// ghost_dir: where to apply ghost cells, for each dimension
// -1 = apply ghost layer to minimum sides of block only,
//  1 = apply ghost layer to maximum sides of block only,
//  0 = apply ghost layer equally to all sides of block
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//   eg., {0, 0, 0, t} would result in t blocks in the 4th dimension
//
//
void Blocking::ComputeBlocking(bool share_face, int ghost, int* ghost_dir, 
			       int* ghost_dim, int64_t *given) {

  // factor data dimensions
  FactorDims(given);

  if(time_starts != NULL)
  {
    delete[] time_starts;
    time_starts = NULL;
  }
  time_starts = new int64_t[lat_size[3]];
  int time_index = 0;

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
  int d[4]; // actual current block size
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

	  // store all blocks in rbb_list
	  rbb_list[gid].min[0] = cur[0];
	  if (i == (int)lat_size[0] - 1) {
	    rbb_list[gid].max[0] = data_size[0] - 1;
	  }
	  else {
	    rbb_list[gid].max[0] = cur[0] + d[0] - gap;
	  }

	  // only store my local blocks in bb_list
	  if (assign->RoundRobin_gid2proc(gid) == rank) {
	    bb_list[lid].min[0] = cur[0];
	    if (i == (int)lat_size[0] - 1) {
	      bb_list[lid].max[0] = data_size[0] - 1;
	    }
	    else {
	      bb_list[lid].max[0] = cur[0] + d[0] - gap;
	    }
	  }

	  // y
	  if (lat_size[1] * block_size[1] > data_size[1] &&
	      (j + 1) * block_size[1] + 
	      (lat_size[1] - j - 1) * (block_size[1] - 1) >=
	      data_size[1])
	    d[1] = block_size[1] - 1;
	  else
	    d[1] = block_size[1];

	  rbb_list[gid].min[1] = cur[1];
	  if (j == (int)lat_size[1] - 1) {
	    rbb_list[gid].max[1] = data_size[1] - 1;
	  }
	  else {
	    rbb_list[gid].max[1] = cur[1] + d[1] - gap;
	  }

	  if (assign->RoundRobin_gid2proc(gid) == rank) {
	    bb_list[lid].min[1] = cur[1];
	    if (j == (int)lat_size[1] - 1) {
	      bb_list[lid].max[1] = data_size[1] - 1;
	    }
	    else {
	      bb_list[lid].max[1] = cur[1] + d[1] - gap;
	    }
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

	    rbb_list[gid].min[2] = cur[2];
	    if (k == (int)lat_size[2] - 1) {
	      rbb_list[gid].max[2] = data_size[2] - 1;
	    }
	    else {
	      rbb_list[gid].max[2] = cur[2] + d[2] - gap;
	    }

	    if (assign->RoundRobin_gid2proc(gid) == rank) {
	      bb_list[lid].min[2] = cur[2];
	      if (k == (int)lat_size[2] - 1) {
		bb_list[lid].max[2] = data_size[2] - 1;
	      }
	      else {
		bb_list[lid].max[2] = cur[2] + d[2] - gap;
	      }
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

	    rbb_list[gid].min[3] = cur[3];
	    if (l == (int)lat_size[3] - 1) {
	      rbb_list[gid].max[3] = data_size[3] - 1;
	    }
	    else {
	      rbb_list[gid].max[3] = cur[3] + d[3] - gap;
	    }
	    if (data_size[3] == 1) { // special case for 4D but static
	      rbb_list[gid].max[3] = 0;
	    }

	    if (assign->RoundRobin_gid2proc(gid) == rank) {
	      bb_list[lid].min[3] = cur[3];
	      if (l == (int)lat_size[3] - 1) {
		bb_list[lid].max[3] = data_size[3] - 1;
	      }
	      else {
		bb_list[lid].max[3] = cur[3] + d[3] - gap;
	      }
	      if (data_size[3] == 1) { // special case for 4D but static
		bb_list[lid].max[3] = 0;
	      }
	    }
	    time_starts[time_index] = rbb_list[gid].min[3];
	  }

	  if (assign->RoundRobin_gid2proc(gid) == rank)
	    lid++;
	  gid++;
	  cur[0] += d[0];
	}
	cur[1] += d[1];
      }
      cur[2] += d[2];
    }
    cur[3] += d[3];
    time_index++;
  }

  // ghost cells
  if (ghost > 0)
    ApplyGhost(ghost, ghost_dir, ghost_dim);

  // debug: print the bb and rbb list
#if 0
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if(rank == 0)
   {
     for (int i = 0; i < tot_b; i++) {
       gid = i;
       fprintf(stderr, "proc[%i] rbb_list[%d] gid = %d  "
	       "min = [%lld %lld %lld %lld] max = [%lld %lld %lld %lld]\n",
	       rank, i, gid, rbb_list[i].min[0], rbb_list[i].min[1],
	       rbb_list[i].min[2], rbb_list[i].min[3], rbb_list[i].max[0],
	       rbb_list[i].max[1], rbb_list[i].max[2], rbb_list[i].max[3]);
     }
   }
   for (int i = 0; i < nb; i++) {
     gid = assign->RoundRobin_lid2gid(i);
     fprintf(stderr, "proc[%i] bb_list[%d] gid = %d  "
	     "min = [%lld %lld %lld %lld] max = [%lld %lld %lld %lld]\n", rank,
 	     i, gid, bb_list[i].min[0], bb_list[i].min[1], bb_list[i].min[2],
 	     bb_list[i].min[3], bb_list[i].max[0], bb_list[i].max[1], 
 	     bb_list[i].max[2], bb_list[i].max[3]);
   }
#endif

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
  int64_t block_dim[MAX_DIM]; // current block size
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
	fprintf(stderr,"Unable to block the volume with given[%d] = %ld "
		"dimension. Please provide different 'given' constraints and "
		"rerun.\n", i, given[i]);
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

}
//---------------------------------------------------------------------------
//
// adds the ghost layer to the block bounds
//
// ghost: ghost layer per side
// ghost_dir: where to apply ghost cells, one int per dimension
// -1 = apply ghost layer to minimum sides of block only,
//  1 = apply ghost layer to maximum sides of block only,
//  0 = apply ghost layer equally to all sides of block
// ghost_dim: 
//  which dimensions to apply a ghost layer. a 1 means to apply a ghost layer
//  to that dimension, while a 0 means do not apply a ghost layer in that
//  dimension
//
void Blocking::ApplyGhost(int ghost, int* ghost_dir, int* ghost_dim) {

  for (int i = 0; i < nb; i++)  // for each block
  {
    for(int j=0; j<4; j++)      // for each dimension
    {
      if(dim > j)
      {
	if (ghost_dim[j] == 1)
	{

	  if (ghost_dir[j] < 1)
	    bb_list[i].min[j] = max(bb_list[i].min[j] - ghost, (int64_t)0);

	  if (ghost_dir[j] > -1)
	    bb_list[i].max[j] = min(bb_list[i].max[j] + ghost, data_size[j]-1);
	}
      }
    }
  }
}
//---------------------------------------------------------------------------
//
// gets all neighbor blocks of a local block
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
    Gid2Indices(assign->RoundRobin_lid2gid(lid), mi, mj);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	gid = Indices2Gid(mi + i, mj + j);
	if (gid >= 0) {
	  neigh.gid = gid;
	  neigh.proc = assign->RoundRobin_gid2proc(gid);
	  neighbors.push_back(neigh);
	}
      }
    }
    break;
  case 3:
    Gid2Indices(assign->RoundRobin_lid2gid(lid), mi, mj, mk);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	for (k = -1; k <= 1; k++) {
	  gid = Indices2Gid(mi + i, mj + j, mk + k);
	  if (gid >= 0) {
	    neigh.gid = gid;
	    neigh.proc = assign->RoundRobin_gid2proc(gid);
	    neighbors.push_back(neigh);
	  }
	}
      }
    }
    break;
  case 4:
    Gid2Indices(assign->RoundRobin_lid2gid(lid), mi, mj, mk, ml);
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	for (k = -1; k <= 1; k++) {
	  for (l = -1; l <= 1; l++) {
	    gid = Indices2Gid(mi + i, mj + j, mk + k, ml + l);
	    if (gid >= 0) {
	      neigh.gid = gid;
	      neigh.proc = assign->RoundRobin_gid2proc(gid);
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
// returns global id of (first) neighboring block found containing the point
//
// gid: global id of current block
// pt: point (up to number of dimensions specified when creating the class
// ghost: ghost layer per side
// ghost_dir: 
// -1 = apply ghost layer to minimum sides of block only,
//  1 = apply ghost layer to maximum sides of block only,
//  0 = apply ghost layer equally to all sides of block
//
// returns: neighboring block global id or -1 if out of the domain
//
int Blocking::Pt2NeighGid(int gid, float *pt, int ghost, int ghost_dir) {

  int i, j, k, l; // indices of current block, intialized to quiet compiler
  i = j = k = l = 0;

  if (dim > 3)
    Gid2Indices(gid, i, j, k, l); 
  else if (dim > 2)
    Gid2Indices(gid, i, j, k); 
  else
    Gid2Indices(gid, i, j); 

  int l0, l1, k0, k1;
  l0 = l1 = 0;
  k0 = k1 = 0;
  if (dim > 3) {
    l0 = -1;
    l1 = 1;
  }
  if (dim > 2) {
    k0 = -1;
    k1 = 1;
  }

  for (int nl = l0; nl <= l1; nl++) {
    for (int nk = k0; nk <= k1; nk++) { 
      for (int nj = -1; nj <= 1; nj++)  {
	for (int ni = -1; ni <= 1; ni++)  {
	  int bi[4] = {i + ni, j + nj, k + nk, l + nl};
	  if (dim > 3 && (nl + l < 0 || nl + l > (int)lat_size[3] - 1))
	    continue; 
	  else if (dim > 2 && (nk + k < 0 || nk + k > (int)lat_size[2] - 1)) 
	    continue;
	  else if (dim > 1 && (nj + j < 0 || nj + j > (int)lat_size[1] - 1))
	    continue;
	  else if (ni + i < 0 || ni + i > (int)lat_size[0] - 1) 
	    continue;
	  else if (IsIn(pt, bi, ghost, ghost_dir)) {
	    if (dim > 3)
	      return(Indices2Gid(ni + i, nj + j, nk + k, nl + l));
	    else if (dim > 2)
	      return(Indices2Gid(ni + i, nj + j, nk + k));
	    else
	      return(Indices2Gid(ni + i, nj + j));
	  }
	}
      }
    }
  }

  return(-1); 

}
//--------------------------------------------------------------------------
//
// checks if the point pt is in the block indexed by bi
//
// pt: point (x, y, z, ...)
// bi: block index (i, j, k, ...)
// ghost: ghost layer per side
// ghost_dir: 
// -1 = apply ghost layer to minimum sides of block only,
//  1 = apply ghost layer to maximum sides of block only,
//  0 = apply ghost layer equally to all sides of block
//
// returns: true / false
//
bool Blocking::IsIn(float *pt, int *bi, int ghost, int ghost_dir) {

  int lghost, rghost; // left and right side ghost (min side and max side)
  lghost = rghost = ghost;
  if (ghost_dir == -1)
    rghost = 0;
  if (ghost_dir == 1)
    lghost = 0;

  // find the global id of the block being tested
  int gid;
  if (dim > 3)
    gid = Indices2Gid(bi[0], bi[1], bi[2], bi[3]);
  else if (dim > 2)
    gid = Indices2Gid(bi[0], bi[1], bi[2]);
  else
    gid = Indices2Gid(bi[0], bi[1]);

  // the bounds of the block (not counting ghost cells)
  bb_t* bounds = &(rbb_list[gid]);

  for (int i = 0; i < dim; i++) {
    if (pt[i] < 0 || pt[i] > data_size[i]-1)
      return false;
    if (bi[i] < 0 || bi[i] >= (int)lat_size[i])
      return false;
    if(pt[i] < bounds->min[i] || pt[i] > bounds->max[i])
      return false;
  }

  return(true); 

}
//--------------------------------------------------------------------------
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
  if (tsize == 1 || tb == 1 || starts[3] == time_starts[g])
    return true;

  return false;

}
//-----------------------------------------------------------------------
