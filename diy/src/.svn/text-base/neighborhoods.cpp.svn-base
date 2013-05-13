//---------------------------------------------------------------------------
//
// neighborhoods class
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
//--------------------------------------------------------------------------

#include "neighborhoods.hpp"

//--------------------------------------------------------------------------
//
// constructs and initializes my neighborhoods
// using round robin assignment and implicit computation of neighborhoods
//
// did: domain id
// blocking: pointer to blocking class
// assignment: pointer to round robin assignment class
// comm: MPI commnicator
// wrap: whether wraparound neighbors are used (ignored for now, placeholder
//  for later development)
// nhdr: optional number of header counts
//
Neighborhoods::Neighborhoods(int did, Blocking *blocking, 
			     Assignment *assignment, 
			     MPI_Comm comm, bool wrap, int nhdr) {

  wrap = wrap; // quiet compiler warning

  this->did = did;
  this->comm = comm;
  this->blocking = blocking;
  this->assign = assignment;
  this->nhdr = nhdr;
  rem_ids = NULL; // unused in this version of the constructor
  num_rem_ids = NULL;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // init my blocks list
  for (int i = 0; i < assign->NumBlks(); i++)  {
    bl_t block;
    vector<struct gb_t> gbs;
    block.gid = blocking->Lid2Gid(i);
    blocking->GetNeighbors(i, gbs);
    for (vector<struct gb_t>::iterator gi = gbs.begin(); gi != gbs.end(); 
	 gi++) {
      ne_t neighbor;
      neighbor.gb = *gi;
      neighbor.gb.neigh_dir = 0x00;
      neighbor.wrap_dir = 0x00;
      block.neighbors.push_back(neighbor);
    }
    blocks.push_back(block);
  }

  discovery = false;

  // discover neighbor block bounds with a neighbor exchange
  tag = 0; // set the tag for the neighbor bounds exchdange
  GetNeighborBounds();
  tag = 0; // reset the tag for general use

}
//--------------------------------------------------------------------------
//
// constructs and initializes my neighborhoods
// using process order assignment and explicit listing of neighborhoods
// and optional neighbor discovery
//
// did: domain id
// blocking: pointer to blocking class
// assignment: pointer to process order assignment class
// rem_ids: remote ids used for neighbor discovery
// num_rem_ids: number of remote ids for each local block
// vids: local vertex ids for each local block that needs neighbor discovery
// num_vids: number of vids for each local block
// neighbors: neighbor lists for each of my local blocks, in lid order
//  lists can be partial, containing only known information, as long as 
//  the rem_data lists contain enough information for DIY to discover the rest
//  of the unknown neighbors. Unknown (remote) block bounds can be 
//  uninitialized.
// num_neighbors: number of neighbors known so far for each of my 
//  local blocks, in lid order
// comm: MPI commnicator
// wrap: whether wraparound neighbors are used
// nhdr: optional number of header counts
//
Neighborhoods::Neighborhoods(int did, Blocking *blocking, 
			     Assignment *assignment, 
			     ri_t **rem_ids, int *num_rem_ids, int **vids, 
			     int *num_vids, gb_t **neighbors, 
			     int *num_neighbors, MPI_Comm comm, bool wrap,
			     int nhdr) {

  this->did = did;
  this->comm = comm;
  this->blocking = blocking;
  this->assign = assignment;
  this->nhdr = nhdr;
  this->rem_ids = rem_ids;
  this->num_rem_ids = num_rem_ids;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // init my blocks list
  for (int i = 0; i < assign->NumBlks(); i++)  {
    bl_t block;
    block.gid = blocking->Lid2Gid(i);
    for (int j = 0; j < num_neighbors[i]; j++) {
      ne_t neighbor;
      neighbor.gb.gid = neighbors[i][j].gid;
      neighbor.gb.proc = neighbors[i][j].proc;
      for (int k = 0; k < blocking->GetDim(); k++) {
	neighbor.gb.bb.min[k] = neighbors[i][j].bb.min[k];
	neighbor.gb.bb.max[k] = neighbors[i][j].bb.max[k];
      }
      neighbor.wrap_dir = 0x00; // initally
      if (wrap)
	neighbor.gb.neigh_dir = neighbors[i][j].neigh_dir;
      else
	neighbor.gb.neigh_dir = 0x00; 
      block.neighbors.push_back(neighbor);
    }
    blocks.push_back(block);
  }

  // ensure that there is an entry in the neighbors for each process in the
  // remote data, so there can be an exchange of information
  // add dummy neighbor blocks as needed

  if (rem_ids && num_rem_ids && vids && num_vids) { // neighbor discovery

    discovery = true;

    // for all local blocks
    for (int i = 0; i < assign->NumBlks(); i++) {

      // for all remote data pertaining to current local block
      for (int k = 0; k < num_rem_ids[i]; k++) {

	// for all neighbors of current local block
	int j;
	for (j = 0; j < (int)blocks[i].neighbors.size(); j++) {
	  if (rem_ids[i][k].proc == blocks[i].neighbors[j].gb.proc) {
	    break;
	  }
	}
	if (j == (int)blocks[i].neighbors.size()) {
	  ne_t neighbor;
	  neighbor.gb.gid = -1;
	  neighbor.gb.proc = rem_ids[i][k].proc;
	  blocks[i].neighbors.push_back(neighbor);
	}

      }

    }

  }
  else
    discovery = false;

  // discover neighbor blocks with a neighbor exchange
  tag = 0; // set the tag for the neighbor bounds exchange
  GetNeighbors(vids, num_vids);
  tag = 0; // reset the tag for general use

  // adjust bounds of wrapped neighbors w.r.t. local block
  if (wrap) 
    WrapNeighbors();

  // debug
  // DEPRECATED: specific tests for moab example
//   if (rank == 0) { // print only rank 0 to reduce output
//     for (int i = 0; i < assign->NumBlks(); i++) {
//       bb_t bb;
//       blocking->BlockBounds(i, &bb);
//       fprintf(stderr, "gid = %d mins [%.1f %.1f %.1f] maxs [%.1f %.1f %.1f] "
// 	      "neighbors = : ", blocks[i].gid,
// 	      bb.min[0],
// 	      bb.min[1],
// 	      bb.min[2],
// 	      bb.max[0],
// 	      bb.max[1],
// 	      bb.max[2]);

//       int n = 0; // number of remote neighbors
//       for (int j = 0; j < blocks[i].neighbors.size(); j++) {
// 	// print only remote neighbors
// 	if (rank != blocks[i].neighbors[j].gb.proc) {
// 	  fprintf(stderr, "neigh_gid %d neigh_proc %d "
// 		  "neigh_mins %.1f %.1f %.1f neigh_maxs %.1f %.1f %.1f : ", 
// 		  blocks[i].neighbors[j].gb.gid, blocks[i].neighbors[j].gb.proc,
// 		  blocks[i].neighbors[j].gb.bb.min[0],
// 		  blocks[i].neighbors[j].gb.bb.min[1],
// 		  blocks[i].neighbors[j].gb.bb.min[2],
// 		  blocks[i].neighbors[j].gb.bb.max[0],
// 		  blocks[i].neighbors[j].gb.bb.max[1],
// 		  blocks[i].neighbors[j].gb.bb.max[2]);
// 	  assert(blocks[i].neighbors[j].gb.gid > 256);
// 	  assert(blocks[i].neighbors[j].gb.bb.min[0] >= 
// 		 bb.min[0] - 0.5);
// 	  assert(blocks[i].neighbors[j].gb.bb.min[1] >= 
// 		 bb.min[1] - 0.5);
// 	  assert(blocks[i].neighbors[j].gb.bb.min[2] >= 
// 		 bb.min[2] - 0.5);
// 	  assert(blocks[i].neighbors[j].gb.bb.max[0] <= 
// 		 bb.max[0] + 0.5);
// 	  assert(blocks[i].neighbors[j].gb.bb.max[1] <= 
// 		 bb.max[1] + 0.5);
// 	  assert(blocks[i].neighbors[j].gb.bb.max[2] <= 
// 		 bb.max[2] + 0.5);
// 	  n++;
// 	}
//       }
//       fprintf(stderr, " number of remote neighbors = %d\n\n", n);
//       assert(n <= 9);
//     }
//   }

}
//--------------------------------------------------------------------------
//
// destructor
//
Neighborhoods::~Neighborhoods() {

  for (vector<bl_t>::iterator bi = blocks.begin(); bi != blocks.end(); bi++)
    bi->neighbors.clear();
  blocks.clear();

}
//---------------------------------------------------------------------------
//
// handles wrapped neighbors
//
// neighbors: neighbor lists for each of my local blocks, in lid order
//  lists can be partial, containing only known information, as long as 
//  the rem_data lists contain enough information for DIY to discover the rest
//  of the unknown neighbors. Unknown (remote) block bounds can be 
//  uninitialized.
// num_neighbors: number of neighbors known so far for each of my 
//  local blocks, in lid order
//
void Neighborhoods::WrapNeighbors() {

  int dim = blocking->GetDim();
  float pt[DIY_MAX_DIM]; // temporary point
  float my_size[DIY_MAX_DIM]; // size of local block
  float neigh_size[DIY_MAX_DIM]; // size of neighbor block

  // for all blocks
  for (int i = 0; i < assign->NumBlks(); i++) {

    // for all neighbors of current block
    for (int j = 0; j < (int)blocks[i].neighbors.size(); j++) {

      ne_t *neigh = &(blocks[i].neighbors[j]);
      bb_t bb;
      blocking->BlockBounds(i, &bb);
      for (int d = 0; d < dim; d++) {
	my_size[d] = bb.max[d] - bb.min[d];
	neigh_size[d] = neigh->gb.bb.max[d] - neigh->gb.bb.min[d];
      }

      // for each dimension,
      // construct theoretical point in center of neighbor, according to
      // neighbor direction, and if point dos not lie in neighbor bounds,
      // turn on wrapping bit for that dimension
      // note that neighbor directions are w.r.t. the sender (neighbor, not me)

      // x dimension
      if ((neigh->gb.neigh_dir & DIY_X1) == DIY_X1) {
	pt[0] = bb.min[0] - 0.5f * my_size[0];
	// todo: add ghost?

	// actual neighbor to my max side, but I see him to my min
	if (pt[0] < neigh->gb.bb.min[0] || pt[0] > neigh->gb.bb.max[0]) {
	  neigh->wrap_dir |= DIY_X0;
	  neigh->gb.bb.min[0] = bb.min[0] - neigh_size[0];
	  neigh->gb.bb.max[0] = bb.min[0];
	}
      }
      // actual neighbor is to my min side, but I see him to my max
      else if ((neigh->gb.neigh_dir & DIY_X0) == DIY_X0) {
	pt[0] = bb.max[0] + 0.5f * my_size[0];
	if (pt[0] < neigh->gb.bb.min[0] || pt[0] > neigh->gb.bb.max[0]) {
	  neigh->wrap_dir |= DIY_X1;
	  neigh->gb.bb.max[0] = bb.max[0] + neigh_size[0];
	  neigh->gb.bb.min[0] = bb.max[0];
	}
      }

      // y dimension
      if ((neigh->gb.neigh_dir & DIY_Y1) == DIY_Y1) {
	pt[1] = bb.min[1] - 0.5f * my_size[1];
	if (pt[1] < neigh->gb.bb.min[1] || pt[1] > neigh->gb.bb.max[1]) {
	  neigh->wrap_dir |= DIY_Y0;
	  neigh->gb.bb.min[1] = bb.min[1] - neigh_size[1];
	  neigh->gb.bb.max[1] = bb.min[1];
	}
      }
      else if ((neigh->gb.neigh_dir & DIY_Y0) == DIY_Y0) {
	pt[1] = bb.max[1] + 0.5f * my_size[1];
	if (pt[1] < neigh->gb.bb.min[1] || pt[1] > neigh->gb.bb.max[1]) {
	  neigh->wrap_dir |= DIY_Y1;
	  neigh->gb.bb.max[1] = bb.max[1] + neigh_size[1];
	  neigh->gb.bb.min[1] = bb.max[1];
	}
      }

      // z dimension
      if (dim > 2) { 
	if ((neigh->gb.neigh_dir & DIY_Z1) == DIY_Z1) {
	  pt[2] = bb.min[2] - 0.5f * my_size[2];
	  if (pt[2] < neigh->gb.bb.min[2] || pt[2] > neigh->gb.bb.max[2]) {
	    neigh->wrap_dir |= DIY_Z0;
	    neigh->gb.bb.min[2] = bb.min[2] - neigh_size[2];
	    neigh->gb.bb.max[2] = bb.min[2];
	  }
	}
	else if ((neigh->gb.neigh_dir & DIY_Z0) == DIY_Z0) {
	  pt[2] = bb.max[2] + 0.5f * my_size[2];
	  if (pt[2] < neigh->gb.bb.min[2] || pt[2] > neigh->gb.bb.max[2]) {
	    neigh->wrap_dir |= DIY_Z1;
	    neigh->gb.bb.max[2] = bb.max[2] + neigh_size[2];
	    neigh->gb.bb.min[2] = bb.max[2];
	  }
	}
      }

      // t dimension
      if (dim > 3) {
	if ((neigh->gb.neigh_dir & DIY_T1) == DIY_T1) {
	  pt[3] = bb.min[3] - 0.5f * my_size[3];
	  if (pt[3] < neigh->gb.bb.min[3] || pt[3] > neigh->gb.bb.max[3]) {
	    neigh->wrap_dir |= DIY_T0;
	    neigh->gb.bb.min[3] = bb.min[3] - neigh_size[3];
	    neigh->gb.bb.max[3] = bb.min[3];
	  }
	}
	else if ((neigh->gb.neigh_dir & DIY_T0) == DIY_T0) {
	  pt[3] = bb.max[3] + 0.5f * my_size[3];
	  if (pt[3] < neigh->gb.bb.min[3] || pt[3] > neigh->gb.bb.max[3]) {
	    neigh->wrap_dir |= DIY_T1;
	    neigh->gb.bb.max[3] = bb.max[3] + neigh_size[3];
	    neigh->gb.bb.min[3] = bb.max[3];
	  }
	}
      }      

    } // for all neighbors

  } // for all blocks

}
//--------------------------------------------------------------------------
//
// enqueues an item for sending to a neighbor
// will send to self block if neigh_gid equals gid of self
//
// lid: local id of my block
// item: item to be enqueued (char * pointer can point to anything, does not
// need to be chars)
// size: size of item in bytes
// neigh_gid: global id of neighboring destination block
// hdr: if nhdr > 0, pointer tp header counts
// (additional quantities of subitems within the item)
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
void Neighborhoods::EnqueueItem(int lid, char *item, size_t size, 
				int neigh_gid, int *hdr,
				void (*TransformItem)(char *, 
						      unsigned char)) {

  // todo: provide a non-copy mode for very large items

  int n; // local neighbor number of neigh_gid

  // copy item so that it persists
  char *p = new char[size];
  memcpy(p, item, size);

  // find the neighbor gid in my neighboring blocks
  for (n = 0; n < (int)blocks[lid].neighbors.size(); n++) {
    if (neigh_gid == blocks[lid].neighbors[n].gb.gid)
      break;
  }
  assert(n < (int)blocks[lid].neighbors.size()); // sanity

  // transform item if wrapping
  if (TransformItem && blocks[lid].neighbors[n].wrap_dir)
    TransformItem(p, blocks[lid].neighbors[n].wrap_dir);

  // enqueue payload
  blocks[lid].neighbors[n].items.push_back(p);

  // enqueue header
  int j = (int)(blocks[lid].neighbors[n].items.size());
  if ((int)blocks[lid].neighbors[n].hdr.size() < j)
    blocks[lid].neighbors[n].hdr.resize(j);
  for (int i = 0; i < nhdr; i++) 
    blocks[lid].neighbors[n].hdr[j].push_back(hdr[i]);

}
//------------------------------------------------------------------------
//
// enqueues an item for sending to one neighbor identified by a direction
// no direction is defined for sending to self block
//
// lid: local id of my block
// item: item to be enqueued (char * pointer can point to anything, does not
// need to be chars)
// size: size of item in bytes
// hdr: if nhdr > 0, pointer tp header counts
// (additional quantities of subitems within the item)
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
// dir: neighbor direction
//
void Neighborhoods::EnqueueItemDir(int lid, char *item, size_t size, int *hdr, 
				   void (*TransformItem)
				   (char *, unsigned char),
				   unsigned char neigh_dir) {

  // todo: provide a non-copy mode for very large items

  // find the neighbor gid in my neighboring blocks
  for (int n = 0; n < (int)blocks[lid].neighbors.size(); n++) {

    bb_t bb;
    blocking->BlockBounds(lid, &bb);
    int num_dims = blocking->GetDim();
    unsigned char test_mask = 0x03;
    bool reject = false; // reject this neighbor

    for (int d = 0; d < num_dims; d++) {
      unsigned char dir = neigh_dir & test_mask;
      for (int i = 0; i < d; i++)
	dir >>= 2;
      if ((dir == 0x00 &&
	   (bb.min[d] != blocks[lid].neighbors[n].gb.bb.min[d] ||
	    bb.max[d] != blocks[lid].neighbors[n].gb.bb.max[d])) ||
	  (dir == 0x01 && bb.min[d] < blocks[lid].neighbors[n].gb.bb.max[d]) ||
	  (dir == 0x02 && bb.max[d] > blocks[lid].neighbors[n].gb.bb.min[d])) {
	reject = true;
	break;
      }
      test_mask <<= 2;
    }

    if (reject) // go on to the next neighbor
      continue;

    // copy item so that it persists
    char *p = new char[size];
    memcpy(p, item, size);

    // transform item if wrapping
    if (TransformItem && blocks[lid].neighbors[n].wrap_dir)
      TransformItem(p, blocks[lid].neighbors[n].wrap_dir);

    // enqueue payload
    blocks[lid].neighbors[n].items.push_back(p);

    // enqueue header
    int j = (int)(blocks[lid].neighbors[n].items.size());
    if ((int)blocks[lid].neighbors[n].hdr.size() < j)
      blocks[lid].neighbors[n].hdr.resize(j);
    for (int i = 0; i < nhdr; i++) 
      blocks[lid].neighbors[n].hdr[j].push_back(hdr[i]);
  }

}
//------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors
// sending to self block dependings on a parameter, skips by default
//
// lid: local id of my block
// item: item to be enqueued (char * pointer can point to anything, does not
// need to be chars)
// size: size of item in bytes
// hdr: if nhdr > 0, pointer tp header counts
// (additional quantities of subitems within the item)
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
// symmetrical: enqueue in both directions in each axis
//  false = enqueue only to minimum neighbor (eg. left) in each axis
// self: whether to send to self
//
void Neighborhoods::EnqueueItemAll(int lid, char *item, size_t size, 
				   int *hdr, 
				   void (*TransformItem)(char *, 
							 unsigned char),
				   bool symmetrical, bool self) {

  // todo: provide a non-copy mode for very large items

  // find the neighbor gid in my neighboring blocks
  for (int n = 0; n < (int)blocks[lid].neighbors.size(); n++) {

    // skip self
    if (!self && blocks[lid].neighbors[n].gb.gid == DIY_Gid(did, lid))
      continue;

    bool symmetry_ok = true;

    if (!symmetrical) {
      bb_t bb;
      blocking->BlockBounds(lid, &bb);
      for (int d = 0; d < blocking->GetDim(); d++) {
	if (bb.max[d] <= blocks[lid].neighbors[n].gb.bb.min[d]) {
	  symmetry_ok = false;
	  break;
	}
      }
    }

    if (!symmetry_ok)
      continue;

    // copy item so that it persists
    char *p = new char[size];
    memcpy(p, item, size);

    // transform item if wrapping
    if (TransformItem && blocks[lid].neighbors[n].wrap_dir)
      TransformItem(p, blocks[lid].neighbors[n].wrap_dir);

    // enqueue payload
    blocks[lid].neighbors[n].items.push_back(p);

    // enqueue header
    int j = (int)(blocks[lid].neighbors[n].items.size());
    if ((int)blocks[lid].neighbors[n].hdr.size() < j)
      blocks[lid].neighbors[n].hdr.resize(j);
    for (int i = 0; i < nhdr; i++) 
      blocks[lid].neighbors[n].hdr[j].push_back(hdr[i]);

  }

}
//------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors near enough to receive it
// skips sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// size: size of item in bytes
// pt: point near the destination block
// dist: threshold distance for blocks from pt:
//   blocks less than or equal to dist of the pt will be destinations
// hdr: pointer to header (or NULL)
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
// symmetrical: (default true) enqueue in both directions in each axis
//  false = enqueue only to minimum neighbor (eg. left) in each axis
//
void Neighborhoods::EnqueueItemAllNear(int lid, char *item, size_t size,
				       float *pt, float dist, int *hdr,
				       void (*TransformItem)(char *, 
							     unsigned char),
				       bool symmetrical) {

  int dim = blocking->GetDim();
  float dir[DIY_MAX_DIM]; // offset direction
  float new_pt[DIY_MAX_DIM]; // offset point
  bb_t bb; // current block bounds
  int d; // current dimension

  DIY_Block_bounds(did, lid, &bb);

  // debug
//   fprintf(stderr, "checking site %.2lf %.2lf %.2lf\n", pt[0], pt[1], pt[2]);

  // todo: provide a non-copy mode for very large items

  // for all neighbors of lid
  for (int n = 0; n < (int)blocks[lid].neighbors.size(); n++) {

    // compute normalized vector from my block to the neighbor block
    // based on difference between mins
    int sum_d = 0; // number of nonzero direction components
    for (d = 0; d < dim; d++) {
      dir[d] = blocks[lid].neighbors[n].gb.bb.min[d] - bb.min[d];
      if (dir[d] > 0.0) {
	dir[d] = 1.0f;
	sum_d++;
      }
      else if (dir[d] < 0.0) {
	dir[d] = -1.0f;
	sum_d++;
      }
    }
    if (sum_d == 0) // self
      continue; // next neighbor

    // new point is offset from old point by dist in direction of vector
    for(d = 0; d < dim; d++)
      new_pt[d] = pt[d] + dir[d] * dist;

    // check if neighbor is near enough
    for (d = 0; d < dim; d++) {
      // if shifted point did not move into or past the neighbor, break and
      // proceed to next neighbor
      // note dist can be large enough to shift the point beyond the neighbor
      // that means the point was definitely near enough to neighbor
      if ((pt[d] < blocks[lid].neighbors[n].gb.bb.min[d] &&
	   new_pt[d] < blocks[lid].neighbors[n].gb.bb.min[d]) ||  
	  (pt[d] > blocks[lid].neighbors[n].gb.bb.max[d] &&
	   new_pt[d] > blocks[lid].neighbors[n].gb.bb.max[d]))
	break;
    }
    if (d < dim)
      continue; // next neighbor

    // check if neighbor passes symmetry test
    if (!symmetrical) {

      bool pass = false;

      if (d == 2) {
	if ((dir[0] <= 0 && dir[1] <= 0) ||
	    (dir[0] == 1 && dir[1] == -1))
	  pass = true;
      }

      if (d == 3) {
	if ((dir[0] <= 0 && dir[1] <= 0 && dir[2] <= 0) ||
	    (dir[0] == 1 && dir[1] == -1 && dir[2] <= 0) ||
	    (dir[1] == 1 && dir[2] == -1 && dir[0] <= 0) ||
	    (dir[2] == 1 && dir[0] == -1 && dir[1] <= 0))
	  pass = true;
      }

      //todo:  need to check if the 4D case is right
      if (d == 4) {
	if ((dir[0] <= 0 && dir[1] <= 0 && dir[2] <= 0  && dir[3] <= 0) ||
	    (dir[0] == 1 && dir[1] == -1 && dir[2] <= 0 && dir[3] <= 0) ||
	    (dir[1] == 1 && dir[2] == -1 && dir[3] <= 0 && dir[0] <= 0) ||
	    (dir[2] == 1 && dir[3] == -1 && dir[0] <= 0 && dir[1] <= 0) ||
	    (dir[3] == 1 && dir[0] == -1 && dir[1] <= 0 && dir[2] <= 0))
	  pass = true;
      }

      if (!pass)
	continue; // next neighbor

    }

    // debug
//     fprintf(stderr, "lid %d site %.2lf %.2lf %.2lf to gid %d "
// 	    "in dir %.2lf %.2lf %.2lf\n",
// 	    lid, pt[0], pt[1], pt[2], 
// 	    blocks[lid].neighbors[n].gb.gid,
// 	    dir[0], dir[1], dir[2]);

    // copy item so that it persists
    char *p = new char[size];
    memcpy(p, item, size);

    // transform item if wrapping
    if (TransformItem && blocks[lid].neighbors[n].wrap_dir)
      TransformItem(p, blocks[lid].neighbors[n].wrap_dir);

    // enqueue payload
    blocks[lid].neighbors[n].items.push_back(p);

    // enqueue header
    int j = (int)blocks[lid].neighbors[n].items.size();
    if ((int)blocks[lid].neighbors[n].hdr.size() < j)
      blocks[lid].neighbors[n].hdr.resize(j);
    for (int i = 0; i < nhdr; i++) 
      blocks[lid].neighbors[n].hdr[j].push_back(hdr[i]);


  } // for all neighbors

}
//------------------------------------------------------------------------
//
// exchanges items with all neighbors
//
// items: received items for each of my blocks [lid][outoput] (output)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// suggested value: 0.5-0.75
// ItemDtype: pointer to user-supplied function that creates an MPI datatype 
//  for an item to be sent or received
// discovery (optional): whether item discovery is used (default false)
//
// side effects: allocates vector of vectors to hold items
//
// returns: total number of payload items received
//
int Neighborhoods::ExchangeNeighbors(vector<vector<char *> > &items, float wf,
				     void (*ItemDtype)(MPI_Datatype*),
				     bool discovery) {

  // total number of neighbor blocks
  nn = 0;
  for (vector<bl_t>::iterator bi = blocks.begin(); bi != blocks.end(); bi++)
    nn += (int)bi->neighbors.size();
  PackMessages();
  PostMessages(ItemDtype);
  TestMessages(wf, ItemDtype);

  // copy received items from list to output vector
  int tot_npr = ListToVector(items, discovery);

  // clear the queued items
  for (vector<bl_t>::iterator bi = blocks.begin(); bi != blocks.end(); bi++) {
    for (vector<ne_t>::iterator ni = bi->neighbors.begin(); 
	 ni != bi->neighbors.end(); ni++) {
      ni->items.clear();
      ni->hdr.clear();
    }
  }

  return tot_npr;

}
//---------------------------------------------------------------------------
//
// packs items into packages for sending
// a package is one or more items going to the same process
// one package = 2 messages (a counts message and a payload message)
//
void Neighborhoods::PackMessages() {
  
  vector<pp_t>::iterator pi; // packages iterator

  // for all my blocks
  for (vector<bl_t>::iterator bi = blocks.begin(); bi != blocks.end(); bi++) {

    // for all the neighboring blocks of my blocks
    for (vector<ne_t>::iterator ni = (bi->neighbors).begin(); 
	 ni != (bi->neighbors).end(); 
	 ni++) {

      // look for destination process in all messages packed so far
      for (pi = pps.begin(); pi != pps.end(); pi++) {
	if (pi->proc == (ni->gb).proc && !pi->posted)
	  break;
      }

      // dest. proc. not found; add a new packed message
      if (pi == pps.end()) {
	pp_t pp; // a new package
	pp.proc = (ni->gb).proc;
	pp.posted = false;
	pp.c = (int*)malloc(sizeof(int));
	assert(pp.c != NULL);
	pp.c[0] = 0;
	pp.cn = 1;
	pps.push_back(pp); // invalidates iterator pi
	pi = pps.end() - 1; // reset pi to last item
      }

      // add count information

      // todo: realloc more than one int at a time
      pi->c = (int*)realloc(pi->c, ((pi->cn) + 2) * sizeof(int));

      assert(pi->c != NULL);
      pi->c[pi->cn] = (ni->gb).gid;
      pi->c[(pi->cn)+1] = (int)((ni->items).size());
      pi->cn = (pi->cn) + 2;
      pi->c[0]++;

      // add payload and header to the message
      int j = 0;
      for(vector<char *>::iterator ii = ni->items.begin(); 
	  ii != ni->items.end(); ii++) {
	pi->p.push_back(*ii); // payload
	for (vector<int>::iterator hi = (ni->hdr)[j].begin(); 
	     hi != (ni->hdr)[j].end(); hi++) {

	  // todo: realloc more than one int at a time
	  pi->c = (int*)realloc(pi->c, ((pi->cn) + 1) * sizeof(int));

	  pi->c[pi->cn] = *hi;
	  pi->cn = (pi->cn) + 1;
	}
	j++;
      }

    }
  }

}
//---------------------------------------------------------------------------
//
// posts counts-sends, points-sends, and count-receives messages
// these get posted first and don't depend on whether anything arrived
//
// ItemDtype: pointer to user-supplied function that creates an MPI datatype 
//  for an item to be sent or received
//
void Neighborhoods::PostMessages(void (*ItemDtype)(MPI_Datatype *)) {

  MPI_Request req;
  pl_t pl; // one payload message
  ct_t ct; // one count message
  int *rcv_ct; // one receive-counts message

  // post counts-sends, points-sends, counts-receives
  for (vector<pp_t>::iterator pi = pps.begin(); pi != pps.end(); pi++) {

    if (!pi->posted) {

      // counts-sends
      MPI_Isend(pi->c, pi->cn, MPI_INT, pi->proc, tag * 2, comm, &req);
      pi->posted = true;
      ct.req = req;
      ct.proc = pi->proc;
      ct.tag = tag * 2;
      ct.done = false;
      ct.c = NULL; // unused for sending
      send_cts.push_back(ct);

      // payload-sends
      if (pi->p.size() > 0) { // at least one item to send
	MPI_Datatype *itype = new MPI_Datatype;
	ItemDtype(itype);
	MPI_Datatype *mtype = SendMsgDtype(pi->c, &(pi->p)[0], itype);
	MPI_Isend(MPI_BOTTOM, 1, *mtype, pi->proc, tag * 2 + 1, comm, &req);
	pl.req = req;
	pl.proc = pi->proc;
	pl.done = false;
	send_pts.push_back(pl);
	MPI_Type_free(mtype);
	delete mtype;
	MPI_Type_free(itype);
	delete itype;
      }

      // counts-receives
      rcv_ct = new int[nn * (nhdr + 2) + 1];
      MPI_Irecv(rcv_ct, nn * (nhdr + 2) + 1, MPI_INT, pi->proc, tag * 2, comm, 
		&req);
      ct.req = req;
      ct.proc = pi->proc;
      ct.tag = tag * 2;
      ct.done = false;
      ct.c = rcv_ct;
      recv_cts.push_back(ct);

    }

  }

  tag++;

}
//---------------------------------------------------------------------------
//
// tests count-receive messages and posts payload-receive messages
// for those counts that arrived
//
// wf: wait_factor for nonblocking communication [0.0-1.0]
// ItemDtype: pointer to user-supplied function that creates an MPI datatype 
//  for an item to be sent or received
//
void Neighborhoods::TestMessages(float wf, 
				 void (*ItemDtype)(MPI_Datatype *)) {

  int npr; // number of received points from each process
  list<ct_t>::iterator ct_it; // request list iterators
  list<pl_t>::iterator pl_it; // request list iterators
  int p; // process number
  char *rcv_p = NULL; // one payload-receive
  int i, j, k;
  MPI_Request *reqs; // pending requests
  int *arr; // requests that arrived
  MPI_Status *stats; // statuses for arrivals
  int narr; // number of requests that arrived
  MPI_Status stat;
  int tot_narr = 0; // total number counts-receive messages arrived this round
  int nreqs; // number of requests
  if (!assign->GetStaticMode()) // override wf for dynamic repartitioning
    wf = 1.0f;
  int min_arr = (int)(wf * (int)pps.size()); // wait for this number of 
                                        // counts-receives
                                        // to arrive in this round

  if (recv_cts.size() > 0) {

    reqs = new MPI_Request[recv_cts.size()];

    arr = new int[recv_cts.size()];
    stats = new MPI_Status[recv_cts.size()];

    // wait for enough items in count-receive list to arrive
    while (tot_narr < min_arr) {

      nreqs = 0;

      for (ct_it = recv_cts.begin(); ct_it != recv_cts.end(); ct_it++) {
	if (!ct_it->done)
	  reqs[nreqs++] = ct_it->req;
      }

      if (nreqs) {

	MPI_Waitsome(nreqs, reqs, &narr, arr, stats);
	// post payload-receive for counts that arrived
	ct_it = recv_cts.begin();
	j = 0;

	for (i = 0; i < narr; i++) {

	  while (ct_it->done || j < arr[i]) { 
	    if (!ct_it->done)
	      j++;
	    ct_it++;
	  }

	  ct_it->done = true;
	  ct_it->tag = stats[i].MPI_TAG;

	  // count number of items expected
	  npr = 0;
	  for (k = 0; k < (ct_it->c)[0]; k++)
	    npr += (ct_it->c)[k * (nhdr + 2) + 2];

	  // post payload-receive
	  if (npr > 0) { // at least one point is expected

	    p = ct_it->proc;
	    MPI_Datatype *itype = new MPI_Datatype;
	    ItemDtype(itype);
	    MPI_Datatype *mtype = RecvMsgDtype(&(ct_it->c)[0], rcv_p, itype);
	    MPI_Recv(rcv_p, 1, *mtype, p, ct_it->tag + 1, comm, &stat);
	    pl_t pt; // one payload-receive message
	    pt.req = 0;
	    pt.done = true;
	    pt.proc = p;
	    pt.tag = ct_it->tag + 1; // matching tag for payload-receive
	    pt.p = rcv_p;
	    MPI_Aint lb, extent;
	    MPI_Type_get_extent(*itype, &lb, &extent);
	    pt.item_size = (int)extent;
	    recv_pts.push_back(pt);
	    MPI_Type_free(mtype);
	    delete mtype;
	    MPI_Type_free(itype);
	    delete itype;

	  } // if npr > 0

	  ct_it++;
	  j++;

	} // for i < narr

      } // if nreqs

      tot_narr += narr;

    } // tot_narr < min_narr

    delete[] reqs;
    delete[] arr;
    delete[] stats;

  } // recv_cts.size() > 0

}
//---------------------------------------------------------------------------
//
// flushes exchange with neighbors
//
// items: received items for each of my blocks [lid] (output)
// ItemDtype: pointer to user-supplied function that creates an MPI datatype 
//  for an item to be sent or received
// discovery (optional): whether item discovery is used (default false)
//
// returns: total number of payload items received
//
int Neighborhoods::FlushNeighbors(vector<vector<char *> > &items, 
				  void (*ItemDtype)(MPI_Datatype *),
				  bool discovery) {

  char* rcv_p; // one payload-receive
  MPI_Status stat;
  list<ct_t>::iterator ct_it; // request list iterators
  int npr; // number of points received
  int p; // process id
  int i;
  MPI_Request *reqs; // pending requests
  MPI_Status *stats; // statuss for arrivals
  int narr; // number of requests that arrived

  // wait for all items in counts-send list
  reqs = new MPI_Request[send_cts.size()];
  i = 0;

  for (ct_it = send_cts.begin(); ct_it != send_cts.end(); ct_it++) {
    if (!ct_it->done)
      reqs[i++] = ct_it->req;
  }

  stats = new MPI_Status[i];
  MPI_Waitall(i, reqs, stats);
  delete[] reqs;
  delete[] stats;

  // wait for all items in the counts-receive list
  reqs = new MPI_Request[recv_cts.size()];
  i = 0;

  for (ct_it = recv_cts.begin(); ct_it != recv_cts.end(); ct_it++) {
    if (!ct_it->done)
      reqs[i++] = ct_it->req;
  }

  stats = new MPI_Status[i];
  MPI_Waitall(i, reqs, stats);
  narr = i;

  // look for counts-receives that arrived and post points-receives for them
  ct_it = recv_cts.begin();

  for (i = 0; i < narr; i++) {

    while (ct_it->done)
      ct_it++;

    // count number of points expected
    npr = 0;

    for (int j = 0; j < (ct_it->c)[0]; j++)
      npr += (ct_it->c)[j * (nhdr + 2) + 2];

    // post points-receive
    if (npr > 0) { // at least one point is expected

      ct_it->done = true;
      ct_it->tag = stats[i].MPI_TAG;
      p = ct_it->proc;
      MPI_Datatype *itype = new MPI_Datatype;
      ItemDtype(itype);
      MPI_Datatype *mtype = RecvMsgDtype(&(ct_it->c)[0], rcv_p, itype);
      MPI_Recv(rcv_p, 1, *mtype, p, ct_it->tag + 1, comm, &stat);
      pl_t pt; // one point-receive message
      pt.req = 0; // not needed when doing blocking receive
      pt.proc = p;
      pt.tag = ct_it->tag + 1; // matching tag for point-receive
      pt.done = true;
      pt.p = rcv_p;
      MPI_Aint lb, extent;
      MPI_Type_get_extent(*itype, &lb, &extent);
      pt.item_size = (int)extent;
      recv_pts.push_back(pt);
      MPI_Type_free(mtype);
      delete mtype;
      MPI_Type_free(itype);
      delete itype;

    } // npr > 0

    ct_it++;

  } // for

  delete[] reqs;
  delete[] stats;

  // wait for all items in payloads-send list
  reqs = new MPI_Request[send_pts.size()];
  i = 0;
  for (list<pl_t>::iterator pl_it = send_pts.begin(); 
       pl_it != send_pts.end(); pl_it++) {
    if (!pl_it->done)
      reqs[i++] = pl_it->req;
  }
  stats = new MPI_Status[i];
  MPI_Waitall(i, reqs, stats);
  delete[] reqs;
  delete[] stats;

  // copy counts and points from receive list to output array
  int tot_npr = ListToVector(items, discovery);

  // cleanup
  for (vector<pp_t>::iterator ppi = pps.begin(); ppi != pps.end(); ppi++) {
    for (vector<char *>::iterator pi = ppi->p.begin(); pi != ppi->p.end();
	 pi++)
      delete[] *pi; // item payload, new'ed back when item was enqueued
    free(ppi->c);
    ppi->c = NULL;
    ppi->cn = 0;
    ppi->p.clear();
  }
  pps.clear();

  for (list<ct_t>::iterator ci = recv_cts.begin(); ci != recv_cts.end(); ci++)
    delete[] ci->c; // array of received counts
  for (list<pl_t>::iterator pl_it = recv_pts.begin(); pl_it != recv_pts.end(); 
       pl_it++) {
    if(pl_it->p != NULL)
      delete[] pl_it->p; // array of received items
  }

  send_cts.clear();
  send_pts.clear();
  recv_cts.clear();
  recv_pts.clear();

  tag = 0;
  return tot_npr;

}
//---------------------------------------------------------------------------
//
// (shallow) copies payloads from list to output vector
//
// items: vector of items for each of my blocks [lid][item]
// discovery: whether item discovery is used
//
// side effects: allocates vector of vectors to hold items
//
// returns: total number of items received
//
int Neighborhoods::ListToVector(vector<vector<char *> >& items, 
				bool discovery) {

  list<ct_t>::iterator ct_it;
  int tot_npr = 0; // total number of points received
  int lid, gid; // local and global block index

  items.resize(assign->NumBlks()); // start with a vector for each block

  list<pl_t>::iterator pl_it = recv_pts.begin();
  while (recv_pts.size() > 0 && pl_it != recv_pts.end()) {
    if (pl_it->done) {
      // find the matching count request
      for (ct_it = recv_cts.begin(); ct_it != recv_cts.end(); ct_it++) {
	if (ct_it->proc == pl_it->proc && ct_it->tag + 1 == pl_it->tag)
	  break;
      }
      assert(ct_it != recv_cts.end() && ct_it->done); // sanity
      // traverse the counts message
      int n = 0;
      for (int j = 0; j < (ct_it->c)[0]; j++) {

	if (ct_it->c[j * 2 + 2] > 0) { // there is something to copy

	  gid = (ct_it->c)[2 * j + 1]; // global and local block id of dest.

	  // special discovery version of Gid2Lid is in neighbor class,
	  // otherwise, plain version in blocking class
	  if (discovery) {

	    gb_t *gb = (gb_t *)&(pl_it->p[j * pl_it->item_size]);
	    vector<int> lids;
	    int ret_val = DiscoverGid2Lid(gb, lids);
	    assert(ret_val >= 0);

	    // copy items from recv_pts to >= 1 local block of output items
	    for (int k = 0; k < (ct_it->c)[j * (nhdr + 2) + 2]; k++) {
	      char *pp = &(pl_it->p[n]);
	      for (int i = 0; i < (int)lids.size(); i++)
		items[lids[i]].push_back(pp);
	      tot_npr++;
	      n += pl_it->item_size;
	    }

	  } // discovery mode

	  else { // normal mode

	    lid = blocking->Gid2Lid(gid);
	    assert(lid >= 0);

	    // copy items from recv_pts to one local block of output items
	    for (int k = 0; k < (ct_it->c)[j * (nhdr + 2) + 2]; k++) {
	      char *pp = &(pl_it->p[n]);
	      items[lid].push_back(pp);
	      tot_npr++;
	      n += pl_it->item_size;
	    }

	  } // normal mode

	} // there is something new to copy

      } // traverse the counts message

      // delete the messages

      // individual payload message needs to be deleted by caller
      pl_it = recv_pts.erase(pl_it);

      // counts message deleted here, not needed anymore
      delete[] ct_it->c;
      ct_it = recv_cts.erase(ct_it);

    } // done

    else
      pl_it++;

  } // recv_pts

  return tot_npr;

}
//---------------------------------------------------------------------------
//
// neighbor discovery version of searching local blocks for received 
//   global block
//  If gid is not found, checks local vids against incoming vid and 
//   expected process source against actual process source. If these match,
//   returns local lid of receiving block
// 
// gb: incoming global block
// lids: one or more local blocks that are expecting the received global block
//
// returns: 1 upon success, -1 if not found
//
int Neighborhoods::DiscoverGid2Lid(gb_t *gb, vector<int> &lids) {

  int i, j, m;

  lids.clear(); // just in case the caller did not send us an empty list

  for (i = 0; i < (int)blocks.size(); i++) {

    // first chance: check if I have the incoming gid
    if (blocks[i].gid == gb->gid) {
      lids.push_back(i);
      continue;
    }

    // second chance: check if incoming process and vid are in my memote ids
    bool found = false;
    for (j = 0; j < gb->num_vids; j++) {
      found = false;
      for (m = 0; m < num_rem_ids[i]; m++) {
	if (rem_ids[i][m].vid == gb->vids[j] && 
	    rem_ids[i][m].proc == gb->proc) {
	  lids.push_back(i);
	  found = true;
	  break;
	}
      }
      if (found)
	break;
    }

    // update the neighbor gid, perhaps adding a new neighbor
    if (found) {

      // check if a dummy neighbor exists that can be updated
      int nn;
      for (nn = 0; nn < (int)blocks[i].neighbors.size(); nn++) {
	if (blocks[i].neighbors[nn].gb.gid == -1 &&
	    blocks[i].neighbors[nn].gb.proc == gb->proc) {
	  blocks[i].neighbors[nn].gb.gid = gb->gid;
	  break;
	}
      }

      // add a new neighbor
      if (nn == (int)blocks[i].neighbors.size()) {
	ne_t new_neigh;
	new_neigh.gb.gid = gb->gid;
	new_neigh.gb.proc = gb->proc;
	blocks[i].neighbors.push_back(new_neigh);
      }

    } // update the neighbor gid

  }

  return(lids.size() ? 1 : -1);

}
//---------------------------------------------------------------------------
//
// makes MPI datatype for a payloads-receive message and allocates the receive
//  buffer
//
// cts: pointer to counts message
// pts: pointer to payloads message
//
// side effects: allocates points message if the counts message has >= 1 item
//               allocates MPI datatype
//               commits MPI datatype
//
// returns: pointer to MPI datatype
//
MPI_Datatype* Neighborhoods::RecvMsgDtype(int *cts, char* &pts,
					  MPI_Datatype *itype) {

  int npr = 0; // number of items received
  for (int i = 0; i < cts[0]; i++)
    npr += cts[i * 2 + 2];

  MPI_Datatype *mtype = new MPI_Datatype;
  MPI_Type_contiguous(npr, *itype, mtype);
  MPI_Type_commit(mtype);
  MPI_Aint extent; // datatype size in bytes
  MPI_Type_extent(*mtype, &extent);
  if (extent > 0) // allocate receive buffer
    pts = new char[extent];

  return mtype;

}
//-----------------------------------------------------------------------
//
// makes an MPI datatype for a payloads-send message
//
// cts: pointer to counts message
// pts: pointer to payloads message
// itype: datatype of a single item
//
// side effects: allocates MPI datatype
//               commits MPI datatype
//
// returns: pointer to MPI datatype
//
//
MPI_Datatype* Neighborhoods::SendMsgDtype(int *cts, char **pts, 
					  MPI_Datatype *itype) {

  int nps = 0; // number of items being sent
  for (int i = 0; i < cts[0]; i++)
    nps += cts[i * 2 + 2];

  // lengths and displacements array
  int *lengths = new int[nps];
  MPI_Aint *disps = new MPI_Aint[nps];
  for (int i = 0; i < nps; i++) {
    lengths[i] = 1;
    MPI_Get_address((void *)pts[i], &disps[i]);
  }

  MPI_Datatype *mtype = new MPI_Datatype; // datatype for entire message
  MPI_Type_create_hindexed(nps, lengths, disps, *itype, mtype);
  MPI_Type_commit(mtype);

  delete[] lengths;
  delete[] disps;

  return mtype;

}
//-----------------------------------------------------------------------
//
// exchanges all block bounds among all neighbors
// neighbor bounds always exclude ghost
//
void Neighborhoods::GetNeighborBounds() {

  int nblocks = assign->NumBlks();
  gb_t gb; // global bounds object
  gb_t *pgb; // pointer to a global bounds object

  for (int b = 0; b < nblocks; b++)  {

    // prepare a global block for sending
    gb.gid = blocking->Lid2Gid(b); // my gid
    gb.proc = rank; // my pid
    gb.num_vids = 0; // unused
    bb_t bounds; // my block bounds
    blocking->NoGhostBlockBounds(b, &bounds);
    for (int i = 0; i < blocking->GetDim(); i++) {
      gb.bb.min[i] = bounds.min[i];
      gb.bb.max[i] = bounds.max[i];
    }
    for (int i = blocking->GetDim(); i < DIY_MAX_DIM; i++) {
      gb.bb.min[i] = 0; // initialize the unused dims to be safe
      gb.bb.max[i] = 0;
    }

    // enqueue the item to send to all neighbors
    EnqueueItemAll(b, (char *)&gb, sizeof(gb_t), NULL, NULL, true, true);

  }

  vector<vector <char *> > items; // received items

  ExchangeNeighbors(items, 1.0f, &Nbhds_ItemType);

  assert(nblocks == (int)items.size()); // sanity

  // store received block bounds
  for (int b = 0; b < (int)items.size(); b++) {
    for (int n = 0; n < (int)items[b].size(); n++) {
      pgb = (gb_t *)items[b][n];
      // debug
//       fprintf(stderr, "b = %d n = %d gid = %d vid = %d p = %d "
// 	      "mins = %.1f %.1f %.1f %.1f maxs = %.1f %.1f %.1f %.1f\n", 
// 	      b, n, pgb->gid, pgb->vid, pgb->proc, 
// 	      pgb->bb.min[0], pgb->bb.min[1], pgb->bb.min[2], pgb->bb.min[3], 
// 	      pgb->bb.max[0], pgb->bb.max[1], pgb->bb.max[2], pgb->bb.max[3]);
      SetNeighborBounds(b, pgb);

    }
  }

  // cleanup
  // todo: this fails when nb > 1 because several items may have been allocated
  // in one message, and can't be freed inidividually
  // todo: find a better solution
  if (assign->NumBlks() == 1) {
    for (int b = 0; b < (int)items.size(); b++) {
      for (int n = 0; n < (int)items[b].size(); n++) {
	delete[] items[b][n];
      }
    }
  }

  // precautionary, should be nothing left to receive
  FlushNeighbors(items, &Nbhds_ItemType);

}
//-----------------------------------------------------------------------
//
// sets up a valid neighbors structure
//
// vids: local vertex ids for each local block that needs neighbor discovery
// num_vids: number of vids for each local block
//
void Neighborhoods::GetNeighbors(int **vids, int *num_vids) {

  int nblocks = assign->NumBlks();
  gb_t gb; // global bounds object
  gb_t *pgb; // pointer to a global bounds object

  for (int b = 0; b < nblocks; b++)  {

    // prepare a global block for sending
    gb.gid = blocking->Lid2Gid(b); // my gid
    gb.proc = rank; // my (sending) proc
    if (vids && num_vids) { // neighbor discovery is used
      assert(num_vids[b] <= DIY_MAX_VIDS); // todo: dynamic allocate?
      for (int i = 0; i < num_vids[b]; i++)
	gb.vids[i] = vids[b][i];
      gb.num_vids = num_vids[b];
    }
    else // no neighbor discovery
      gb.num_vids = 0;
    bb_t bounds; // my block bounds
    blocking->BlockBounds(b, &bounds);
    for (int i = 0; i < blocking->GetDim(); i++) {
      gb.bb.min[i] = bounds.min[i];
      gb.bb.max[i] = bounds.max[i];
    }
    for (int i = blocking->GetDim(); i < DIY_MAX_DIM; i++) {
      gb.bb.min[i] = 0; // initialize the unused dims to be safe
      gb.bb.max[i] = 0;
    }

    // enqueue the item to send to all neighbors
    EnqueueItemAll(b, (char *)&gb, sizeof(gb_t), NULL, NULL, true, true);

    // debug
//     if (rank == 0) {
//       if (num_rem_ids[b] == 10) {
// 	fprintf(stderr, "gid = %d mins = %.1f %.1f %.1f maxs = %.1f %.1f %.1f\n", 
// 		gb.gid, gb.bb.min[0], gb.bb.min[1], gb.bb.min[2], 
// 		gb.bb.max[0], gb.bb.max[1], gb.bb.max[2]);
// 	fprintf(stderr, "%d vids: ", gb.num_vids);
// 	for(int j = 0; j < gb.num_vids; j++)
// 	  fprintf(stderr, "%d ", gb.vids[j]);
// 	fprintf(stderr, "\n");
// 		for(int j = 0; j < num_rem_ids[b]; j++)
// 		  fprintf(stderr, "[vid, proc] = [%d, %d] ", 
// 			  rem_ids[b][j].vid, rem_ids[b][j].proc);
// 		fprintf(stderr, "\n");
//       }
//     }
  }

  vector<vector <char *> > items; // received items

  ExchangeNeighbors(items, 1.0f, &Nbhds_ItemType, discovery);

  assert(nblocks == (int)items.size()); // sanity

  // store received block bounds
  for (int b = 0; b < (int)items.size(); b++) {
    for (int n = 0; n < (int)items[b].size(); n++) {
      pgb = (gb_t *)items[b][n];

      // debug
//       fprintf(stderr, "b = %d n = %d gid = %d p = %d "
// 	      "mins = %.1f %.1f %.1f %.1f maxs = %.1f %.1f %.1f %.1f\n", 
// 	      b, n, pgb->gid, pgb->proc, 
// 	      pgb->bb.min[0], pgb->bb.min[1], pgb->bb.min[2], pgb->bb.min[3], 
// 	      pgb->bb.max[0], pgb->bb.max[1], pgb->bb.max[2], pgb->bb.max[3]);

      SetNeighborBounds(b, pgb);

    }
  }

  // precautionary, should be nothing left to receive
  FlushNeighbors(items, &Nbhds_ItemType, true);

}
//-----------------------------------------------------------------------
//
// sets block bounds for a neighboring block of my local block
//
// lid: local block id
// gb: pointer to received global block
//
void Neighborhoods::SetNeighborBounds(int lid, gb_t *gb) {

  int n; // neighbor number

  // Note that we don't break out of this loop after finding a matching
  // neighbor, because several can match when wraparound neighbors are used.
  for (n = 0; n < (int)blocks[lid].neighbors.size(); n++) {
    if (gb->gid == blocks[lid].neighbors[n].gb.gid) {

      assert(gb->proc == blocks[lid].neighbors[n].gb.proc); // sanity

      for (int j = 0; j < DIY_MAX_DIM; j++) {
	blocks[lid].neighbors[n].gb.bb.min[j] = gb->bb.min[j];
	blocks[lid].neighbors[n].gb.bb.max[j] = gb->bb.max[j];
      }

    }
  }

}
//---------------------------------------------------------------------------
//
// returns global id of (first) neighboring block found containing the point
//
// lid: local id of current block
// pt: point (up to number of dimensions specified when creating the class)
//
// returns: neighboring block global id or -1 if out of the domain
//
int Neighborhoods::Pt2NeighGid(int lid, float *pt) {

  bool is_in; // the point is tentatively in the block

  for (int n = 0; n < (int)blocks[lid].neighbors.size(); n++) {

    is_in = true;
    for (int d = 0; d < blocking->GetDim(); d++) {

      if (pt[d] < blocks[lid].neighbors[n].gb.bb.min[d] || 
	  pt[d] > blocks[lid].neighbors[n].gb.bb.max[d]) {
	is_in = false;
	break;
      }

    }

    if (is_in)
      return(blocks[lid].neighbors[n].gb.gid);

  }

  return(-1); 

}
//--------------------------------------------------------------------------
//
// makes MPI datatype for sending and receiving one item
//
// type: pointer to MPI datatype
//
static void Nbhds_ItemType(MPI_Datatype *type) {

  // datatype for block bounds
  MPI_Datatype btype;
  struct map_block_t bmap[] = {
    { MPI_FLOAT, OFST, DIY_MAX_DIM, offsetof(struct bb_t, min) },
    { MPI_FLOAT, OFST, DIY_MAX_DIM, offsetof(struct bb_t, max) },
  };
  DIY_Create_struct_datatype(0, 2, bmap, &btype);

  // datatype for global block
  struct map_block_t map[] = {
    { MPI_INT,           OFST, 1,            offsetof(struct gb_t, gid)      },
    { MPI_INT,           OFST, DIY_MAX_VIDS, offsetof(struct gb_t, vids)     },
    { MPI_INT,           OFST, 1,            offsetof(struct gb_t, num_vids) },
    { MPI_INT,           OFST, 1,            offsetof(struct gb_t, proc)     },
    { MPI_UNSIGNED_CHAR, OFST, 1,            offsetof(struct gb_t, neigh_dir)},
    { btype,             OFST, 1,            offsetof(struct gb_t, bb)       },
};
  DIY_Create_struct_datatype(0, 6, map, type); 
  DIY_Destroy_datatype(&btype);

}
//-----------------------------------------------------------------------
//
// finds neighbors that intersect bounds +/- an extra amount t
//
// lid: local block id
// bounds: target bounds
// t: additional extension on all sides of bounds
// num_intersect (output) number of intersecting neighbors found
// gids_intersect (output) the intersecting neighbor block gids
//
void Neighborhoods::BoundsIntersectNeighbors(int lid, bb_t cell_bounds, 
					     float t, int *num_intersect, 
					     int *gids_intersect) {

  *num_intersect = 0;

  for (int n = 0; n < (int)blocks[lid].neighbors.size(); n++) {

    bool is_intersect = false;
    for (int d = 0; d < blocking->GetDim(); d++) {

      if (! ((cell_bounds.min[d] < blocks[lid].neighbors[n].gb.bb.min[d]-t && 
	      cell_bounds.max[d] < blocks[lid].neighbors[n].gb.bb.min[d]-t) ||
	     (blocks[lid].neighbors[n].gb.bb.min[d]-t < cell_bounds.min[d] &&
	      blocks[lid].neighbors[n].gb.bb.max[d]+t < cell_bounds.min[d])) ) {
	is_intersect = true;
	break;
      }

    }

    if (is_intersect) {
      gids_intersect[*num_intersect] = blocks[lid].neighbors[n].gb.gid;
      (*num_intersect)++;
    }  

  }

}
//-----------------------------------------------------------------------
