//---------------------------------------------------------------------------
//
// diy wrappers, callable from C, C++, and Fortran
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
#include "diy.h"
#include "blocking.hpp"
#include "io.hpp"
#include "merge.hpp"
#include "swap.hpp"
#include "bil.h"
#include "neighborhoods.hpp"
#include "comm.hpp"
#include "util.hpp"

// diy global data
static int64_t sizes[DIY_MAX_DIM]; // grid size (extents.max - extents.min + 1)
static int tb; // total number of blocks in the domain
static int nb; // my number of blocks
static int maxb; // max. number of blocks in any process
static int dim; // number of dimensions
static MPI_Comm comm; // MPI communicator
static int rank; // MPI rank of this process
static int groupsize; // MPI groupsize
static int rounds; // number of merge rounds
static int kvs[DIY_MAX_R]; // merge k-values

Assignment *assign; // assignment object (abstract base class)
RoundRobinAssignment *round_robin_assign; // round robin assignment object
ProcOrderAssignment *proc_order_assign; // process order assignment object
Blocking *blocking; // blocking object
Merge *merging; // merge object
Swap *swapping; // swap object
IO *io; // output I/O object
Neighborhoods *nbhds; // neighborhoods object
Comm *cc; // communication object

// DIY datatypes are aliases for MPI datatypes
// DIY doesn't distinguish between signed and unsigned, lengths are the same
// DIY supports short and long, but their use is not recommended because lengths
// can vary across architectures
DIY_Datatype DIY_BYTE = MPI_BYTE; /* 1 byte */
DIY_Datatype DIY_SHORT = MPI_SHORT; /* 2 bytes */
DIY_Datatype DIY_INT = MPI_INT; /* 4 bytes */
DIY_Datatype DIY_LONG = MPI_LONG; /* 4 bytes */
DIY_Datatype DIY_LONG_LONG = MPI_LONG_LONG_INT; /* 8 bytes */
DIY_Datatype DIY_FLOAT = MPI_FLOAT; /* 4 bytes */
DIY_Datatype DIY_DOUBLE = MPI_DOUBLE; /* 8 bytes */
DIY_Datatype DIY_LONG_DOUBLE = MPI_LONG_DOUBLE; /* 16 bytes */

/* DIY_BOTTOM is an alias for MPI_BOTTOM */
void *DIY_BOTTOM = MPI_BOTTOM;

// vector version of results of neighbor exchange
vector<vector< char *> > items_v;

// vector versions of block compresseion / decompression buffers
vector<unsigned char> comp_buf_v;
vector<unsigned char> decomp_buf_v;

//--------------------------------------------------------------------------
//
// most functions return an error code
// 0 indicates succes
// > 0 indicates various errors, meaning is function specific
//
// currently all functions return 0 (success)
// todo: add error checking
//
//--------------------------------------------------------------------------
//
// Initializes DIY
//
// dim: number of dimensions (2, 3, or 4) (input)
// data_size: data size in each dimension (only for structured data), pass NULL
//  for other cases
// num_threads: number of threads DIY is allowed to use (>= 1)
//  ignored if DIY is built with openmp disabled (enabled by default)
// comm: MPI cpmmunicator
//
// returns: error code
//
int DIY_Init(int dim, int *data_size, int num_threads, MPI_Comm comm) {

  nbhds = NULL;
  ::dim = dim;
  ::comm = comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);
  rounds = 0;

  if (data_size) {
    for (int i = 0; i < dim; i++)
      sizes[i] = data_size[i];
  }

  BIL_Init(comm);
  cc = new Comm(comm);
  merging = new Merge(comm);
  swapping = new Swap(comm);

#ifdef OMP
  omp_set_num_threads(num_threads);
#endif

  return 0;

}
//--------------------------------------------------------------------------
//
// Decomposes the domain
//
// block_order: ROUND_ROBIN_ORDER or CONTIGUOUS_ORDER numbering of 
//   global block idss to processes (input)
// glo_num__blocks: total number of blocks in the global domain (input)
//   pass 0 or anything for CONTIGUOUS_ORDER (unused)
// loc_num_blocks: local number of blocks on this process (output)
// share_face: whether neighboring blocks share a common face or are
//  separated by a gap of one unit
// ghost: ghost layer for each dimension and side (min, max)
//  each entry can be 0 (no ghost) or > 0 (this many ghost cells per side)
//   {x min side ghost, x max side ghost, y min side ghost, y max side ghost...}
// given: constraints on the blocking entered as an array where
//   0 implies no constraint in that direction and some value n > 0 is a given
//   number of blocks in a given direction
//   eg., {0, 0, 0, t} would result in t blocks in the 4th dimension
//
//
// returns: error code
//
int DIY_Decompose(int block_order, int glo_num_blocks, int *loc_num_blocks, 
		  int share_face, int *ghost, int *given) {

  if (block_order == ROUND_ROBIN_ORDER)
    assign = 
      new RoundRobinAssignment(glo_num_blocks, nb, maxb, comm);
  else
    assign = 
      new ProcOrderAssignment(glo_num_blocks, nb, maxb, comm);
  *loc_num_blocks = nb;
  ::tb = glo_num_blocks;

  int64_t given64[DIY_MAX_DIM]; // int64_t version of given
  for (int i = 0; i < dim; i++)
    given64[i] = given[i];

  // decomposing a new domain
  blocking = new Blocking(dim, tb, sizes, share_face, ghost, 
			  given64, assign, comm);
  nbhds = new Neighborhoods(blocking, assign, comm, false);
  io = new IO(dim, tb, maxb, comm);

  return 0;

}
//--------------------------------------------------------------------------
//
// Describes the already decomposed domain
//
// loc_num_blocks: local number of blocks on this process
// gids: global ids of my local blocks
// bounds: block bounds (extents) of my local blocks
// rem_ids: remote ids used for neighbor discovery (pass NULL if 
//  neighbor discovery not needed at all)
// num_rem_ids: number of remote ids for each local block (pass 0 for any
//  blocks that don't need neighbor discovery)
// vids: local vertex ids for each local block that needs neighbor discovery
//  (pass NULL if not needed at all)
// num_vids: number of vids for each local block (pass 0 for any blocks that
//  don't need neighbor discovery)
// neighbors: neighbor lists for each of my local blocks, in lid order
//  lists can be partial, containing only known information, as long as 
//  the rem_data lists contain enough information for DIY to discover the rest
//  of the unknown neighbors. Unknown (remote) block bounds can be 
//  uninitialized. If wrapping is used, (see below) a neighbor direction 
//  must be provided for each neighbor
// num_neighbors: number of neighbors known so far for each of my 
//  local blocks, in lid order
// wrap: whether wraparound neighbors are used (0 = no wraparound neighbors
//  were provided, 1 = provided some wraparound neighbors)
//
// returns: error code
//
int DIY_Decomposed(int loc_num_blocks, int *gids, struct bb_t *bounds, 
		   struct ri_t **rem_ids, int *num_rem_ids, int **vids, 
		   int *num_vids, struct gb_t **neighbors, 
		   int *num_neighbors, int wrap) {

  assign = new ExistingAssignment(loc_num_blocks, maxb, tb, comm);

  nb = loc_num_blocks;

  // existing decomposition version of blocking
  blocking = new Blocking(dim, tb, gids, bounds, assign, comm);

  // discovery version of constructing the neighborhoods
  nbhds = new Neighborhoods(blocking, assign, rem_ids, num_rem_ids, vids,
			    num_vids, neighbors, num_neighbors, comm, wrap);
  io = new IO(dim, tb, maxb, comm);

  return 0;

}
//--------------------------------------------------------------------------
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
// returns: error code
//
int DIY_Block_starts_sizes(int lid, int *starts, int *sizes) {

  int64_t starts64[DIY_MAX_DIM];
  int64_t sizes64[DIY_MAX_DIM];

  blocking->BlockStartsSizes(lid, starts64, sizes64);
  for (int i = 0; i < dim; i++) {
    starts[i] = starts64[i];
    sizes[i] = sizes64[i];
  }

  return 0;

}
//--------------------------------------------------------------------------
//
// block bounds including ghost
// for blocks consisting of continuous spatiotemporal regions
//
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
// returns: error code
//
int DIY_Block_bounds(int lid, struct bb_t *bounds) {


  blocking->BlockBounds(lid, bounds);

  return 0;

}
//--------------------------------------------------------------------------
//
// block bounds excluding ghost
// for blocks consisting of continuous spatiotemporal regions
//
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
// returns: error code
//
int DIY_No_ghost_block_bounds(int lid, struct bb_t *bounds) {


  blocking->NoGhostBlockBounds(lid, bounds);

  return 0;

}
//--------------------------------------------------------------------------
//
// finds the time block to which a local block belongs
//
// lid: local block id
// time_block: time block containing lid (output)
//
// return: error code
//
int DIY_In_time_block(int lid, int *time_block) {

  int64_t lat_nblocks[DIY_MAX_DIM]; // number of blocks in each dimension
  blocking->NumLatBlocks(lat_nblocks);

  for (int g = 0; g < tb; g++) {
    if (blocking->InTimeBlock(g, lid, sizes[3], lat_nblocks[3])) {
      *time_block = g;
      return 0;
    }
  }

  return 1;     

}
//-----------------------------------------------------------------------
//
// posts a data read for a block
//
// block_starts: starting grid indices for minimum corner of block
// block_sizes: size of block in each dimension (number of vertices)
// file_name: input file name
// var_type: input datatype
// buffer: pointer to void * data buffer
//
// returns: error code
//
int DIY_Add_data_raw(int *block_starts, int *block_sizes, char *file_name, 
		     DIY_Datatype var_type, void** buffer) {

  // reverse order for BIL
  int bil_data_size[3];
  int bil_min[3];
  int bil_size[3];
  assert(dim >= 3); // don't know if bil does 2D, probably not
  int ofst = (dim == 3 ? 0 : 1); // offset for reversing order
  for (int i = 0; i < dim - ofst; i++) {
    bil_data_size[i] = sizes[dim - ofst - i - 1];
    bil_min[i] = block_starts[dim - ofst - i - 1];
    bil_size[i] = block_sizes[dim - ofst - i - 1];
  }

  // BIL always posts 3d blocks (1st arg.), even when dim = 4, 
  // it posts a sequence of 3d blocks
  // todo: does it support 2d blocks?
  BIL_Add_block_raw(3, bil_data_size, bil_min, bil_size, file_name, 
		    var_type, buffer);

  return 0;

}
//--------------------------------------------------------------------------
//
// executes a parallel data read of all blocks
//
// note: performs an MPI_Barrier afterwards to eliminate any process
//   skew before proceeding
//
// returns: error code
//
int DIY_Read_data_all() {

  BIL_Read();
  MPI_Barrier(comm);

  return 0;

}
//--------------------------------------------------------------------------
//
// configurable in-place merge reduction
//
// blocks: pointers to input/output blocks, results in first num_blocks_out
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// num_rounds: number of rounds
// k_values: radix (group size) for each round
// reduce_func: pointer to merging or reduction function
// create_func: pointer to function that creates item
// destroy_func: pointer to function the destroys item
// type_func: pointer to function that creates MPI datatype for item 
//   returns the base address associated with the datatype
// num_blocks_out: number of output blocks (output)
//
// side effects: allocates output items and array of pointers to them
//  if not reducing in-place
//
// returns: error code
//
int DIY_Merge_blocks(char **blocks, int **hdrs, int num_rounds, int *k_values,
		     void (*reduce_func)(char **, int *, int), 
		     char *(*create_func)(int *),
		     void (*destroy_func)(void *),
		     void*(*type_func)(void*, DIY_Datatype*),
		     int *num_blocks_out) {

  rounds = num_rounds;
  for (int i = 0; i < num_rounds; i++)
    kvs[i] = k_values[i];

  *num_blocks_out = merging->MergeBlocks(blocks, hdrs, num_rounds, k_values,
					 cc, assign, reduce_func,
					 create_func, destroy_func, type_func);

  return 0;

}
//--------------------------------------------------------------------------
//
// configurable in-place swap reduction
//
// blocks: pointers to input/output blocks
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// num_elems: number of elements in a block
// num_rounds: number of rounds
// k_values: radix (group size) for each round
// starts: start of result in each block (output)
// sizes: number of elements in result in each block (output)
//  starts and sizes are allocated by the caller
// reduce_func: pointer to reduction function
// recv_create_func: pointer to function that creates received item 
//   with given number of elements (less than original item)
// recv_destroy_func: pointer to function that destroys received item
//  part of a total number of parts in the item
// send_type_func: pointer to function that creates MPI datatype for sending
//   a subset of the item starting at an element and having a given number
//   of elements (less than the original item)
//   returns the base address associated with the datatype
// recv_type_func: pointer to function that creates MPI datatype for receiving
//   a subset of the item with a given number of elements 
//   (less than the original item)
//   returns the base address associated with the datatype
//
// returns: error code
//
int DIY_Swap_blocks(char **blocks, int **hdrs, int num_elems,
		    int num_rounds, int *k_values, int *starts, int *sizes,
		    void (*reduce_func)(char **, int *, int, int), 
		    char *(*recv_create_func)(int *, int),
		    void (*recv_destroy_func)(void *),
		    void*(*send_type_func)(void*, DIY_Datatype*, int, int),
		    void*(*recv_type_func)(void*, DIY_Datatype*, int)) {

  rounds = num_rounds;
  for (int i = 0; i < num_rounds; i++)
    kvs[i] = k_values[i];

  swapping->SwapBlocks(blocks, hdrs, num_rounds, k_values, num_elems, 
		       starts, sizes, cc, assign, reduce_func,
		       recv_create_func, recv_destroy_func, send_type_func, 
		       recv_type_func);

  return 0;

}
//--------------------------------------------------------------------------
//
// initializes parallel writing of analysis blocks
//
// filename: output filename
// compress: whether to compress output (0 = normal, 1 = compress)
//   (1: zlib's default compression level 6 is applied blockwise)
//
// returns: eerror code
//
int DIY_Write_open_all(char *filename, int compress) {

  MPI_Barrier(comm); // synchronize to clear any skew before proceeding
  io->WriteAnaInit(filename, compress);

  return 0;

}
//----------------------------------------------------------------------------
//
// writes all analysis blocks in parallel with all other processes
//
// blocks: array of pointers to analysis blocks
// num_blocks: number of blocks
// hdrs: headers, one per analysis block (NULL if not used)
// num_hdr_elems; number of header elements (0 if not used), 
//   same for all headers
// type_func: pointer to function that creates MPI datatype for item 
//   returns the base address associated with the datatype
//
// returns: error code
//
int DIY_Write_blocks_all(void **blocks, int num_blocks, int **hdrs,
			 int num_hdr_elems, 
			 void*(*type_func)(void*, int, DIY_Datatype*)) {

  int tot_nb_merged = tb; // max possible number of merged blocks
  for (int i = 0; i < rounds; i++)
    tot_nb_merged /= kvs[i];

  // max number of blocks in any process after possible merging
  int max_blocks = maxb < tot_nb_merged ? maxb : tot_nb_merged;

  io->WriteAllAna(blocks, num_blocks, max_blocks, hdrs, num_hdr_elems, 
		  type_func);

  return 0;

}
//----------------------------------------------------------------------------
//
// finalizes parallel writing of analysis blocks
//
// returns: error code
//
int DIY_Write_close_all() {

  MPI_Barrier(comm); // flushes all I/O
  io->WriteAnaFinalize();

  return 0;

}
//----------------------------------------------------------------------------
//
// initializes parallel reading of analysis blocks
//
// filename: input filename
// swap_bytes: whether to swap bytes for endian conversion
//  only applies to reading the headers and footer
//  user must swap bytes manually for datatypes because they are custom
// compress: whether to compress output (0 = normal, 1 = compress)
//   (1: zlib's default compression level 6 is applied blockwise)
//
// returns: eerror code
//
int DIY_Read_open_all(char *filename, int swap_bytes, int compress) {

  MPI_Barrier(comm); // synchronize to clear any skew before proceeding
  io->ReadAnaInit(filename, swap_bytes, compress);

  return 0;

}
//----------------------------------------------------------------------------
//
// reads all analysis blocks in parallel with all other processes
//
// blocks: pointer to array of pointers for analysis blocks being read (output)
//  DIY will allocate blocks for you
// num_blocks: number of local blocks read (output)
// hdrs: headers, one per analysis block, allocated by caller
//   (pass NULL if not used)
// create_type_func: pointer to function that takes a block local id, 
//   block header, and creates (allocates) a block and creates an MPI datatype 
//   for it. Returns the base address associated with the datatype
//
// returns: error code
//
int DIY_Read_blocks_all(void ***blocks, int *num_blocks, int **hdrs,
			void* (*create_type_func)(int, int*, DIY_Datatype*)) {

  *num_blocks = io->ReadAllAna(*blocks, hdrs, create_type_func);

  return 0;

}
//----------------------------------------------------------------------------
//
// finalizes parallel reading of analysis blocks
//
// returns: error code
//
int DIY_Read_close_all() {

  MPI_Barrier(comm); // flushes all I/O
  io->ReadAnaFinalize();

  return 0;

}
//----------------------------------------------------------------------------
//
// exchanges items with all neighbors
//
// items: pointer to received items for each of my blocks [lid][item] (output)
// num_items: number of items for each block (allocated by user)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// suggested value: 0.1
// RecvItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message  and
//   creates an MPI datatype for the payloads message
// SendItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and a payloads message and
//   creates an MPI datatype for the payloads message
//
// side effects: allocates items and array of pointers to them
//
// returns: error code
//
int DIY_Exchange_neighbors(void ***items, int *num_items, float wf,
			   DIY_Datatype* (*RecvItemDtype)(int *),
			   DIY_Datatype* (*SendItemDtype)(int *, char**)) {

  // init / clear the items vector
  for (int i = 0; i < assign->NumBlks(); i++) {
    if (i + 1 > (int)items_v.size()) {
      vector<char *> v;
      items_v.push_back(v);
    }
    else
      items_v[i].clear();
  }

  nbhds->ExchangeNeighbors(items_v, wf, RecvItemDtype, SendItemDtype);

  for (int i = 0; i < assign->NumBlks(); i++) {
    num_items[i] = items_v[i].size();
    if (items_v[i].size() > 0)
      items[i] = (void **)(&(items_v[i][0]));
  }

  return 0;

}
//----------------------------------------------------------------------------
//
// flushes exchange with neighbors
//
// items: pointer to received items for each of my blocks [lid][item] (output)
// num_items: number of items for each block (allocated by user)
// RecvItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and
//   creates an MPI datatype for the payloads message
//
// side effects: allocates items and array of pointers to them
//
// returns: error code
//
int DIY_Flush_neighbors(void ***items, int *num_items,
			DIY_Datatype* (*RecvItemDtype)(int *)) {

  // init / clear the items vector
  for (int i = 0; i < assign->NumBlks(); i++) {
    if (i + 1 > (int)items_v.size()) {
      vector<char *> v;
      items_v.push_back(v);
    }
    else
      items_v[i].clear();
  }

  nbhds->FlushNeighbors(items_v, RecvItemDtype);
  for (int i = 0; i < assign->NumBlks(); i++) {
    num_items[i] = items_v[i].size();
    if (items_v[i].size() > 0)
      items[i] = (void **)(&(items_v[i][0]));
  }

  return 0;

}
//----------------------------------------------------------------------------
//
// finds neighbors that intersect bounds +/- extension t
//
// lid: local block id
// bounds: target bounds
// t: additional extension on all sides of bounds
// num_intersect (output) number of intersecting neighbors found
// gids_intersect (output) the intersecting neighbor block gids
//
// returns: error code
//
int DIY_Bounds_intersect_neighbors(int lid, bb_t cell_bounds, float t, 
				   int *num_intersect, int *gids_intersect) {

  nbhds->BoundsIntersectNeighbors(lid, cell_bounds, t, num_intersect, 
				  gids_intersect);

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to neighbors given their global block ids
// reflexive: sends to self block if dest_gids includes global id of self
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// dest_gids: array of gids of neighbors to send to
// num_gids: the number of neighbors to send to
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
int DIY_Enqueue_item_gids(int lid, void *item, int *hdr,
			  int item_size, int *dest_gids, int num_gids,
			  void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_gids; i++)
    nbhds->EnqueueItem(lid, (char *)item, item_size, dest_gids[i], hdr,
		       TransformItem);

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to neighbors given a destination point in each
//  neighbor
// reflexive: sends to self block if points are inside bounds of self
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// dest_pts: points in the destination blocks, by which the destinations can 
// be identified. Points have dimension d and are listed as follows:
//  eg, for dim = 4, x,y,z,t,x,y,z,t,....
// num_dest_pts: number of destination blocks
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
int DIY_Enqueue_item_points(int lid, void *item, int *hdr,
			    int item_size, float *dest_pts, int num_dest_pts,
			    void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_dest_pts; i++) {

    // global id of neighboring block where dest_pt should go
    int neigh_gid = nbhds->Pt2NeighGid(lid, &dest_pts[i * dim]);

    if (neigh_gid >= 0)
      nbhds->EnqueueItem(lid, (char *)item, item_size, neigh_gid, hdr, 
			 TransformItem);

  }

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to one or more neighbors given
//  directionss from the enumeration of possible neighbors. Each direction can
//  be a bitwise OR of several directions
// not reflexive: no direction is defined for sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// neigh_dirs: destination neighbor(s)
// num_neigh_dirs: number of neighbors
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
int DIY_Enqueue_item_dirs(int lid, void *item, int *hdr,
			  int item_size, unsigned char *neigh_dirs, 
			  int num_neigh_dirs,
			  void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_neigh_dirs; i++) {

    nbhds->EnqueueItemDir(lid, (char *)item, item_size, hdr, TransformItem,
			  neigh_dirs[i]);

  }

  return 0;

}
//----------------------------------------------------------------------------
//
// DEPRECATED
//
// Jingyuan's version
//
// enqueues an item for sending to one or more neighbors given a mask array
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// neigh_mask: destination neighbor(s)
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
// int DIY_Enqueue_item_mask(int lid, void *item, int *hdr,
// 			   int item_size, int *neigh_mask,
// 			   void (*TransformItem)(char *, unsigned char)) {

//   nbhds->EnqueueItemMask(lid, (char *)item, item_size, hdr, TransformItem,
// 			 neigh_mask);


//   return 0;

// }
//----------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors
// not reflexive: skips sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// dest_pt: point in the destination block, by which the destination can
//  be identified
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
int DIY_Enqueue_item_all(int lid, void *item, int *hdr, int item_size,
			 void (*TransformItem)(char *, unsigned char)) {

  nbhds->EnqueueItemAll(lid, (char *)item, item_size, hdr, TransformItem,
			true, false);

  return 0;

}
//----------------------------------------------------------------------------
//
// DEPRECATED
//
// enqueues an item for sending to all neighbors that are to one side
//  (eg., left, bottom, rear) of my block
// not reflexive: skips sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// dest_pt: point in the destination block, by which the destination can
//  be identified
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
// int DIY_Enqueue_item_half(int lid, void *item, int *hdr, int item_size,
// 			  void (*TransformItem)(char *, unsigned char)) {

//   nbhds->EnqueueItemAll(lid, (char *)item, item_size, hdr, TransformItem,
// 			false, false);

//   return 0;

// }
//----------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors near enough to receive it
// not reflexive: skips sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// near_pt: point near the destination block
// near_dist: blocks less than or equal to near_dist of the near_pt will be 
//   destinations for the enqueued item . If an item is sent to more than one 
//   neighbor that share faces, it is also sent to the diagonal neighbor 
//   sharing a line or point
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
int DIY_Enqueue_item_all_near(int lid, void *item, int *hdr,
			      int item_size, float *near_pt, float near_dist,
			      void (*TransformItem)(char *, unsigned char)) {

  nbhds->EnqueueItemAllNear(lid, (char *)item, item_size,
			    near_pt, near_dist, hdr, TransformItem, 
			    true);

  return 0;

}
//----------------------------------------------------------------------------
//
// DEPRECATED
//
// enqueues an item for sending to all neighbors near enough to receive it
//  that are to one side (eg., left, bottom, rear) of my block
// not reflexive: skips sending to self block
//
// lid: local id of my block
// item: item to be enqueued
// hdr: pointer to header (or NULL)
// size: size of item in bytes
// near_pt: point near the destination block
// near_dist: blocks less than or equal to near_dist of the near_pt will be 
//   destinations for the enqueued item . If an item is sent to more than one 
//   neighbor that share faces, it is also sent to the diagonal neighbor 
//   sharing a line or point
// TransformItem: pointer to function that transforms the item before
//  enqueueing to a wraparound neighbor, given the wrapping direction
//  (pass NULL if wrapping is unused)
//
// returns: error code
//
// int DIY_Enqueue_item_half_near(int lid, void *item, int *hdr,
// 			       int item_size, float *near_pt, float near_dist,
// 			       void (*TransformItem)(char *, unsigned char)) {

//   nbhds->EnqueueItemAllNear(lid, (char *)item, item_size,
// 			    near_pt, near_dist, hdr, TransformItem, false,
//                          false);

//   return 0;

// }
//----------------------------------------------------------------------------
//
// checks whether all processes are done working via a global reduction
//
// done: whether my local process is done (1 = done, 0 = still working)
//
// returns: whether all processes are done (1 = all done, 0 = still working)
//  (not an error code, unlike most functions)
//
int DIY_Check_done_all(int done) {

  int working = (done ? 0 : 1);
  int sum;
  MPI_Allreduce(&working, &sum, 1, MPI_INT, MPI_SUM, comm);

  return (sum ? 0 : 1);

}
//----------------------------------------------------------------------------
//
// finalizes DIY
//
// returns: error code
//
int DIY_Finalize() {

  delete assign;
  delete io;
  delete merging;
  delete swapping;
  delete blocking;
  BIL_Finalize();
  if (nbhds)
    delete nbhds;

  return 0;

}
//--------------------------------------------------------------------------
//
// creates a vector datatype
//
// num_elems: number of elements in the vector
// stride: number of elements between start of each element (usually 1)
// base_type: data type of vector elements
// type: new (output) data type
//
// returns: error code
//
int DIY_Create_vector_datatype(int num_elems, int stride, 
			       DIY_Datatype base_type, DIY_Datatype *type) {

  MPI_Type_vector(num_elems, 1, stride, base_type, type);
  MPI_Type_commit(type);

  return 0;

}
//--------------------------------------------------------------------------
//
// creates a subarray datatype
//
// num_dims: number of dimensions in the subarray
// full_size: full sizes of the array ([x][y][z] order)
// start_pos: starting indices of the array ([x][y][z] order)
// sub_size: desired sizes of subarray ([x][y][z] order)
// base_type: data type of array elements
// type: new (output) data type
//
// returns: error code
//
int DIY_Create_subarray_datatype(int num_dims, int *full_size, int *sub_size,
				 int *start_pos, DIY_Datatype base_type,
				 DIY_Datatype *type) {

  // fortran order below is not a bug: I always want [x][y][z] order
  MPI_Type_create_subarray(num_dims, full_size, sub_size, start_pos,
			   MPI_ORDER_FORTRAN, base_type, type);
  MPI_Type_commit(type);

  return 0;

}
//--------------------------------------------------------------------------
//
// creates a structure datatype
//
// basse_addr: base address added to relative OFST displacements
// num_map_blocks: number of map blocks
// map: typemap with num_blocks blocks
// type: new (output) datatype
//
// returns: error code
//
int DIY_Create_struct_datatype(DIY_Aint base_addr, int num_map_blocks, 
			struct map_block_t *map, DIY_Datatype *type) {

  vector<map_block_t>tmap(map, map + num_map_blocks);
  CreateDtype(base_addr, &tmap, type);
  MPI_Type_commit(type);

  return 0;

}
//--------------------------------------------------------------------------
//
// destroys a datatype
//
// type: datatype
//
// returns: error code
//
int DIY_Destroy_datatype(DIY_Datatype *type) {

  MPI_Type_free(type);

  return 0;

}
//--------------------------------------------------------------------------
//
// returns an DIY_Aint address given a pointer or address
//
// addr: pointer or address
//
// returns: DIY address (not an error code, unlike most functions)
//
DIY_Aint DIY_Addr(void *addr) {

  MPI_Aint p;
  MPI_Get_address(addr, &p);
  return p;

}
//----------------------------------------------------------------------------
//
// returns the global block identification number (gid)
//   given a local block number
//
// block_num: local block number
//
// returns: global block ID (not an error code, unlike most functions)
//
int DIY_Gid(int block_num) {

  return blocking->Lid2Gid(block_num);

}
//----------------------------------------------------------------------------
//
// block compression
//
// addr: address of start of datatype
// dtype: MPI datatype
// comm: MPI communicator
// comp_buf: pointer to compressed buffer, datatype DIY_BYTE (output)
// comp_size: size of compressed buffer in bytes
//
// side effects: allocates comp_buf, will be freed automatically the next 
// time this function is called, user need not free comp_buf,
// but should use up the result before calling again, because it will disappear
//
// returns: error code
//
int DIY_Compress_block(void* addr, DIY_Datatype dtype, MPI_Comm comm,
		       unsigned char **comp_buf, int *comp_size) {

  comp_buf_v.clear();
  CompressBlock(addr, dtype, comm, &comp_buf_v, comp_size);
  *comp_buf = &comp_buf_v[0];
  return 0;

}
//----------------------------------------------------------------------------
//
// block decompression
//
// in_buf: input block buffer (DIY_BYTE datatype)
// in_size: input size in bytes
// decomp_buf: decompressed buffer
// decomp_size: decompressed size in bytes (output)
//
// side effects: allocates comp_buf, will be freed automatically the next 
// time this function is called, user need not free comp_buf,
// but should use up the result before calling again, because it will disappear
//
// returns: error code
//
int DIY_Decompress_block(unsigned char* in_buf, int in_size, 
			 unsigned char **decomp_buf, int *decomp_size) {

  decomp_buf_v.clear();
  DecompressBlockToBuffer(in_buf, in_size, &decomp_buf_v, decomp_size);
  *decomp_buf = &decomp_buf_v[0];
  return 0;

}
//----------------------------------------------------------------------------
