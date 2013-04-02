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
static int sizes[DIY_MAX_DIM]; // grid size (extents.max - extents.min + 1)
static int tot_blocks = 0; // total number of blocks in all domains
static vector <int> tb; // total number of blocks in each domain
static vector <int> nb; // my local number of blocks in ach domain
static vector <int> maxb; // max. num. of blocks in any process in each domain
static int dim; // number of dimensions
static MPI_Comm comm; // MPI communicator
static int rank; // MPI rank of this process
static int groupsize; // MPI groupsize
static int rounds; // number of merge rounds
static int kvs[DIY_MAX_R]; // merge k-values
static int num_dids = 0; // number of decompostion ids
static int diy_initialized = 0; // whether initialized already

bool dtype_absolute_address = false; // addresses in current datatype
                                     // are absolute w.r.t. MPI_BOTTOM
                                     // or relative w.r.t. base address

vector <Assignment *> assign; // assignment object (abstract base class)
RoundRobinAssignment *round_robin_assign; // round robin assignment object
ProcOrderAssignment *proc_order_assign; // process order assignment object
vector <Blocking *> blocking; // blocking object
Merge *merging; // merge object
Swap *swapping; // swap object
vector <IO *> io; // output I/O object
vector <Neighborhoods *> nbhds; // neighborhoods object
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
// all functions return an error code or a useful value that doubles as en
//  error code such that >= 0 indicates succes and < 0 indicates various errors
// currently few or no functions actually return < 0, indicating an error
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

  if (diy_initialized) // prevent multiple initialiations
    return 0;

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
  merging = new Merge(tot_blocks, comm);
  swapping = new Swap(tot_blocks, comm);

#ifdef OMP
  omp_set_num_threads(num_threads);
#else
  num_threads = num_threads; // quiet compiler warning
#endif

  diy_initialized = 1;
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
// returns: id of this domain (< 0 if error)
//
int DIY_Decompose(int block_order, int glo_num_blocks, int *loc_num_blocks, 
		  int share_face, int *ghost, int *given) {

  Assignment *assign; // assignment object (abstract base class)
  Blocking *blocking; // blocking object
  IO *io; // output I/O object
  Neighborhoods *nbhds; // neighborhoods object
  int nblocks; // local number of blocks in current domain
  int mblocks; // max number of blocks in any process in current domain

  if (block_order == ROUND_ROBIN_ORDER)
    assign = 
      new RoundRobinAssignment(tot_blocks, glo_num_blocks, nblocks, mblocks, 
			       comm);
  else
    assign = 
      new ProcOrderAssignment(tot_blocks, glo_num_blocks, nblocks, mblocks, 
			      comm);
  ::assign.push_back(assign);

  nb.push_back(nblocks);
  *loc_num_blocks = nblocks;
  tb.push_back(glo_num_blocks);
  maxb.push_back(mblocks);

  // decomposing a new domain
  blocking = new Blocking(tot_blocks, num_dids, dim, glo_num_blocks, sizes, 
			  share_face, ghost, given, assign, comm);
  nbhds = new Neighborhoods(num_dids, blocking, assign, comm, false);
  io = new IO(num_dids, dim, glo_num_blocks, mblocks, comm);
  ::blocking.push_back(blocking);
  ::nbhds.push_back(nbhds);
  ::io.push_back(io);

  num_dids++;
  tot_blocks += glo_num_blocks;

  return(num_dids - 1);;

}
//--------------------------------------------------------------------------
//
// Describes the already decomposed domain
//
// loc_num_blocks: local number of blocks on this process
// gids: global ids of my local blocks (unique across all domains)
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
// returns: id of domain (< 0 if error)
//
int DIY_Decomposed(int loc_num_blocks, int *gids, struct bb_t *bounds, 
		   struct ri_t **rem_ids, int *num_rem_ids, int **vids, 
		   int *num_vids, struct gb_t **neighbors, 
		   int *num_neighbors, int wrap) {

  Assignment *assign; // assignment object (abstract base class)
  Blocking *blocking; // blocking object
  IO *io; // output I/O object
  Neighborhoods *nbhds; // neighborhoods object
  int glo_num_blocks; // total number of blocks
  int mblocks; // max number of blocks in any process in current domain

  assign = new ExistingAssignment(tot_blocks, loc_num_blocks, mblocks,
				  glo_num_blocks, comm);
  ::assign.push_back(assign);

  nb.push_back(loc_num_blocks);
  tb.push_back(glo_num_blocks);
  maxb.push_back(mblocks);

  // existing decomposition version of blocking
  blocking = new Blocking(tot_blocks, num_dids, dim, glo_num_blocks, gids, 
			  bounds, assign, comm);
  ::blocking.push_back(blocking);

  // discovery version of constructing the neighborhoods
  nbhds = new Neighborhoods(num_dids, blocking, assign, rem_ids, 
			    num_rem_ids, vids,
			    num_vids, neighbors, num_neighbors, comm, wrap);
  ::nbhds.push_back(nbhds);

  io = new IO(num_dids, dim, glo_num_blocks, mblocks, comm);
  ::io.push_back(io);

  num_dids++;
  tot_blocks += glo_num_blocks;

  return(num_dids - 1);

}
//--------------------------------------------------------------------------
//
// Reads a decomposition from a file
//
// filename: input filename
// swap_bytes: whether to swap bytes for endian conversion
// glo_num__blocks: total number of blocks in the global domain (output)
// loc_num_blocks: local number of blocks on this process (output)
//
// returns: id of this domain (< 0 if error)
//
int DIY_Read_decomposed(char *filename, int swap_bytes,
			int *glo_num_blocks, int *loc_num_blocks) {

  int tot_blocks, nblocks; // global and local number of blocks
  int given[3] = {0, 0, 0}; // no constraints on decomposition in {x, y, z}
  int ghost[6] = {0, 0, 0, 0, 0, 0}; // -x, +x, -y, +y, -z, +z ghost
  MPI_File fd; // file descriptor

  int retval = MPI_File_open(comm, (char *)filename, MPI_MODE_RDONLY,
			     MPI_INFO_NULL, &fd);
  assert(retval == MPI_SUCCESS);

  IO::ReadInfo(fd, swap_bytes, comm, tot_blocks, nblocks);

  MPI_File_close(&fd);

  // only works for regular grids with share_face = 1
  // todo: need a complete solution where  everything needed to create a 
  // decomposition is included in the file
  int did = DIY_Decompose(CONTIGUOUS_ORDER, tot_blocks, &nblocks, 1, 
			  ghost, given);
  *glo_num_blocks = tot_blocks;
  *loc_num_blocks = nblocks;

  return did;

}
//--------------------------------------------------------------------------
//
// block starts and sizes
// for blocks consisting of discrete, regular grid points
//
// did: domain id
// lid: local block id
// starts: pointer to allocated array of starting block extents (output), index
//  of starting grid point (not cell) in each direction  
// sizes: pointer to allocated array of block sizes (output), number of grid
//  points (not cells) in each direction
//
// returns: error code
//
int DIY_Block_starts_sizes(int did, int lid, int *starts, int *sizes) {

  blocking[did]->BlockStartsSizes(lid, starts, sizes);

  return 0;

}
//--------------------------------------------------------------------------
//
// block bounds including ghost
// for blocks consisting of continuous spatiotemporal regions
//
// did: domain id
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
// returns: error code
//
int DIY_Block_bounds(int did, int lid, struct bb_t *bounds) {

  blocking[did]->BlockBounds(lid, bounds);

  return 0;

}
//--------------------------------------------------------------------------
//
// block bounds excluding ghost
// for blocks consisting of continuous spatiotemporal regions
//
// did: domain id
// lid: local block id
// bounds; pointer to a block bounds structure (output), allocated or 
//  declared  by caller
//
// returns: error code
//
int DIY_No_ghost_block_bounds(int did, int lid, struct bb_t *bounds) {

  blocking[did]->NoGhostBlockBounds(lid, bounds);

  return 0;

}
//--------------------------------------------------------------------------
//
// finds the time block to which a local block belongs
//
// did: domain id
// lid: local block id
// time_block: time block containing lid (output)
//
// return: error code
//
int DIY_In_time_block(int did, int lid, int *time_block) {

  int lat_nblocks[DIY_MAX_DIM]; // number of blocks in each dimension
  blocking[did]->NumLatBlocks(lat_nblocks);

  for (int g = 0; g < tb[did]; g++) {
    if (blocking[did]->InTimeBlock(g, lid, sizes[3], lat_nblocks[3])) {
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
// sends an item to a block (asynchronous)
//
// did: domain id
// lid: local block id
// item: item(s) to be sent
// count: number of items
// datatype: item datatype
// dest_gid: destination global block id
//
// returns: error code
//
int DIY_Send(int did, int lid, void *item, int count, DIY_Datatype datatype, 
	     int dest_gid) {

#ifdef _MPI3
  int my_gid = DIY_Gid(did, lid);
  cc->RmaSend(item, count, datatype, my_gid, dest_gid, assign);
#else
  did = did; // quiet compiler warning
  lid = lid; // ditto
  cc->Send(item, count, datatype, dest_gid, assign);
#endif    

  return 0;

}
//--------------------------------------------------------------------------
//
// receives an item from a block (asynchronous)
//
// did: domain id
// lid: local block id
// items: items to be received (output, array af pointers allocated by caller)
// count: number of items received (output)
// wait: whether to wait for one or more items to arrive (0 or 1)
// datatype: item datatype
// src_gids: source global block ids (output, array allocated by caller)
//  only valid if MPI-3 is used, otherwise filled with -1 values
// sizes: size of each item received in datatypes (not bytes)
//  (output, array allocated by caller)
//
// returns: error code
//
int DIY_Recv(int did, int lid, void **items, int *count, int wait,
	     DIY_Datatype datatype, int *src_gids, int *sizes) {

  int my_gid = DIY_Gid(did, lid);

#ifdef _MPI3
  *count = cc->RmaRecv(my_gid, items, datatype, src_gids, wait, assign[did], 
		       sizes);
#else
  *count = cc->Recv(my_gid, items, datatype, wait, sizes);
  for (int i = 0; i < *count; i++)
    src_gids[i] = -1; // only valid for RMA version
#endif

  return 0;

}
//--------------------------------------------------------------------------
//
// flushes asynchronous sending and receiving for all local blocks
//  (collective, must be called by all processes)
//
// barrier: whether to issue a barrier (0 or 1)
//  recommended if more more sends / receives to follow
//
// returns: error code
//
int DIY_Flush_send_recv(int barrier) {

#ifdef _MPI3
  cc->RmaFlushSendRecv(barrier);
#else
  cc->FlushSendRecv(barrier);
#endif

  return 0;

}
//--------------------------------------------------------------------------
//
// configurable in-place merge reduction
//
// did: domain id
// blocks: pointers to input/output blocks, results in first num_blocks_out
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// num_rounds: number of rounds
// k_values: radix (group size) for each round
// reduce_func: pointer to merging or reduction function
// create_func: pointer to function that creates item
// destroy_func: pointer to function the destroys item
// type_func: pointer to function that creates DIY datatype for item 
// num_blocks_out: number of output blocks (output)
//
// side effects: allocates output items and array of pointers to them
//  if not reducing in-place
//
// returns: error code
//
int DIY_Merge_blocks(int did, char **blocks, int **hdrs, int num_rounds, 
		     int *k_values,
		     void (*reduce_func)(char **, int *, int, int *), 
		     char *(*create_func)(int *),
		     void (*destroy_func)(void *),
		     void (*type_func)(void*, DIY_Datatype*, int *),
		     int *num_blocks_out) {

  rounds = num_rounds;
  for (int i = 0; i < num_rounds; i++)
    kvs[i] = k_values[i];

  *num_blocks_out = merging->MergeBlocks(did, blocks, hdrs, num_rounds, 
					 k_values, cc, assign[did], reduce_func,
					 create_func, destroy_func, type_func);

  return 0;

}
//--------------------------------------------------------------------------
//
// configurable in-place swap reduction
//
// did: domain id
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
// send_type_func: pointer to function that creates DIY datatype for sending
//   a subset of the item starting at an element and having a given number
//   of elements (less than the original item)
//   returns the base address associated with the datatype
// recv_type_func: pointer to function that creates DIY datatype for receiving
//   a subset of the item with a given number of elements 
//   (less than the original item)
//
// returns: error code
//
int DIY_Swap_blocks(int did, char **blocks, int **hdrs, int num_elems,
		    int num_rounds, int *k_values, int *starts, int *sizes,
		    void (*reduce_func)(char **, int *, int, int), 
		    char *(*recv_create_func)(int *, int),
		    void (*recv_destroy_func)(void *),
		    void*(*send_type_func)(void*, DIY_Datatype*, int, int),
		    void (*recv_type_func)(void*, DIY_Datatype*, int)) {

  rounds = num_rounds;
  for (int i = 0; i < num_rounds; i++)
    kvs[i] = k_values[i];

  swapping->SwapBlocks(did, blocks, hdrs, num_rounds, k_values, num_elems, 
		       starts, sizes, cc, assign[did], reduce_func,
		       recv_create_func, recv_destroy_func, send_type_func, 
		       recv_type_func);

  return 0;

}
//--------------------------------------------------------------------------
//
// initializes parallel writing of analysis blocks
//
// did: domain id
// filename: output filename
// compress: whether to compress output (0 = normal, 1 = compress)
//   (1: zlib's default compression level 6 is applied blockwise)
//
// returns: eerror code
//
int DIY_Write_open_all(int did, char *filename, int compress) {

  MPI_Barrier(comm); // synchronize to clear any skew before proceeding
  io[did]->WriteAnaInit(filename, compress);

  return 0;

}
//----------------------------------------------------------------------------
//
// writes all analysis blocks in parallel with all other processes
//
// did: domain id
// blocks: array of pointers to analysis blocks
// num_blocks: number of blocks
// hdrs: headers, one per analysis block (NULL if not used)
// type_func: pointer to function that creates DIY datatype for item 
//
// returns: error code
//
int DIY_Write_blocks_all(int did, void **blocks, int num_blocks, int **hdrs,
			 void (*type_func)(void*, int, int, DIY_Datatype*)) {

  int tot_nb_merged = tb[did]; // max possible number of merged blocks
  for (int i = 0; i < rounds; i++)
    tot_nb_merged /= kvs[i];

  // max number of blocks in any process after possible merging
  int max_blocks = maxb[did] < tot_nb_merged ? maxb[did] : tot_nb_merged;

  io[did]->WriteAllAna(blocks, num_blocks, max_blocks, hdrs, type_func);

  return 0;

}
//----------------------------------------------------------------------------
//
// finalizes parallel writing of analysis blocks
//
// did: domain id
//
// returns: error code
//
int DIY_Write_close_all(int did) {

  MPI_Barrier(comm); // flushes all I/O
  io[did]->WriteAnaFinalize();

  return 0;

}
//----------------------------------------------------------------------------
//
// initializes parallel reading of analysis blocks
//
// did: domain id
// filename: input filename
// swap_bytes: whether to swap bytes for endian conversion
//  only applies to reading the headers and footer
//  user must swap bytes manually for datatypes because they are custom
// compress: whether file is compressed (0 = normal, 1 = compressed)
// glo_num__blocks: total number of blocks in the global domain (output)
// loc_num_blocks: local number of blocks on this process (output)
//
// returns: error code
//
int DIY_Read_open_all(int did, char *filename, int swap_bytes, int compress,
		      int *glo_num_blocks, int *loc_num_blocks) {

  int tot_blocks, nblocks; // global and local number of blocks

  MPI_Barrier(comm); // synchronize to clear any skew before proceeding
  io[did]->ReadAnaInit(filename, swap_bytes, compress, tot_blocks, nblocks);
  *glo_num_blocks = tot_blocks;
  *loc_num_blocks = nblocks;

  return 0;

}
//----------------------------------------------------------------------------
//
// reads all analysis blocks in parallel with all other processes
//
// did: domain id
// blocks: pointer to array of pointers for analysis blocks being read (output)
//  DIY will allocate blocks for you
// hdrs: headers, one per analysis block, allocated by caller
//   (pass NULL if not used)
// create_type_func: pointer to a function that allocates a block and
//   creates an DIY datatype for it.
//
// returns: error code
//
int DIY_Read_blocks_all(int did, void ***blocks, int **hdrs,
			void* (*create_type_func)(int, int, int*, 
						  DIY_Datatype*)) {

  io[did]->ReadAllAna(*blocks, hdrs, create_type_func);

  return 0;

}
//----------------------------------------------------------------------------
//
// finalizes parallel reading of analysis blocks
//
// did: domain id
//
// returns: error code
//
int DIY_Read_close_all(int did) {

  MPI_Barrier(comm); // flushes all I/O
  io[did]->ReadAnaFinalize();

  return 0;

}
//----------------------------------------------------------------------------
//
// exchanges items with all neighbors
//
// did: domain id
// items: pointer to received items for each of my blocks [lid][item] (output)
// num_items: number of items for each block (allocated by user)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// suggested value: 0.1
// ItemDtype: pointer to user-supplied function that creates a DIY datatype 
//  for an item to be sent or received
//
// side effects: allocates items and array of pointers to them
//
// returns: error code
//
int DIY_Exchange_neighbors(int did, void ***items, int *num_items, float wf,
			   void (*ItemDtype)(DIY_Datatype *)) {

  // init / clear the items vector
  for (int i = 0; i < assign[did]->NumBlks(); i++) {
    if (i + 1 > (int)items_v.size()) {
      vector<char *> v;
      items_v.push_back(v);
    }
    else
      items_v[i].clear();
  }

  nbhds[did]->ExchangeNeighbors(items_v, wf, ItemDtype);

  for (int i = 0; i < assign[did]->NumBlks(); i++) {
    num_items[i] = (int)items_v[i].size();
    if (items_v[i].size() > 0)
      items[i] = (void **)(&(items_v[i][0]));
  }

  return 0;

}
//----------------------------------------------------------------------------
//
// flushes exchange with neighbors
//
// did: domain id
// items: pointer to received items for each of my blocks [lid][item] (output)
// num_items: number of items for each block (allocated by user)
// ItemDtype: pointer to user-supplied function that creates a DIY datatype 
//  for an item to be sent or received
//
// side effects: allocates items and array of pointers to them
//
// returns: error code
//
int DIY_Flush_neighbors(int did, void ***items, int *num_items,
			void (*ItemDtype)(DIY_Datatype *)) {

  // init / clear the items vector
  for (int i = 0; i < assign[did]->NumBlks(); i++) {
    if (i + 1 > (int)items_v.size()) {
      vector<char *> v;
      items_v.push_back(v);
    }
    else
      items_v[i].clear();
  }

  nbhds[did]->FlushNeighbors(items_v, ItemDtype);
  for (int i = 0; i < assign[did]->NumBlks(); i++) {
    num_items[i] = (int)items_v[i].size();
    if (items_v[i].size() > 0)
      items[i] = (void **)(&(items_v[i][0]));
  }

  return 0;

}
//----------------------------------------------------------------------------
//
// finds neighbors that intersect bounds +/- extension t
//
// did: domain id
// lid: local block id
// bounds: target bounds
// t: additional extension on all sides of bounds
// num_intersect (output) number of intersecting neighbors found
// gids_intersect (output) the intersecting neighbor block gids
//
// returns: error code
//
int DIY_Bounds_intersect_neighbors(int did, int lid, bb_t cell_bounds, float t, 
				   int *num_intersect, int *gids_intersect) {

  nbhds[did]->BoundsIntersectNeighbors(lid, cell_bounds, t, num_intersect, 
				       gids_intersect);

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to neighbors given their global block ids
// reflexive: sends to self block if dest_gids includes global id of self
//
// did: domain id
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
int DIY_Enqueue_item_gids(int did, int lid, void *item, int *hdr,
			  int item_size, int *dest_gids, int num_gids,
			  void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_gids; i++)
    nbhds[did]->EnqueueItem(lid, (char *)item, item_size, dest_gids[i], hdr,
			    TransformItem);

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to neighbors given a destination point in each
//  neighbor
// reflexive: sends to self block if points are inside bounds of self
//
// did: domain id
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
int DIY_Enqueue_item_points(int did, int lid, void *item, int *hdr,
			    int item_size, float *dest_pts, int num_dest_pts,
			    void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_dest_pts; i++) {

    // global id of neighboring block where dest_pt should go
    int neigh_gid = nbhds[did]->Pt2NeighGid(lid, &dest_pts[i * dim]);

    if (neigh_gid >= 0)
      nbhds[did]->EnqueueItem(lid, (char *)item, item_size, neigh_gid, hdr, 
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
// did: domain id
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
int DIY_Enqueue_item_dirs(int did, int lid, void *item, int *hdr,
			  int item_size, unsigned char *neigh_dirs, 
			  int num_neigh_dirs,
			  void (*TransformItem)(char *, unsigned char)) {

  for (int i = 0; i < num_neigh_dirs; i++) {

    nbhds[did]->EnqueueItemDir(lid, (char *)item, item_size, hdr, 
			       TransformItem, neigh_dirs[i]);

  }

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors
// not reflexive: skips sending to self block
//
// did: domain id
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
int DIY_Enqueue_item_all(int did, int lid, void *item, int *hdr, int item_size,
			 void (*TransformItem)(char *, unsigned char)) {

  nbhds[did]->EnqueueItemAll(lid, (char *)item, item_size, hdr, TransformItem,
			     true, false);

  return 0;

}
//----------------------------------------------------------------------------
//
// enqueues an item for sending to all neighbors near enough to receive it
// not reflexive: skips sending to self block
//
// did: domain id
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
int DIY_Enqueue_item_all_near(int did, int lid, void *item, int *hdr,
			      int item_size, float *near_pt, float near_dist,
			      void (*TransformItem)(char *, unsigned char)) {

  nbhds[did]->EnqueueItemAllNear(lid, (char *)item, item_size,
				 near_pt, near_dist, hdr, TransformItem, 
				 true);

  return 0;

}
//----------------------------------------------------------------------------
//
// checks whether all processes are done working via a global reduction
//
// done: whether my local process is done (1 = done, 0 = still working)
//
// returns: whether all processes are done (1 = all done, 0 = still working,
//  < 0 if error)
//
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

  if (!diy_initialized) // prevent finalizing without initializing
    return 0;

  delete merging;
  delete swapping;
  delete cc;

  BIL_Finalize();

  for (int i = 0; i < num_dids; i++) {
    delete assign[i];
    delete io[i];
    delete blocking[i];
    delete nbhds[i];
  }
  if (num_dids) {
    assign.clear();
    blocking.clear();
    nbhds.clear();
    io.clear();
  }

  diy_initialized = 0;
  num_dids = 0;
  tot_blocks = 0;

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
  dtype_absolute_address = false;

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
  dtype_absolute_address = false;

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
  dtype_absolute_address = true;

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
// returns: DIY address (<0 if error)
//
DIY_Aint DIY_Addr(void *addr) {

  MPI_Aint p;
  MPI_Get_address(addr, &p);
  return p;

}
//----------------------------------------------------------------------------
//
// returns the global block identification number (gid)
//   given a local block id
//
// did: domain id
// lid: local block id
//
// returns: global block ID (< 0 if error)
//
int DIY_Gid(int did, int lid) {

  return blocking[did]->Lid2Gid(lid);

}
//----------------------------------------------------------------------------
//
// block compression
//
// addr: address of start of datatype
// dtype: DIY datatype
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
//
// Get total number of domains so far
//
// returns: number of domains (< 0 if error)
//
int DIY_Num_dids() {

  return num_dids;

}
//----------------------------------------------------------------------------
//
// Get total number of blocks in all domains
//
// returns: number of blocks (< 0 if error)
//
int DIY_Tot_num_gids() {

  return tot_blocks;

}
//----------------------------------------------------------------------------
//
// Get global number of blocks in one domain
//
// did: domain id
//
// returns: number of blocks (< 0 if error)
//
int DIY_Num_gids(int did) {

  return tb[did];

}
//----------------------------------------------------------------------------
//
// Get local number of blocks in one domain
//
// did: domain id
//
// returns: number of blocks (< 0 if error)
//
int DIY_Num_lids(int did) {

  return nb[did];

}
//----------------------------------------------------------------------------
//
// Get starting gid in one domain (assuming gids numbered consecutively
//  across domains)
//
// did: domain id
//
// returns: starting gid (< 0 if error)
//
int DIY_Start_gid(int did) {

  return assign[did]->StartGid();

}
//----------------------------------------------------------------------------
//
// Build K-D tree (prototype)
//
// did: domain id
// pts: point locations to be indexed in kd-tree (for now)
// loc_num_pts: global number of points
// glo_num_pts: global number of points
// num_levels: number of tree levels, counting root
// num_bins: number of histogram bins at all levels
//
// returns: error code
//
int DIY_Build_tree(int did, float *pts, int loc_num_pts, int glo_num_pts, 
		   int num_levels, int num_bins) {

  blocking[did]->BuildTree(pts, loc_num_pts, glo_num_pts, num_levels, num_bins);

  return 0;

}
//----------------------------------------------------------------------------
//
// Search K-D tree (prototype)
//
// did: domain id
// pt: target point
// leaf: (output) leaf node containing the point
//
// returns: error code
//
int DIY_Search_tree(int did, float *pt, struct leaf_t *leaf) {

  int gid = blocking[did]->SearchTree(pt, 0);
  blocking[did]->GetLeaf(gid, leaf);

  return 0;

}
//----------------------------------------------------------------------------
