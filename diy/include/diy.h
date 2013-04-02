/*---------------------------------------------------------------------------
 *
 * diy C, C++ interface
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * (C) 2011 by Argonne National Laboratory.
 * See COPYRIGHT in top-level directory.
 *
--------------------------------------------------------------------------*/

#ifndef _DIY
#define _DIY

#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <memory.h>
#include "mpi.h"
#ifdef OMP
#include <omp.h>
#endif

/* DIY datatypes aliases for MPI datatypes
   DIY doesn't distinguish between signed and unsigned, lengths are the same
   DIY supports short and long, but their use is not recommended because lengths
   can vary across architectures */
typedef MPI_Datatype DIY_Datatype;
typedef MPI_Aint DIY_Aint;
extern void* DIY_BOTTOM;
extern DIY_Datatype DIY_BYTE; /* 1 byte */
extern DIY_Datatype DIY_SHORT; /* 2 bytes */
extern DIY_Datatype DIY_INT; /* 4 bytes */
extern DIY_Datatype DIY_LONG; /* 4 bytes */
extern DIY_Datatype DIY_LONG_LONG; /* 8 bytes */
extern DIY_Datatype DIY_FLOAT; /* 4 bytes */
extern DIY_Datatype DIY_DOUBLE; /* 8 bytes */
extern DIY_Datatype DIY_LONG_DOUBLE; /* 16 bytes */

/* number of header elements and size of header */
#define DIY_MAX_HDR_ELEMENTS 256

#define DIY_MAX_DIM 10 /* maximum number of dimensions */
#define DIY_MAX_R 64 /* maximum number of rounds */
#define DIY_MAX_VIDS 8 /* maximum number of local vertex ids per block */

/* RMA window size (in items, where 1 item = 3 ints = 12 bytes) */
#define DIY_RMA_MAX_ITEMS 256

/* block order */
#define ROUND_ROBIN_ORDER 0
#define CONTIGUOUS_ORDER 1

/* displacement types */
#define OFST 0
#define ADDR 1

/* neighbor direction enumeration
   used to identify direction of one neighbor block for both regular and 
   wrapround neighbors
   can be bitwise ORed, eg., maximum-side neighboring block in 3 dimensions 
   would be DIY_X1 | DIY_Y1 | DIY_Z1 
   each use identifies exactly one block
   eg. DIY_X0 is the (one) left neighbor block, not the entire left plane */
#define DIY_X0    0x01 /* minimum-side x (left) neighbor */
#define DIY_X1    0x02 /* maximum-side x (right) neighbor */
#define DIY_Y0    0x04 /* minimum-side y (bottom) neighbor */
#define DIY_Y1    0x08 /* maximum-side y (top) neighbor */
#define DIY_Z0    0x10 /* minimum-side z (back) neighbor */
#define DIY_Z1    0x20 /* maximum-side z (front)neighbor */
#define DIY_T0    0x30 /* minimum-side t (earlier) neighbor */
#define DIY_T1    0x40 /* maximum-side t (later) neighbor */

/* block bounds */
struct bb_t { 
  float min[DIY_MAX_DIM];
  float max[DIY_MAX_DIM];
};

/* one global block */
struct gb_t {
  int gid; /* global block id */
  int vids[DIY_MAX_VIDS]; /* optional local ids for vertices in the block on
			     the owner process, used to discover neighbors */
  int num_vids; /* number of vids */
  int proc; /* process to which block is assigned */
  unsigned char neigh_dir; /* neighbor direction, if this block is used as a
			      neighbor and its direction from the local block
			      is needed */
  struct bb_t bb; /* block bounds */
};

/* remote ids used for aligning neighboring blocks during their discovery */
struct ri_t {
  int vid; /* vertex id on remote process */
  int proc; /* remote process */
};

/* typemap block for creating custom datatypes */
struct map_block_t {
  DIY_Datatype base_type; /* existing datatype used to create this one */
  int disp_type; /* diplacement is relative OFST or absolute ADDR */
  int count;  /* count of each element */
  DIY_Aint disp; /* displacement of each element in bytes
		 OFSTs are from the start of the type and ADDRs are from 0x */
};

/* leaf node from a tree */
struct leaf_t {
  int gid; /* node global id */
  struct bb_t bounds; /* bounds of region correponding to this tree node */
  int proc; /* process to which this node is block is or will be assigned */
};

/* ----------------------------------------------------------------------- */

/*
  Public API

  all functions return an error code or a useful value that doubles as en
  error code such that >= 0 indicates succes and < 0 indicates various errors
  currently few or no functions actually return < 0, indicating an error
  todo: add error checking
*/

/* ----------------------------------------------------------------------- */

/*
  Initializes DIY

  dim: number of dimensions (2, 3, or 4) (input)
  data_size: data size in each dimension (only for structured data), pass NULL
    for other cases
  num_threads: number of threads DIY is allowed to use (>= 1)
    ignored if DIY is built with openmp disabled (enabled by default)
  comm: MPI cpmmunicator

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Init(int dim, int *data_size, int num_threads,
	     MPI_Comm comm);

/* ----------------------------------------------------------------------- */

/*
  Decomposes the domain

  block_order: ROUND_ROBIN_ORDER or CONTIGUOUS_ORDER numbering of
  global block ids to processes (input)
  glo_num__blocks: total number of blocks in the global domain (input)
    pass 0 or anything for CONTIGUOUS_ORDER (unused)
  loc_num_blocks: local number of blocks on this process (output)
  share_face: whether neighboring blocks share a common face or are
  separated by a gap of one unit
  ghost: ghost layer for each dimension and side (min, max)
   each entry can be 0 (no ghost) or > 0 (this many ghost cells per side)
   {x min side ghost, x max side ghost, y min side ghost, y max side ghost...}
  given: constraints on the blocking entered as an array where
  0 implies no constraint in that direction and some value n > 0 is a given
  number of blocks in a given direction
  eg., {0, 0, 0, t} would result in t blocks in the 4th dimension

  returns: id of domain (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Decompose(int block_order, int glo_num_blocks, int *loc_num_blocks, 
		  int share_face, int *ghost, int *given);

/* ----------------------------------------------------------------------- */

/*
  Describes the already decomposed domain

  loc_num_blocks: local number of blocks on this process
  gids: global ids of my local blocks (unique across all domains)
  bounds: block bounds (extents) of my local blocks
  rem_ids: remote ids used for neighbor discovery (pass NULL if 
  neighbor discovery not needed at all)
  num_rem_ids: number of remote ids for each local block (pass 0 for any
  blocks that don't need neighbor discovery)
  vids: local vertex ids for each local block that needs neighbor discovery (pass
  NULL if not needed at all)
  num_vids: number of vids for each local block (pass 0 for any blocks that
  don't need neighbor discovery)
  neighbors: neighbor lists for each of my local blocks, in lid order
  lists can be partial, containing only known information, as long as 
  the rem_data lists contain enough information for DIY to discover the rest
  of the unknown neighbors. Unknown (remote) block bounds can be 
  uninitialized. If wrapping is used, (see below) a neighbor direction
  must be provided for each neighbor
  num_neighbors: number of neighbors for each of my local blocks, in lid order
  wrap: whether wraparound neighbors are used (0 = no wraparound neighbors
  were provided, 1 = provided some wraparound neighbors)

  returns: id of domain (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Decomposed(int loc_num_blocks, int *gids, struct bb_t *bounds, 
		   struct ri_t **rem_ids,
		   int *num_rem_ids, int **vids, int *num_vids,
		   struct gb_t **neighbors, int *num_neighbors, int wrap);

/* ----------------------------------------------------------------------- */
/*
  Initializes parallel reading of analysis blocks and decomposes domain
  accordingliy

  filename: input filename
  swap_bytes: whether to swap bytes for endian conversion
  glo_num__blocks: total number of blocks in the global domain (output)
  loc_num_blocks: local number of blocks on this process (output)

  returns: id of this domain (< 0 if error)
*/
#ifdef __cplusplus
extern "C"
#endif
int DIY_Read_decomposed(char *filename, int swap_bytes,
			int *glo_num_blocks, int *loc_num_blocks);

/* ----------------------------------------------------------------------- */
/*
  Finds block starts and sizes
  for blocks consisting of discrete, regular grid points

  did: domain id
  lid: local block id
  starts: pointer to allocated array of starting block extents (output), index
  of starting grid point (not cell) in each direction
  sizes: pointer to allocated array of block sizes (output), number of grid
  points (not cells) in each direction

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Block_starts_sizes(int did, int lid, int *starts, int *sizes);

/* ----------------------------------------------------------------------- */

/*
  Finds block bounds including ghost
  for blocks consisting of continuous spatiotemporal regions

  did: domain id
  lid: local block id
  bounds; pointer to a block bounds structure (output), allocated or 
  declared  by caller

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Block_bounds(int did, int lid, struct bb_t *bounds);

/* ----------------------------------------------------------------------- */
/*
  Finds block bounds excluding ghost
  for blocks consisting of continuous spatiotemporal regions

  did: domain id
  lid: local block id
  bounds; pointer to a block bounds structure (output), allocated or 
  declared  by caller

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_No_ghost_block_bounds(int did, int lid, struct bb_t *bounds);

/* ----------------------------------------------------------------------- */

/*
  Finds the time block to which a local block belongs

  did: domain id
  lid: local block id
  time_block: time block containing lid (output)

  return: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_In_time_block(int did, int lid, int *time_block);

/* ----------------------------------------------------------------------- */

/*
  Posts a data read for a block

  block_starts: starting grid indices for minimum corner of block
  block_sizes: size of block in each dimension (number of vertices)
  file_name: input file name
  var_type: input datatype
  buffer: pointer to void * data buffer

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Add_data_raw(int *block_starts, int *block_sizes, char *file_name, 
		     DIY_Datatype var_type, void** buffer);

/* ----------------------------------------------------------------------- */

/*
  Executes a parallel data read of all blocks

  note: performs an MPI_Barrier afterwards to eliminate any process
  skew before proceeding

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Read_data_all();

/* ----------------------------------------------------------------------- */

/*
  sends an item to a block (asynchronous)

  did: domain id
  lid: local block id
  item: item(s) to be sent
  count: number of items
  datatype: item datatype
  dest_gid: destination global block id

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Send(int did, int lid, void *item, int count, DIY_Datatype datatype, 
	     int dest_gid);

/* ----------------------------------------------------------------------- */

/*
  receives an item from a block (asynchronous)

  did: domain id
  lid: local block id
  items: items to be received (output, array of pointers allocated by caller)
  count: number of items received (output)
  wait: whether or not to wait for one or more items to arrive (0 or 1)
  datatype: item datatype
  src_gids: source global block ids (output, array allocated by caller)
  only valid if MPI-3 is used, otherwise filled with -1 values
  sizes: size of each item received in datatypes (not bytes)
   (output, array allocated by caller)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Recv(int did, int lid, void **items, int *count, int wait,
	     DIY_Datatype datatype, int *src_gids, int *sizes);

/* ----------------------------------------------------------------------- */

/*
  flushes asynchronous sending and receiving
  (collective, must be called by all processes)

  barrier: whether to issue a barrier (0 or 1)
  recommended if more sends and receives to follow

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Flush_send_recv(int barrier);

/* ----------------------------------------------------------------------- */
/* DEPRECATED */
/* /\* */
/*   flushes asynchronous sending and receiving */
/*   (collective, must be called by all processes) */

/*   items: items to be received (output, array of pointers allocated by caller) */
/*   count: number of items received (output) */
/*   datatype: item datatype */
/*   src_gids: source global block ids (output, array allocated by caller) */
/*   dest_gids: my destination global block ids (output, array allocated by caller) */

/*   returns: error code */
/* *\/ */

/* #ifdef __cplusplus */
/* extern "C" */
/* #endif */
/* int DIY_Flush_send_recv(void **items, int *count, DIY_Datatype datatype,  */
/* 			int *src_gids, int *dest_gids); */

/* /\* ----------------------------------------------------------------------- *\/ */

/*
  Configurable in-place merge reduction

  did: domain id
  blocks: pointers to input/output blocks, result in in_blocks[0]
  hdrs: pointers to input headers (optional, pass NULL if unnecessary)
  num_rounds: number of rounds
  k_values: radix (group size) for each round
  reduce_func: pointer to merging or reduction function
  create_func: pointer to function that creates item
  destroy_func: pointer to function that destroys item
  type_func: pointer to function that creates DIY datatype for item 
  num_blocks_out: number of output blocks (output)

  side effects: allocates output items and array of pointers to them
    if not reducing in-place

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Merge_blocks(int did, char **blocks, int **hdrs, int num_rounds, 
		     int *k_values,
		     void (*reduce_func)(char **, int *, int, int *), 
		     char *(*create_func)(int *),
		     void (*destroy_func)(void *),
		     void (*type_func)(void *, DIY_Datatype*, int *),
		     int *num_blocks_out);

/* ----------------------------------------------------------------------- */

/*
  Configurable in-place swap reduction

  did: domain id
  blocks: pointers to input/output blocks
  hdrs: pointers to input headers (optional, pass NULL if unnecessary)
  num_elems: number of elements in a block
  num_rounds: number of rounds
  k_values: radix (group size) for each round
  starts: start of result in each block (output)
  sizes: number of elements in result in each block (output)
  starts and sizes are allocated by the caller
  reduce_func: pointer to reduction function
  recv_create_func: pointer to function that creates received item 
  with given number of elements (less than original item)
  recv_destroy_func: pointer to function that desroys received item
  part of a total number of parts in the item
  send_type_func: pointer to function that creates DIY datatype for sending
  a subset of the item starting at an element and having a given number
  of elements (less than the original item)
  returns the base address associated with the datatype
  recv_type_func: pointer to function that creates DIY datatype for receiving
  a subset of the item with a given number of elements 
  (less than the original item)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Swap_blocks(int did, char **blocks, int **hdrs, int num_elems,
		    int num_rounds, int *k_values, int *starts, int *sizes,
		    void (*reduce_func)(char **, int *, int, int), 
		    char *(*recv_create_func)(int *, int),
		    void (*recv_destroy_func)(void *),
		    void*(*send_type_func)(void*, DIY_Datatype*, int, int),
		    void (*recv_type_func)(void*, DIY_Datatype*, int));

/* ----------------------------------------------------------------------- */

/*
  Initializes parallel writing of analysis blocks

  did: domain id
  filename: output filename
  compress: whether to compress output (0 = normal, 1 = compress)
  (1: zlib's default compression level 6 is applied blockwise)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Write_open_all(int did, char *filename, int compress);

/* ----------------------------------------------------------------------- */

/*
  Writes all analysis blocks in parallel with all other processes

  did: domain id
  blocks: array of pointers to analysis blocks
  num_blocks: number of blocks
  hdrs: headers, one per analysis block (NULL if not used)
  type_func: pointer to function that creates DIY datatype for item 

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Write_blocks_all(int did, void **blocks, int num_blocks, int **hdrs,
			 void (*type_func)(void*, int, int, DIY_Datatype*));

/* ----------------------------------------------------------------------- */

/*
  Finalizes parallel writing of analysis blocks

  did: domain id

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Write_close_all(int did);

/* ----------------------------------------------------------------------- */

/*
  Initializes parallel reading of analysis blocks

  did: domain id
  filename: output filename
  swap_bytes: whether to swap bytes for endian conversion
  only applies to reading the headers and footer
  user must swap bytes manually for datatypes because they are custom
  compress: whether file is compressed (0 = normal, 1 = compressed)
  glo_num__blocks: total number of blocks in the global domain (output)
  loc_num_blocks: local number of blocks on this process (output)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Read_open_all(int did, char *filename, int swap_bytes, int compress,
		      int *glo_num_blocks, int *loc_num_blocks);

/* ----------------------------------------------------------------------- */

/*
  Reads all analysis blocks in parallel with all other processes

  did: domain id
  blocks: pointer to array of pointers for analysis blocks being read (output)
  DIY will allocate blocks for you
  hdrs: headers, one per analysis block, allocated by caller
  (pass NULL if not used)
  create_type_func: pointer to a function that allocates a block and
  creates a DIY datatype for it.

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Read_blocks_all(int did, void ***blocks, int **hdrs,
			void* (*create_type_func)(int, int, int *, 
						  DIY_Datatype*));

/* ----------------------------------------------------------------------- */

/*
  Finalizes parallel reading of analysis blocks

  did: domain id

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Read_close_all(int did);

/* ----------------------------------------------------------------------- */

/*
  Exchanges items with all neighbors

  did: domain id
  items: pointer to received items for each of my blocks [lid][item] (output)
  num_items: number of items for each block (allocated by user)
  wf: wait_factor for nonblocking communication [0.0-1.0]
  0.0 waits the minimum (1 message per round)
  1.0 waits the maximum (all messages per round)
  suggested value: 0.1
  ItemDtype: pointer to user-supplied function that creates a DIY datatype 
  for an item to be sent or received

  side effects: allocates items and array of pointers to them

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Exchange_neighbors(int did, void ***items, int *num_items, float wf,
			   void (*ItemDtype)(DIY_Datatype *));

/* ----------------------------------------------------------------------- */

/*
  Flushes exchange with neighbors

  did: domain id
  items: pointer to received items for each of my blocks [lid][item] (output)
  num_items: number of items for each block (allocated by user)
  ItemDtype: pointer to user-supplied function that creates a DIY datatype 
  for an item to be sent or received

  side effects: allocates items and array of pointers to them

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Flush_neighbors(int did, void ***items, int *num_items,
			void (*ItemDtype)(DIY_Datatype *));

/* ----------------------------------------------------------------------- */

/*
  finds neighbors that intersect bounds +/- extension t

  did: domain id
  lid: local block id
  bounds: target bounds
  t: additional extension on all sides of bounds
  num_intersect (output) number of intersecting neighbors found
  gids_intersect (output) the intersecting neighbor block gids

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Bounds_intersect_neighbors(int did, int lid, struct bb_t cell_bounds, 
				   float t, 
				   int *num_intersect, int *gids_intersect);

/* ----------------------------------------------------------------------- */

/*
  Enqueues an item for sending to neighbors given their global block ids
  Reflexive: sends to self block if dest_gids includes global id of self

  did: domain id
  lid: local id of my block
  item: item to be enqueued
  hdr: pointer to header (or NULL)
  size: size of item in bytes
  dest_gids: array of gids of neighbors to send to
  num_gids: the number of neighbors to send to
  TransformItem: pointer to function that transforms the item before
  enqueueing to a wraparound neighbor, given the wrapping direction
  (pass NULL if wrapping is unused)

  returns: error code

*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Enqueue_item_gids(int did, int lid, void *item, int *hdr,
			  int item_size, int *dest_gids, int num_gids,
			  void (*TransformItem)(char *, unsigned char));

/* ----------------------------------------------------------------------- */

/*
  Enqueues an item for sending to a neighbor given a destination point in each
  neighbor
  Reflexive: sends to self block if points are inside bounds of self

  did: domain id
  lid: local id of my block
  item: item to be enqueued
  hdr: pointer to header (or NULL)
  size: size of item in bytes
  dest_pts: points in the destination blocks, by which the destinations can 
  be identified. Points have dimension d and are listed as follows:
  eg, for dim = 4, x,y,z,t,x,y,z,t,....
  num_dest_pts: number of destination points
  TransformItem: pointer to function that transforms the item before
  enqueueing to a wraparound neighbor, given the wrapping direction
  (pass NULL if wrapping is unused)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Enqueue_item_points(int did, int lid, void *item, int *hdr,
			    int item_size, float *dest_pts, int num_dest_pts,
			    void (*TransformItem)(char *, unsigned char));

/* ----------------------------------------------------------------------- */

/*
  Enqueues an item for sending to one or more neighbors given
  directions from the enumeration of possible neighbors. Each direction can
  be a bitwise OR of several directions
  Not reflexive: no direction is defined for sending to self block

  did: domain id
  lid: local id of my block
  item: item to be enqueued
  hdr: pointer to header (or NULL)
  size: size of item in bytes
  neigh_dirs: destination neighbor(s)
  num_neigh_dirs: number of neighbors
  TransformItem: pointer to function that transforms the item before
  enqueueing to a wraparound neighbor, given the wrapping direction
  (pass NULL if wrapping is unused)


  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Enqueue_item_dirs(int did, int lid, void *item, int *hdr,
			  int item_size, unsigned char *neigh_dirs, 
			  int num_neigh_dirs,
			  void (*TransformItem)(char *, unsigned char));

/* ----------------------------------------------------------------------- */

/*
  Enqueues an item for sending to all neighbors
  Not reflexive: skips sending to self block

  did: domain id
  lid: local id of my block
  item: item to be enqueued
  hdr: pointer to header (or NULL)
  size: size of item in bytes
  dest_pt: point in the destination block, by which the destination can
  be identified
  TransformItem: pointer to function that transforms the item before
  enqueueing to a wraparound neighbor, given the wrapping direction
  (pass NULL if wrapping is unused)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Enqueue_item_all(int did, int lid, void *item, int *hdr, int item_size,
			 void (*TransformItem)(char *, unsigned char));

/* ----------------------------------------------------------------------- */

/*
  Enqueues an item for sending to all neighbors near enough to receive it
  Not reflexive: skips sending to self block

  did: domain id
  lid: local id of my block
  item: item to be enqueued
  hdr: pointer to header (or NULL)
  size: size of item in bytes
  near_pt: point near the destination block
  near_dist: blocks less than or equal to near_dist of the near_pt will be 
  destinations for the enqueued item . If an item is sent to more than one 
  neighbor that share faces, it is also sent to the diagonal neighbor 
  sharing a line or point
  TransformItem: pointer to function that transforms the item before
  enqueueing to a wraparound neighbor, given the wrapping direction
  (pass NULL if wrapping is unused)

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Enqueue_item_all_near(int did, int lid, void *item, int *hdr, 
			      int item_size, float *near_pt, float near_dist,
			      void (*TransformItem)(char *, unsigned char));


/* ----------------------------------------------------------------------- */

/*
  checks whether all processes are done working via a global reduction

  done: whether my local process is done (1 = done, 0 = still working)

  returns: whether all processes are done (1 = all done, 0 = still working,
  < 0 if error)

*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Check_done_all(int done);

/* ----------------------------------------------------------------------- */

/*
  Finalizes DIY

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Finalize();

/* ----------------------------------------------------------------------- */

/*
  Creates a vector datatype

  num_elems: number of elements in the vector
  stride: number of elements between start of each element (usually 1)
  base_type: data type of vector elements
  type: new (output) data type

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Create_vector_datatype(int num_elems, int stride, 
			       DIY_Datatype base_type, DIY_Datatype *type);

/* ----------------------------------------------------------------------- */

/*
  Creates a subarray datatype

  num_dims: number of dimensions in the subarray
  full_size: full sizes of the array ([x][y][z] order)
  start_pos: starting indices of the array ([x][y][z] order)
  sub_size: desired sizes of subarray ([x][y][z] order)
  base_type: data type of array elements
  type: new (output) data type

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Create_subarray_datatype(int num_dims, int *full_size, int *sub_size,
				 int *start_pos, DIY_Datatype base_type,
				 DIY_Datatype *type);

/* ----------------------------------------------------------------------- */

/*
  Creates a structure datatype

  basse_addr: base address added to relative OFST displacements
  num_map_blocks: number of map blocks
  map: typemap with num_blocks blocks
  type: new (output) datatype

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Create_struct_datatype(DIY_Aint base_addr, int num_blocks, 
			       struct map_block_t *map, DIY_Datatype *type);

/* ----------------------------------------------------------------------- */

/*
  Destroys a datatype

  type: datatype

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Destroy_datatype(DIY_Datatype *type);

/* ----------------------------------------------------------------------- */

/*
  Returns an DIY_Aint address given a pointer or address

  addr: pointer or address

  returns: DIY address (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
DIY_Aint DIY_Addr(void *addr);

/* ----------------------------------------------------------------------- */

/*
  Returns the global block identification number (gid) given a local 
  block number

  did: domain id
  lid: local block id

  returns: global block ID (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Gid(int did, int lid);

/* ----------------------------------------------------------------------- */

/*
  Compresses a block

  addr: address of start of datatype
  dtype: DIY datatype
  comm: MPI communicator
  comp_buf: pointer to compressed buffer, datatype DIY_BYTE (output)
  comp_size: size of compressed buffer in bytes

  side effects: allocates comp_buf, will be freed automatically the next 
  time this function is called, user need not free comp_buf,
  but should use up the result before calling again, because it will disappear

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Compress_block(void* addr, DIY_Datatype dtype, MPI_Comm comm,
		       unsigned char **comp_buf, int *comp_size);

/* ----------------------------------------------------------------------- */

/*
  Decompresses a block

  in_buf: input block buffer (DIY_BYTE datatype)
  in_size: input size in bytes
  decomp_buf: decompressed buffer
  decomp_size: decompressed size in bytes (output)

  side effects: allocates comp_buf, will be freed automatically the next 
  time this function is called, user need not free comp_buf,
  but should use up the result before calling again, because it will disappear

  returns: error code
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Decompress_block(unsigned char* in_buf, int in_size, 
			 unsigned char **decomp_buf, int *decomp_size);

/* ----------------------------------------------------------------------- */

/* 
   Build Kd-tree (prototype)

   not documented until I have a better idea how this will be used 
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Build_tree(int did, float *pts, int loc_num_pts, int glo_num_pts,
		   int num_levels,int num_bins);

/* ----------------------------------------------------------------------- */

/* 
   Search Kd-tree (prototype)

   not documented until I have a better idea how this will be used 
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Search_tree(int did, float *pt, struct leaf_t *leaf);

/* ----------------------------------------------------------------------- */

/*
  Get total number of domains so far

  returns: number of domains (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Num_dids();

/* ----------------------------------------------------------------------- */

/*
  Get total number of blocks in all domains

  returns: number of blocks (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Tot_num_gids();

/* ----------------------------------------------------------------------- */

/*
  Get global number of blocks in one domain

  did: domain id

  returns: number of blocks (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Num_gids(int did);

/* ----------------------------------------------------------------------- */

/*
  Get local number of blocks in one domain

  did: domain id

  returns: number of blocks (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Num_lids(int did);

/* ----------------------------------------------------------------------- */

/*
  Get starting gid in one domain (assuming gids numbered consecutively
  across domains)

  did: domain id

  returns: starting gid (< 0 if error)
*/

#ifdef __cplusplus
extern "C"
#endif
int DIY_Start_gid(int did);

/* ----------------------------------------------------------------------- */

#endif
