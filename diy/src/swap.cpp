//---------------------------------------------------------------------------
//
// swap class
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

#include "swap.hpp"

//----------------------------------------------------------------------------
//
// constructor
//
// start_b: starting block global id (number of blocks in prior domains)
// comm: MPI communicator
//
Swap::Swap(int start_b, MPI_Comm comm) {

  this->start_b = start_b;
  this->comm = comm;

}
//----------------------------------------------------------------------------
//
// computes group number and position within that group for my block
// to participate in the swap communication
//
// group number is 0 to the global number of groups in the current round - 1
// position number is 0 to k value of the current round - 1
//
// inputs:
// cur_r: current round
// kv: vector of k values
// gid: global id of the block
//
// outputs:
// g: group number
// p: position number within the group
//
inline void Swap::GetGrpPos(int cur_r, const int *kv, int gid, 
			     int &g, int &p) {

  int step = 1;

  for (int i = 0; i < cur_r; i++)
    step *= kv[i];

  // the second term in the following expression does not simplify to
  // (gid - start_b) / kv[r]
  // because the division gid / (step * kv[r]) is integer and truncates
  // this is exactly what we want
  g = (gid - start_b) % step + (gid - start_b) / (step * kv[cur_r]) * step;

  p = (gid - start_b) / step % kv[cur_r];

}
//----------------------------------------------------------------------------
//
// gets the global ids of the blocks in my group
//
// inputs:
// kv: vector of k values
// cur_r: current round number (0 to r - 1)
// gid: global id of the block
//
// output:
// partners: global ids of the partners (blocks) in my group, including myself
// grp: number of the group to which I belong
// pos: my position in that group
//
inline void Swap::GetPartners(const int *kv, int cur_r, int gid, 
			      int *partners, int &grp, int &pos) {

  int step = 1; // gids jump by this much in the current round

  GetGrpPos(cur_r, kv, gid, grp, pos);

  for (int r = 0; r < cur_r; r++)
    step *= kv[r];

  partners[0] = gid - pos * step;
  for (int k = 1; k < kv[cur_r]; k++)
    partners[k] = partners[k - 1] + step;

}
//----------------------------------------------------------------------------
//
// radix-k swap
//
// did: decomposition id
// its: pointers to input/output items (reduced in-place)
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// nr: number of rounds
// kv: k vector, radix for each round
// num_elems: number of elements in the item
// starts: start of result in each block (output)
// sizes: number of elements in result in each block (output)
//  starts and sizes are allocated by the caller
// cc: pointer to communicate object
// assign: pointer to assignment object
// reduce_func: pointer to reduction function
// recv_create_func: pointer to function that creates received item
// recv_destroy_func: pointer to function that destroys received item
//   with given number of elements (less than original item)
// send_type_func: pointer to function that creates MPI datatype for sending
//   a subset of the item starting at an element and having a given number
//   of elements (less than the original item)
// recv_type_func: pointer to function that creates MPI datatype for receiving
//   a subset of the item with a given number of elements 
//   (less than the original item)
//
void Swap::SwapBlocks(int did, char **its, int **hdrs, int nr, int *kv, 
		      int num_elems, int *starts, int *sizes,
		      Comm *cc, Assignment *assign,
		      void (*reduce_func)(char **, int *, int, int), 
		      char* (*recv_create_func)(int *, int),
		      void (*recv_destroy_func)(void *),
		      void* (*send_type_func)(void*, MPI_Datatype*, int, int),
		      void (*recv_type_func)(void*, MPI_Datatype*, int)) {

  int p; // process rank
  int nb = assign->NumBlks(); // number of local blocks
  void *addr; // base address of datatype
  MPI_Datatype dtype;
  int grp; // group number
  int pos; // position in group
  int ne, se; // number of elements and starting element of current part of item

  // init
  assert(nr > 0 && nr <= DIY_MAX_R); // sanity

  // for all rounds
  for (int r = 0; r < nr; r++) {

    int send_gids[nb * kv[r]];
    int num_send_gids = 0;    
    int partners[kv[r]]; // gids in my group, including myself

    // all my blocks
    for (int b = 0; b < nb; b++) {

      int gid = DIY_Gid(did, b);
      GetPartners(kv, r, gid, partners, grp, pos);

      // init start and size of this block's active part
      if (r == 0)  {
	starts[b] = 0;
	sizes[b] = num_elems;
      }

      // for all communicating blocks in this group
      for (int i = 0; i < kv[r]; i++) {

	// starting element and number of elements in the subset of the item
	se = i * sizes[b] / kv[r];
	if (i == kv[r] - 1)
	  ne = sizes[b] - (i * sizes[b] / kv[r]);
	else
	  ne = sizes[b] / kv[r]; // number of elements in this part

	// post sends of headers and items
	p = assign->Gid2Proc(partners[i]);
	addr = send_type_func(its[b], &dtype, se, ne);
	send_gids[num_send_gids] = gid;
	// tag is destination block gid
	if (hdrs)
	  cc->SendItem((char *)addr, hdrs[b], p, partners[i], &dtype);
	else
	  cc->SendItem((char *)addr, NULL, p, partners[i], &dtype);
	num_send_gids++;
	MPI_Type_free(&dtype);

	// post receives of headers
	if (i == pos) { // receiver for this part
	  for (int k = 0; k < kv[r]; k++) { // sources
	    p = assign->Gid2Proc(partners[k]); // todo: change this
	    // tag is my block gid
	    if (hdrs)
	      cc->StartRecvItem(p, true);
	    else
	      cc->StartRecvItem(p, false);
	  } // sources
	} // receiver

      } // for all blocks in this group

      // update start and size of this block's active part for next time
      starts[b] += (pos * sizes[b] / kv[r]);
      if (pos == kv[r] - 1)
	sizes[b] = sizes[b] - (pos * sizes[b] / kv[r]);
      else
	sizes[b] = sizes[b] / kv[r];

    } // for all my blocks

    // finish receiving and reduce blocks
    ReduceBlocks(did, its, r, kv, sizes, cc, assign, reduce_func, 
		 recv_create_func, recv_destroy_func, recv_type_func);

  } // for all rounds

}
//----------------------------------------------------------------------------
//
// finish communication and reduce blocks
//
// did: decomposition id
// its: pointers to input and output items (reduced in-place)
// cur_r: current round
// kv: k-values for all rounds
// sz_part: size of active part in each of my blocks
// cc: pointer to communicate object
// assign: pointer to assignment object
// reduce_func: pointer to reduction function
// recv_create_func: pointer to function that creates item for the current
//  part out of a total number of parts in the item
// recv_destroy_func: pointer to function that destroys received item
// recv_type_func: pointer to function that creates MPI datatype for the current
//   part of a total of parts in theitem, and
//
void Swap::ReduceBlocks(int did, char** its, int cur_r, int *kv,
			int *sz_part, Comm *cc, Assignment *assign,
			void (*reduce_func)(char **, int *, int, int), 
			char* (*recv_create_func)(int *, int),
			void (*recv_destroy_func)(void *),
			void (*recv_type_func)(void*, MPI_Datatype*, int)) {

  int nb = assign->NumBlks();
  int n_recv = nb * kv[cur_r];
  int partners[kv[cur_r]];

  // received items
  char *recv_its[n_recv]; // received items
  int recv_gids[n_recv]; // gids of desitination blocks doing the reduction
  int recv_procs[n_recv]; // source proc of each received item

  // max size of active part of a block
  int max_sz_part = -1;
  for (int b = 0; b < nb; b++) {
    if (sz_part[b] > max_sz_part)
      max_sz_part = sz_part[b];
  }

  // finish receiving all items w.t.r my process (for all my blocks)
  cc->FinishRecvItemsSwap(recv_its, recv_gids, recv_procs, max_sz_part, 
			  recv_create_func, recv_type_func); 

  // for my local blocks
  for (int b = 0; b < nb; b++) {

    // pointers to items for reduction, reduce_its[0] is the original (full) 
    // item, others are received (partial) items
    vector<char *>reduce_its;
    reduce_its.reserve(kv[cur_r]);

    // source gids of items for reduction
    vector<int> reduce_gids;
    reduce_gids.reserve(kv[cur_r]);

    // get partners in the group of my current block
    int unused;
    int pos; // position of my block in the group
    int gid = DIY_Gid(did, b);
    GetPartners(kv, cur_r, gid, partners, unused, pos);

    // collect items to be reduced
    for (int j = 0; j < kv[cur_r]; j++) {

      reduce_gids.push_back(partners[j]);

      // find the received item that matches the gid
      // note the received gids are those of my current block, but we need to
      // derive the sending gids from the sending proc knowing that multiple
      // blocks from the same process will arrive in the order they were sent
      for (int i = 0; i < n_recv; i++) {
	if (gid == recv_gids[i] && 
	    assign->Gid2Proc(partners[j]) == recv_procs[i]) {
	  reduce_its.push_back(recv_its[i]);
	  recv_procs[i] = -1; // mark as done; will not match again
	}
      }

    }

    if (pos != 0) { // swap my block to be first in the group
      int temp_gid = reduce_gids[0];
      reduce_gids[0] = reduce_gids[pos];
      reduce_gids[pos] = temp_gid;
      char *temp_it = reduce_its[0];
      reduce_its[0] = reduce_its[pos];
      reduce_its[pos] = temp_it;
    }

    // reduce each block and (shallow) copy out
    reduce_func(&reduce_its[0], &reduce_gids[0], kv[cur_r], sz_part[b]);
    its[b] = reduce_its[0];

    // cleanup
    for (int i = 1; i < kv[cur_r]; i++)
      recv_destroy_func(reduce_its[i]);

  } // for all my local blocks

}
//----------------------------------------------------------------------------
