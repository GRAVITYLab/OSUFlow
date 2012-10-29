//---------------------------------------------------------------------------
//
// merge class
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

#include "merge.hpp"

//----------------------------------------------------------------------------
//
// constructor
//
// comm: MPI communicator
//
Merge::Merge(MPI_Comm comm) {

  this->comm = comm;

}
//----------------------------------------------------------------------------
//
// gets the global ids of the blocks in my group
//
// inputs:
// kv: vector of k values
// cur_r: current round number (0 to r - 1)
// gid: global id of the block
// partners: global ids of the blocks in my group, including myself (output)
//
// returns: true: my gid is the root (maximum) of the group, else false
//
inline bool Merge::GetPartners(const int *kv, int cur_r, 
			       int gid, int *partners) {

  int step = 1; // gids jump by this much in the current round
  int p; // position of this gid in the group
  int r, k;

  for (r = 0; r < cur_r; r++)
    step *= kv[r];

  p = gid / step % kv[cur_r];
  partners[0] = gid - p * step;
  for (k = 1; k < kv[cur_r]; k++)
    partners[k] = partners[k - 1] + step;

  return(gid == partners[kv[cur_r] - 1]);

}
//----------------------------------------------------------------------------
//
// radix-k merge
//
// its: pointers to input/ouput items, results in first number of output items
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// nr: number of rounds
// kv: k vector, radix for each round
// cc: pointer to communicate object
// assign: pointer to assignment object
// merge_func: pointer to merging function
// create_func: pointer to function that creates item
// destroy_func: pointer to function that destroys item
// type_func: pointer to function that creates MPI datatype for item 
//   returns the base address associated with the datatype
//
// side effects: allocates output items and array of pointers to them, if
//   not reducing in-place
//
// returns: number of output items
//
int Merge::MergeBlocks(char **its, int **hdrs, 
		       int nr, int *kv, Comm *cc, Assignment *assign,
		       void (*merge_func)(char **, int *, int), 
		       char * (*create_func)(int *),
		       void (*destroy_func)(void *),
		       void* (*type_func)(void*, MPI_Datatype*)) {

  int rank, groupsize; // MPI usual
  int gid; // global id of current item block
  int p; // process rank
  void *addr; // base address of datatype
  MPI_Datatype dtype; // data type
  int ng; // number of groups this process owns
  int nb = assign->NumBlks(); // number of blocks this process owns
  vector<char *> my_its(its, its + nb);  // copy of its
  vector<bool> done(nb, false); // done items
  vector<int> root_gids; // distinct gids of root blocks

  // init
  assert(nr > 0 && nr <= DIY_MAX_R); // sanity
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // for all rounds
  for (int r = 0; r < nr; r++){

    int n_recv = 0; // number of received blocks by root
    int partners[kv[r]]; // gids of partners in a group
    root_gids.clear();
    root_gids.reserve(kv[r]);

    // all my current blocks must participate in a round (send or receive)
    for (int b = 0; b < nb; b++) {

      if (!done[b]) { // blocks that survived to this round

	gid = assign->Lid2Gid(b);
	bool root = GetPartners(kv, r, gid, partners);

	if (!root) { // nonroots post sends of headers and items
	  p = assign->Gid2Proc(partners[kv[r] - 1]);
	  addr = type_func(my_its[b], &dtype);
	  // tag is source block gid
	  if (hdrs)
	    cc->SendItem((char *)addr, hdrs[b], p, gid, &dtype);
	  else
	    cc->SendItem((char *)addr, NULL, p, gid, &dtype);
	  MPI_Type_free(&dtype);
	  done[b] = true; // nonroot blocks are done after they have been sent
	}

	else { // root posts receives of headers
	  root_gids.push_back(partners[kv[r] - 1]);
	  for (int k = 0; k < kv[r] - 1; k++) { // receive the others
	    p = assign->Gid2Proc(partners[k]);
	    cc->StartRecvItem(p, hdrs);
	    n_recv++;
	  }
	}

      } // blocks that survived to this round

    } // all my current blocks

    // finish receiving all items
    char *recv_its[n_recv]; // received items
    int recv_gids[n_recv]; // (source) gids of the received items
    int recv_procs[n_recv]; // source proc of each received item
    cc->FinishRecvItemsMerge(recv_its, recv_gids, recv_procs, create_func, 
			     type_func); 

    // merge each group
    ng = root_gids.size(); // number of groups this process owns
    for (int j = 0; j < ng; j++) {

      vector<char *>reduce_its; // items ready for reduction in a group
      vector<int>reduce_gids; // gids for reduce_its
      reduce_its.reserve(kv[r]);
      reduce_gids.reserve(kv[r]);

      int lid = assign->Gid2Lid(root_gids[j]);
      reduce_its.push_back(my_its[lid]);
      reduce_gids.push_back(root_gids[j]);

      GetPartners(kv, r, root_gids[j], partners);

      for (int i = 0; i < n_recv; i++) { // collect items for this group
	if (find(partners, partners + kv[r], recv_gids[i]) != 
	    partners + kv[r]) {
	  reduce_its.push_back(recv_its[i]);
	  reduce_gids.push_back(recv_gids[i]);
	}
      }

      merge_func(&reduce_its[0], &reduce_gids[0], kv[r]);
      my_its[lid] = reduce_its[0];

    }

    // cleanup
    if (ng) {
      for (int i = 0; i < n_recv; i++)
	destroy_func(recv_its[i]);
    }

  } // for all rounds

  // move results to the front, swapping them rather than copying so that user
  // can free all items without having duplicated pointers that get freed
  // multiple times
  for (int i = 0; i < ng; i++) {
    char *temp = its[i];
    its[i] = my_its[assign->Gid2Lid(root_gids[i])];
    its[assign->Gid2Lid(root_gids[i])] = temp;
  }

  return ng;

}
//----------------------------------------------------------------------------
//
// asynchronous radix-k merge
//
// its: pointers to input/ouput items, results in first number of output items
// hdrs: pointers to input headers (optional, pass NULL if unnecessary)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// nr: number of rounds
// kv: k vector, radix for each round
// cc: pointer to communicate object
// assign: pointer to assignment object
// merge_func: pointer to merging function
// create_func: pointer to function that creates item
// destroy_func: pointer to function that destroys item
// type_func: pointer to function that creates MPI datatype for item 
//   returns the base address associated with the datatype
//
// side effects: allocates output items and array of pointers to them, if
//   not reducing in-place
//
// returns: number of output items
//
int Merge::AsyncMergeBlocks(char **its, int **hdrs, float wf,
			    int nr, int *kv, Comm *cc, Assignment *assign,
			    void (*merge_func)(char **, int *, int), 
			    char * (*create_func)(int *),
			    void (*destroy_func)(void *),
			    void* (*type_func)(void*, MPI_Datatype*)) {

  int rank, groupsize; // MPI usual
  int gid; // global id of current item block
  int p; // process rank
  void *addr; // base address of datatype
  MPI_Datatype dtype; // data type
  int ng = 0; // number of groups this process owns
  int nb = assign->NumBlks(); // number of blocks this process owns
  vector<char *> my_its(its, its + nb);  // copy of its
  vector<bool> done(nb, false); // done items
  vector<int> root_gids; // distinct gids of root blocks

  // init
  assert(nr > 0 && nr <= DIY_MAX_R); // sanity
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // for all rounds
  for (int r = 0; r < nr; r++){

    int n_recv = 0; // number of received blocks by root
    int partners[kv[r]]; // gids of partners in a group
    root_gids.clear();
    root_gids.reserve(kv[r]);

    // all my current blocks must participate in a round (send or receive)
    for (int b = 0; b < nb; b++) {

      if (!done[b]) { // blocks that survived to this round

	gid = assign->Lid2Gid(b);
	bool root = GetPartners(kv, r, gid, partners);

	if (!root) { // nonroots post sends of headers and items
	  p = assign->Gid2Proc(partners[kv[r] - 1]);
	  addr = type_func(my_its[b], &dtype);
	  // tag is source block gid
	  if (hdrs)
	    cc->SendItem((char *)addr, hdrs[b], p, gid, &dtype);
	  else
	    cc->SendItem((char *)addr, NULL, p, gid, &dtype);
	  MPI_Type_free(&dtype);
	  done[b] = true; // nonroot blocks are done after they have been sent
	}

	else { // root posts receives of headers
	  root_gids.push_back(partners[kv[r] - 1]);
	  for (int k = 0; k < kv[r] - 1; k++) { // receive the others
	    p = assign->Gid2Proc(partners[k]);
	    cc->StartRecvItem(p, hdrs);
	    n_recv++;
	  }
	}

      } // blocks that survived to this round

    } // all my current blocks


    // get and merge one or more items at a time

    wf = (wf < 0.1 ? 0.1 : wf); // clamp minimum wf
    int num_merge_rounds = 1.0 / wf;
    for (int mr = 0; mr < num_merge_rounds; mr++) {

      // finish receiving all items
      char *recv_its[n_recv]; // received items
      int recv_gids[n_recv]; // (source) gids of the received items
      int recv_procs[n_recv]; // source proc of each received item

      if (mr < num_merge_rounds - 1)
	cc->RecvItemsMerge(recv_its, recv_gids, recv_procs, wf, create_func, 
				 type_func); 
      else
	cc->RecvItemsMerge(recv_its, recv_gids, recv_procs, wf, create_func, 
			   type_func); 

      // merge each group
      ng = root_gids.size(); // number of groups this process owns
      for (int j = 0; j < ng; j++) {

	vector<char *>reduce_its; // items ready for reduction in a group
	vector<int>reduce_gids; // gids for reduce_its
	reduce_its.reserve(kv[r]);
	reduce_gids.reserve(kv[r]);

	int lid = assign->Gid2Lid(root_gids[j]);
	reduce_its.push_back(my_its[lid]);
	reduce_gids.push_back(root_gids[j]);

	GetPartners(kv, r, root_gids[j], partners);

	for (int i = 0; i < n_recv; i++) { // collect items for this group
	  if (find(partners, partners + kv[r], recv_gids[i]) != 
	      partners + kv[r]) {
	    reduce_its.push_back(recv_its[i]);
	    reduce_gids.push_back(recv_gids[i]);
	  }
	}

	merge_func(&reduce_its[0], &reduce_gids[0], kv[r]);
	my_its[lid] = reduce_its[0];

      }

      // cleanup
      if (ng) {
	for (int i = 0; i < n_recv; i++)
	  destroy_func(recv_its[i]);
      }

    }

  } // for all rounds

  // move results to the front, swapping them rather than copying so that user
  // can free all items without having duplicated pointers that get freed
  // multiple times
  for (int i = 0; i < ng; i++) {
    char *temp = its[i];
    its[i] = my_its[assign->Gid2Lid(root_gids[i])];
    its[assign->Gid2Lid(root_gids[i])] = temp;
  }

  return ng;

}
//----------------------------------------------------------------------------
