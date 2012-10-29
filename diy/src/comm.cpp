//------------------------------------------------------------------------------
//
// communication class
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

#include "comm.hpp"

//--------------------------------------------------------------------------
//
// constructor
//
// comm: MPI communicator
//
Comm::Comm(MPI_Comm comm) {

  this->comm = comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

}
//----------------------------------------------------------------------------
//
// sends one item to a destination rank
// item: item to be sent
//  specifically, the base address of the displacements in the item's datatype
//  usually the address of the item but could be MPI_BOTTOM if absolute
//  addresses are used in the datatype
// hdr: header containing quantity information
// dest_rank: destination rank
// gid: source or destination gid, depending on how caller wants to use it
//  message gets tagged with gid and gid is retrieved upon receipt
// dtype: MPI datatype of block
//
void Comm::SendItem(char *item, int *hdr, int dest_rank, int gid,
		  MPI_Datatype *dtype) {

  MPI_Request req;
  snd_item si; // sent item

  si.hdr = hdr;
  if (hdr) {
    MPI_Isend(hdr, DIY_MAX_HDR_ELEMENTS, MPI_INT, dest_rank, 2 * gid, comm, 
	      &req);
    si.hreq = req;
    use_header = true;
  }
  else
    use_header = false;
  MPI_Isend(item, 1, *dtype, dest_rank, 2 * gid + 1, comm, &req);

  si.ireq = req;
  snd_items.push_back(si);

}
//----------------------------------------------------------------------------
//
// initiates receiving one item from a source rank
// item is not available until FinishRecvItems is called
//
// src_rank: source rank
// use_header: whether headers are used
//
// returns: pointer to the item
//
void Comm::StartRecvItem(int src_rank, bool use_header) {

  int *hdr = NULL;
  MPI_Request req;
  rcv_hdr rh; // received header

  if (use_header) {
    hdr = new int[DIY_MAX_HDR_ELEMENTS];
    MPI_Irecv(hdr, DIY_MAX_HDR_ELEMENTS, MPI_INT, src_rank, MPI_ANY_TAG, 
	      comm, &req);
    rh.hdr = hdr;
    rh.req = req;
    this->use_header = true;
  }
  else {
    rh.hdr = NULL;
    this->use_header = false;
  }
  rh.proc = src_rank;
  rh.tag = -1;
  rcv_hdrs.push_back(rh);

}
//----------------------------------------------------------------------------
//
// completes item communication for merge
//
// items: array of pointers to items (output)
// gids: gids of received messeages
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// CreateItem: pointer to user-supplied function that takes a header
//   and creates an item, returning a pointer to it
// RecvItemDtype: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::FinishRecvItemsMerge(char **items, int *gids, int *procs,
			      char * (*CreateItem)(int *),
			      void* (*RecvItemDtype)(void*, MPI_Datatype*)) {

  char* (*CreateItemMerge)(int*) = CreateItem;
  char* (*CreateItemSwap)(int*, int) = NULL;
  void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*) = RecvItemDtype;
  void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = NULL;

  return FinishRecvItems(items, gids, procs, 0, CreateItemMerge, 
			 CreateItemSwap, RecvItemDtypeMerge, RecvItemDtypeSwap);

}
//----------------------------------------------------------------------------
//
// item communication for merge
//
// items: array of pointers to items (output)
// gids: gids of received messeages
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// CreateItem: pointer to user-supplied function that takes a header
//   and creates an item, returning a pointer to it
// RecvItemDtype: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::RecvItemsMerge(char **items, int *gids, int *procs, float wf,
			      char * (*CreateItem)(int *),
			      void* (*RecvItemDtype)(void*, MPI_Datatype*)) {

  char* (*CreateItemMerge)(int*) = CreateItem;
  char* (*CreateItemSwap)(int*, int) = NULL;
  void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*) = RecvItemDtype;
  void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = NULL;

  return RecvItems(items, gids, procs, wf, 0, CreateItemMerge, 
		   CreateItemSwap, RecvItemDtypeMerge, RecvItemDtypeSwap);

}
//----------------------------------------------------------------------------
//
// completes item communication for swap
//
// items: array of pointers to items (output)
// gids: gids of received messeages (output)
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// ne: number of elements in the received item (less than entire item)
// CreateItem: pointer to user-supplied function that takes a header and
//   the current part out of total number of parts, and
//   creates an item, returning a pointer to it
// RecvItemDtype: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::FinishRecvItemsSwap(char **items, int *gids, int *procs, int ne,
			     char * (*CreateItem)(int *, int),
			     void* (*RecvItemDtype)(void*, MPI_Datatype*, 
						    int)) {

  char* (*CreateItemMerge)(int*) = NULL;
  char* (*CreateItemSwap)(int*, int) = CreateItem;
  void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*) = NULL;
  void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = RecvItemDtype;
  return FinishRecvItems(items, gids, procs, ne, CreateItemMerge, 
			 CreateItemSwap, RecvItemDtypeMerge, RecvItemDtypeSwap);

}
//----------------------------------------------------------------------------
//
// item communication for swap
//
// items: array of pointers to items (output)
// gids: gids of received messeages (output)
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// wf: wait_factor for nonblocking communication [0.0-1.0]
// 0.0 waits the minimum (1 message per round)
// 1.0 waits the maximum (all messages per round)
// ne: number of elements in the received item (less than entire item)
// CreateItem: pointer to user-supplied function that takes a header and
//   the current part out of total number of parts, and
//   creates an item, returning a pointer to it
// RecvItemDtype: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::RecvItemsSwap(char **items, int *gids, int *procs, float wf, int ne,
		      char * (*CreateItem)(int *, int),
		      void* (*RecvItemDtype)(void*, MPI_Datatype*, 
					     int)) {

  char* (*CreateItemMerge)(int*) = NULL;
  char* (*CreateItemSwap)(int*, int) = CreateItem;
  void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*) = NULL;
  void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = RecvItemDtype;
  return RecvItems(items, gids, procs, wf, ne, CreateItemMerge, 
		   CreateItemSwap, RecvItemDtypeMerge, RecvItemDtypeSwap);

}
//----------------------------------------------------------------------------
//
// completes item communication
//
// items: array of pointers to items (output)
// gids: gids of received messeages (output)
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// ne: number of elements in the received item (less than entire item)
//   used only for swap; pass 0 (or whatever) for mege
// CreateItemMerge: pointer to user-supplied function that takes a header
//   and creates an item, returning a pointer to it
// CreateItemSwap: similar to CreateItemMerge, but function takes two
//  more ints for fraction of total data
// RecvItemDtypeMerge: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
// RecvItemDtypeSwap: similar to RecvItemDtypeMerge, but function takes two
//  more ints for fraction of total data to convert to data type
// Only one of RecvItemDtypeMerge and RecvItemDtypeSwap should be non-null,
//  the other null
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::FinishRecvItems(char **items, int *gids, int *procs, int ne,
			 char* (*CreateItemMerge)(int *),
			 char* (*CreateItemSwap)(int *, int),
			 void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*),
			 void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, 
						    int)) {

  MPI_Request req;
  vector<MPI_Request> reqs, reqs1;
  void *addr; // base address for datatype
  MPI_Datatype dm;
  MPI_Status stats[rcv_hdrs.size()];

  // flush completion of header sends
  if (use_header) {
    for (int i = 0; i < (int)snd_items.size(); i++)
      reqs.push_back(snd_items[i].hreq);
    MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  }

  // flush completion of header receives
  if (use_header) {
    reqs.clear();
    for (int i = 0; i < (int)rcv_hdrs.size(); i++)
      reqs.push_back(rcv_hdrs[i].req);
    MPI_Waitall(reqs.size(), &reqs[0], stats);
    for (int i = 0; i < (int)rcv_hdrs.size(); i++)
      rcv_hdrs[i].tag = stats[i].MPI_TAG;
  }

  // post item receives
  reqs.clear();
  for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

    if (RecvItemDtypeMerge) { // merge
      items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
      addr = RecvItemDtypeMerge((void *)items[i], &dm);
    }
    else { // swap
      items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
      addr = RecvItemDtypeSwap((void *)items[i], &dm, ne);
    }
    if (use_header) {	
      MPI_Irecv(addr, 1, dm, rcv_hdrs[i].proc, rcv_hdrs[i].tag + 1, comm, &req);
      gids[i] = rcv_hdrs[i].tag / 2;
      procs[i] = rcv_hdrs[i].proc;
      delete[] rcv_hdrs[i].hdr;
    }
    else
      MPI_Irecv(addr, 1, dm, rcv_hdrs[i].proc, MPI_ANY_TAG, comm, &req);
    MPI_Type_free(&dm);
    reqs1.push_back(req);

  }

  // flush completion of item sends
  reqs.clear();
  for (int i = 0; i < (int)snd_items.size(); i++)
    reqs.push_back(snd_items[i].ireq);
  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  snd_items.clear();

  // flush completion of item receives
  MPI_Waitall(reqs1.size(), &reqs1[0], stats);
  if (!use_header) {
    for (int i = 0; i < (int)rcv_hdrs.size(); i++) {
      gids[i] = (stats[i].MPI_TAG - 1) / 2;
      procs[i] = stats[i].MPI_SOURCE;
    }
  }

  rcv_hdrs.clear();

  return reqs1.size();

}
//----------------------------------------------------------------------------
//
// item communication
//
// items: array of pointers to items (output)
// gids: gids of received messeages (output)
//   source or destination gid, depending on how caller used it when sending
// procs: sending processes of received messages (output)
// ne: number of elements in the received item (less than entire item)
//   used only for swap; pass 0 (or whatever) for mege
// CreateItemMerge: pointer to user-supplied function that takes a header
//   and creates an item, returning a pointer to it
// CreateItemSwap: similar to CreateItemMerge, but function takes two
//  more ints for fraction of total data
// RecvItemDtypeMerge: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//   and returns the base address associated with the datatype
// RecvItemDtypeSwap: similar to RecvItemDtypeMerge, but function takes two
//  more ints for fraction of total data to convert to data type
// Only one of RecvItemDtypeMerge and RecvItemDtypeSwap should be non-null,
//  the other null
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::RecvItems(char **items, int *gids, int *procs, float wf, int ne,
		    char* (*CreateItemMerge)(int *),
		    char* (*CreateItemSwap)(int *, int),
		    void* (*RecvItemDtypeMerge)(void*, MPI_Datatype*),
		    void* (*RecvItemDtypeSwap)(void*, MPI_Datatype*, 
					       int)) {

  MPI_Request req;
  vector<MPI_Request> reqs, reqs1;
  void *addr; // base address for datatype
  MPI_Datatype dm;
  MPI_Status stats[rcv_hdrs.size()];
  MPI_Status stat;
  int arr[rcv_hdrs.size()]; // requests that arrived
  int narr; // number of requests that arrived
  int tot_narr = 0; // total number counts-receive messages arrived this round
  int min_arr = (int)(wf * rcv_hdrs.size()); // wait for this number of 
                                             // messsages to arrive
  static bool first = true; // first time

  // post item receives only once if not using headers
  if (!use_header && first) {

    for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

      if (RecvItemDtypeMerge) { // merge
	items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
	addr = RecvItemDtypeMerge((void *)items[i], &dm);
      }
      else { // swap
	items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
	addr = RecvItemDtypeSwap((void *)items[i], &dm, ne);
      }
      MPI_Irecv(addr, 1, dm, rcv_hdrs[i].proc, MPI_ANY_TAG, comm, &req);
      MPI_Type_free(&dm);
      reqs1.push_back(req);

    }

    first = false;

  }

  while (tot_narr < min_arr) {

    // using header
    if (use_header) {

      // wait for enough headers to arrive
      reqs.clear();
      for (int i = 0; i < (int)rcv_hdrs.size(); i++)
	reqs.push_back(rcv_hdrs[i].req);
      MPI_Waitsome(reqs.size(), &reqs[0], &narr, arr, stats);
      for (int i = 0; i < (int)rcv_hdrs.size(); i++)
	rcv_hdrs[i].tag = stats[i].MPI_TAG;

      // receive items
      reqs.clear();
      for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

	if (RecvItemDtypeMerge) { // merge
	  items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
	  addr = RecvItemDtypeMerge((void *)items[i], &dm);
	}
	else { // swap
	  items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
	  addr = RecvItemDtypeSwap((void *)items[i], &dm, ne);
	}
	MPI_Recv(addr, 1, dm, rcv_hdrs[i].proc, rcv_hdrs[i].tag + 1, comm, 
		 &stat);
	gids[i] = rcv_hdrs[i].tag / 2;
	procs[i] = rcv_hdrs[i].proc;
	delete[] rcv_hdrs[i].hdr;
	MPI_Type_free(&dm);
	reqs1.push_back(req);

      }

    } // using header

    // no header
    else {

      MPI_Waitsome(reqs.size(), &reqs[0], &narr, arr, stats);

    } // no header

  }

  return reqs1.size();

}
//----------------------------------------------------------------------------
