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

extern bool dtype_absolute_address; // addresses in current datatype
                                     // are absolute w.r.t. MPI_BOTTOM
                                     // or relative w.r.t. base address

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

  // init RMA
  MPI_Alloc_mem((DIY_RMA_MAX_ITEMS * 3 + 1) * sizeof(int), MPI_INFO_NULL, 
		&rma_buf);
  for (int i = 0; i < 3 * DIY_RMA_MAX_ITEMS; i++)
    rma_buf[i] = -1; // init items
  rma_buf[DIY_RMA_MAX_ITEMS * 3] = 0; // init number of items
  MPI_Win_create(rma_buf, (DIY_RMA_MAX_ITEMS * 3 + 1) * sizeof(int), 
		 sizeof(int), MPI_INFO_NULL, comm, &rma_win);

}
//----------------------------------------------------------------------------
//
// destructor
//
Comm::~Comm() {

  MPI_Win_free(&rma_win);
  MPI_Free_mem(rma_buf);

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
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::FinishRecvItemsMerge(char **items, int *gids, int *procs,
			       char * (*CreateItem)(int *),
			       void (*RecvItemDtype)(void*, MPI_Datatype*, 
						     int *)) {

  char* (*CreateItemMerge)(int*) = CreateItem;
  char* (*CreateItemSwap)(int*, int) = NULL;
  void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *) = RecvItemDtype;
  void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = NULL;

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
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::RecvItemsMerge(char **items, int *gids, int *procs, float wf,
			 char * (*CreateItem)(int *),
			 void (*RecvItemDtype)(void*, MPI_Datatype*, int *)) {

  char* (*CreateItemMerge)(int*) = CreateItem;
  char* (*CreateItemSwap)(int*, int) = NULL;
  void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *) = RecvItemDtype;
  void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = NULL;

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
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::FinishRecvItemsSwap(char **items, int *gids, int *procs, int ne,
			     char * (*CreateItem)(int *, int),
			     void (*RecvItemDtype)(void*, MPI_Datatype*, 
						    int)) {

  char* (*CreateItemMerge)(int*) = NULL;
  char* (*CreateItemSwap)(int*, int) = CreateItem;
  void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *) = NULL;
  void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = RecvItemDtype;
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
//
// side effects: allocates space for the new items
//
// returns: number of received items
//
int Comm::RecvItemsSwap(char **items, int *gids, int *procs, float wf, int ne,
		      char * (*CreateItem)(int *, int),
		      void (*RecvItemDtype)(void*, MPI_Datatype*, 
					     int)) {

  char* (*CreateItemMerge)(int*) = NULL;
  char* (*CreateItemSwap)(int*, int) = CreateItem;
  void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *) = NULL;
  void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int) = RecvItemDtype;
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
			  void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, 
						     int *),
			  void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, 
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
    MPI_Waitall((int)reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  }

  // flush completion of header receives
  if (use_header) {
    reqs.clear();
    for (int i = 0; i < (int)rcv_hdrs.size(); i++)
      reqs.push_back(rcv_hdrs[i].req);
    MPI_Waitall((int)reqs.size(), &reqs[0], stats);
    for (int i = 0; i < (int)rcv_hdrs.size(); i++)
      rcv_hdrs[i].tag = stats[i].MPI_TAG;
  }

  // post item receives
  reqs.clear();
  for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

    if (RecvItemDtypeMerge) { // merge
      items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
      RecvItemDtypeMerge((void *)items[i], &dm, rcv_hdrs[i].hdr);
      if (dtype_absolute_address)
	addr = MPI_BOTTOM;
      else
	addr = items[i];
    }
    else { // swap
      items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
      RecvItemDtypeSwap((void *)items[i], &dm, ne);
      if (dtype_absolute_address)
	addr = MPI_BOTTOM;
      else
	addr = items[i];
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
  MPI_Waitall((int)reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  snd_items.clear();

  // flush completion of item receives
  MPI_Waitall((int)reqs1.size(), &reqs1[0], stats);
  if (!use_header) {
    for (int i = 0; i < (int)rcv_hdrs.size(); i++) {
      gids[i] = (stats[i].MPI_TAG - 1) / 2;
      procs[i] = stats[i].MPI_SOURCE;
    }
  }

  rcv_hdrs.clear();

  return (int)reqs1.size();

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
		    void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *),
		    void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, 
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
  int min_arr = (int)(wf * (int)rcv_hdrs.size()); // wait for this number of 
                                                  // messsages to arrive
  static bool first = true; // first time

  // post item receives only once if not using headers
  if (!use_header && first) {

    for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

      if (RecvItemDtypeMerge) { // merge
	items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
	RecvItemDtypeMerge((void *)items[i], &dm, rcv_hdrs[i].hdr);
	if (dtype_absolute_address)
	  addr = MPI_BOTTOM;
	else
	  addr = items[i];
      }
      else { // swap
	items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
	RecvItemDtypeSwap((void *)items[i], &dm, ne);
	if (dtype_absolute_address)
	  addr = MPI_BOTTOM;
	else
	  addr = items[i];
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
      MPI_Waitsome((int)reqs.size(), &reqs[0], &narr, arr, stats);
      for (int i = 0; i < (int)rcv_hdrs.size(); i++)
	rcv_hdrs[i].tag = stats[i].MPI_TAG;

      // receive items
      reqs.clear();
      for (int i = 0; i < (int)rcv_hdrs.size(); i++) {

	if (RecvItemDtypeMerge) { // merge
	  items[i] = CreateItemMerge(rcv_hdrs[i].hdr);
	  RecvItemDtypeMerge((void *)items[i], &dm, rcv_hdrs[i].hdr);
	  if (dtype_absolute_address)
	    addr = MPI_BOTTOM;
	  else
	    addr = items[i];
	}
	else { // swap
	  items[i] = CreateItemSwap(rcv_hdrs[i].hdr, ne);
	  RecvItemDtypeSwap((void *)items[i], &dm, ne);
	  if (dtype_absolute_address)
	    addr = MPI_BOTTOM;
	  else
	    addr = items[i];
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

      MPI_Waitsome((int)reqs.size(), &reqs[0], &narr, arr, stats);

    } // no header

  }

  return (int)reqs1.size();

}
//----------------------------------------------------------------------------
//
// sends items (all items in an epoch need to be same datatype)
//
// item: item(s) to be sent
// count: number of items
// datatype: item datatype
// dest_gid: destination gid
// assign: assignment objects
//
void Comm::Send(void *item, int count, DIY_Datatype datatype, 
		int dest_gid, vector <Assignment*> assign) {

  int did; // domain id

  // find the domain of dest_gid
  for (did = 0; did < DIY_Num_dids(); did++) {
    if (DIY_Start_gid(did) <= dest_gid && 
	dest_gid < DIY_Start_gid(did) + DIY_Num_gids(did))
      break;
  }
  assert(did < DIY_Num_dids());
  MPI_Request req;
  int dest_proc = assign[did]->Gid2Proc(dest_gid);
  MPI_Isend(item, count, datatype, dest_proc, dest_gid, comm, &req);
  rma_reqs.push_back(req);

}
//----------------------------------------------------------------------------
//
// receives items (all of the same datatype)
//
// my_gid: my gid
// item: item(s) to be received
// datatype: item datatype
// src_gids: (output) gids of source blocks
// wait: whether to wait for one or more items to arrive (0 or 1)
// sizes: size of each item received in datatypes (not bytes)
//  (output, array allocated by caller)
//
// returns: number of items received
//
int Comm::Recv(int my_gid, void** &items, DIY_Datatype datatype, int wait,
	       int *sizes) {

  MPI_Status status;
  int flag;
  int num_items = 0; // number of items pending for my gid
  int more_items = 0; // maybe more items
  MPI_Aint lb, extent; // datatype lower bound, extent

  // read the RMA buffer, polling as needed according to wait parameter
  while (!num_items || more_items) {

    MPI_Iprobe(MPI_ANY_SOURCE, my_gid, comm, &flag, &status);
    if (flag) {
      int size; // message size in bytes
      MPI_Get_count(&status, MPI_BYTE, &size);
      items[num_items] = new unsigned char[size];
      int src_proc = status.MPI_SOURCE;
      MPI_Recv(items[num_items], size, MPI_BYTE, src_proc, my_gid, comm, 
	       &status);
      MPI_Type_get_extent(datatype, &lb, &extent);
      sizes[num_items] = (int)(size / extent);
      num_items++;
      more_items = 1;
    }
    else
      more_items = 0;

    if (!wait)
      break;

  }

  return num_items;

}
//----------------------------------------------------------------------------
//
// flushes RMA sends and receives
//
// barrier: whether to issue a barrier (0 or 1)
//   (recommended if more sends and receives to follow)
//
void Comm::FlushSendRecv(int barrier) {

  // flush the nonblocking payload sends
  MPI_Status stats[rma_reqs.size()];
  MPI_Waitall((int)rma_reqs.size(), &rma_reqs[0], stats);

  if (barrier)
    MPI_Barrier(comm);

}
//----------------------------------------------------------------------------

#ifdef _MPI3

//----------------------------------------------------------------------------
//
// sends RMA items (all items in an epoch need to be same datatype)
//
// item: item(s) to be sent
// count: number of items
// datatype: item datatype
// my_gid: my gid
// dest_gid: destination gid
// assign: assignment objects
//
void Comm::RmaSend(void *item, int count, DIY_Datatype datatype, 
		   int my_gid, int dest_gid, vector <Assignment*> assign) {

  int did; // domain id
  MPI_Request req;

  // todo: does not check for overflow of receiving RMA window

  int msg[3]; // msg[0] = src gid, msg[1] = dest gid, msg[2] = item count
  msg[0] = my_gid;
  msg[1] = dest_gid;
  msg[2] = count;

  // find the domain of dest_gid
  for (did = 0; did < DIY_Num_dids(); did++) {
    if (DIY_Start_gid(did) <= dest_gid && 
	dest_gid < DIY_Start_gid(did) + DIY_Num_gids(did))
      break;
  }
  assert(did < DIY_Num_dids());
  int dest_proc = assign[did]->Gid2Proc(dest_gid);

  // grab the lock
  MPI_Win_lock(MPI_LOCK_SHARED, dest_proc, 0, rma_win);

  // get and update the count from the RMA buffer
  int num_items; // number of items before incrementing
  int incr = 1; // increment amount
  MPI_Fetch_and_op(&incr, &num_items, MPI_INT, dest_proc, DIY_RMA_MAX_ITEMS * 3,
		   MPI_SUM, rma_win);
  MPI_Win_flush(dest_proc, rma_win);

  // debug
//   fprintf(stderr, "my_gid = %d dest_gid = %d dest_proc = %d num_items = %d\n", 
// 	  my_gid, dest_gid, dest_proc, num_items);

  // put the item
//   MPI_Put(msg, 3, MPI_INT, dest_proc, num_items * 3, 3, MPI_INT, rma_win);
  MPI_Accumulate(msg, 3, MPI_INT, dest_proc, num_items * 3, 3, MPI_INT,
		 MPI_REPLACE, rma_win);
  MPI_Win_flush(dest_proc, rma_win);
  MPI_Win_unlock(dest_proc, rma_win);

  // send the payload
  MPI_Isend(item, count, datatype, dest_proc, dest_gid, comm, &req);
  rma_reqs.push_back(req);

}
//----------------------------------------------------------------------------
//
// receives pending RMA items (all of the same datatype)
//
// my_gid: my gid
// item: item(s) to be received
// datatype: item datatype
// src_gids: (output) gids of source blocks
// wait: whether to wait for one or more items to arrive (0 or 1)
// assign: assignment class object
// sizes: size of each item received in datatypes (not bytes)
//  (output, array allocated by caller)
//
// returns: number of items received
//
int Comm::RmaRecv(int my_gid, void** &items, DIY_Datatype datatype, 
		  int *src_gids, int wait, Assignment *assign, int *sizes) {

  int loc_rma_buf[3 * DIY_RMA_MAX_ITEMS]; // local version of rma buffer
  MPI_Status status;
  int num_items = 0; // number of items pending for my gid

  // read the RMA buffer, polling as needed according to wait parameter
  while (!num_items) {

    MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, rma_win);
    for (int i = 0; i < DIY_RMA_MAX_ITEMS; i++) {
      int msg[3], unused[3];
      MPI_Get_accumulate(unused, 3, MPI_INT, msg, 3, MPI_INT, rank, 3 * i, 3,
			 MPI_INT, MPI_NO_OP, rma_win);
      MPI_Win_flush(rank, rma_win);
      if (msg[0] >= 0 && msg[1] == my_gid && msg[2] >= 0) {
	loc_rma_buf[3 * num_items]     = msg[0];
	loc_rma_buf[3 * num_items + 1] = msg[1];
	loc_rma_buf[3 * num_items + 2] = msg[2];
	msg[0] = -1; // clear the item
	msg[1] = -1;
	msg[2] = -1;
	MPI_Accumulate(msg, 3, MPI_INT, rank, 3 * i, 3, MPI_INT,
		       MPI_REPLACE, rma_win);
	MPI_Win_flush(rank, rma_win);
	num_items++;
      }
    }
    MPI_Win_unlock(rank, rma_win);

    if (!wait)
      break;

    if (!num_items) { // kick the MPI progress engine
      int flag;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		 &flag, MPI_STATUS_IGNORE);
    }

  }

  // process arrived items
  for (int i = 0; i < num_items; i++) {

    int *msg = &loc_rma_buf[3 * i]; // pointer into local rma buffer
    MPI_Aint lb, extent;
    MPI_Type_get_extent(datatype, &lb, &extent);
    items[i] = new unsigned char[msg[2] * extent];
    src_gids[i] = msg[0];
    int src_proc = assign->Gid2Proc(src_gids[i]);
    MPI_Recv(items[i], msg[2], datatype, src_proc, my_gid, comm, &status);
    sizes[i] = msg[2];

  }

  return num_items;

}
//----------------------------------------------------------------------------
//
// flushes RMA sends and receives
//
// barrier: whether to issue a barrier (0 or 1)
//   (recommended if more sends and receives to follow)
//
void Comm::RmaFlushSendRecv(int barrier) {

  // flush the nonblocking payload sends
  MPI_Status stats[rma_reqs.size()];
  MPI_Waitall(rma_reqs.size(), &rma_reqs[0], stats);

  if (barrier)
    MPI_Barrier(comm);

  // clear the RMA buffer, must be locked
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, rma_win);
  for (int i = 0; i < 3 * DIY_RMA_MAX_ITEMS; i++)
    rma_buf[i] = -1; // mark as unused
  rma_buf[DIY_RMA_MAX_ITEMS * 3] = 0; // reset the number of items
  MPI_Win_unlock(rank, rma_win);

  if (barrier)
    MPI_Barrier(comm);

}
//----------------------------------------------------------------------------

#endif

//----------------------------------------------------------------------------
// //
// // DEPRECATED
// //
// // This version should not be needed because all RMA communication has already
// // been flushed, no new data should arrive. Keep it around for a while just
// // in case this is not true.
// //
// // flushes RMA sends and receives
// //
// // item: item(s) to be received
// // datatype: item datatype
// // src_gids: (output) gids of source blocks (allocated by caller)
// // dest_gids: (output) gids of destination blocks (allocated by caller)
// // assign: assignment class object
// //
// // returns: number of items received
// //
// int Comm::RmaFlushSendRecv(void** &items, DIY_Datatype datatype, 
// 			   int *src_gids, int *dest_gids, Assignment *assign) {

//   int loc_rma_buf[3 * DIY_RMA_MAX_ITEMS]; // local version of rma buffer
//   MPI_Status status;
//   int num_items = 0; // number of items pending for my gid

//   // flush the nonblocking payload sends
//   MPI_Status stats[rma_reqs.size()];
//   MPI_Waitall(rma_reqs.size(), &rma_reqs[0], stats);

//   // wait for everyone to be done
//   MPI_Barrier(comm);

//   // read the RMA buffer
//   MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, rma_win);
//   for (int i = 0; i < DIY_RMA_MAX_ITEMS; i++) {
//     int msg[3], unused[3];
//     MPI_Get_accumulate(unused, 3, MPI_INT, msg, 3, MPI_INT, rank, 3 * i, 3,
// 		       MPI_INT, MPI_NO_OP, rma_win);
//     MPI_Win_flush(rank, rma_win);
//     if (msg[0] >= 0 && msg[1] >= 0  && msg[2] >= 0) {
//       loc_rma_buf[3 * num_items]     = msg[0];
//       loc_rma_buf[3 * num_items + 1] = msg[1];
//       loc_rma_buf[3 * num_items + 2] = msg[2];
//       msg[0] = -1; // clear the item
//       msg[1] = -1;
//       msg[2] = -1;
//       MPI_Accumulate(msg, 3, MPI_INT, rank, 3 * i, 3, MPI_INT,
// 		     MPI_REPLACE, rma_win);
//       MPI_Win_flush(rank, rma_win);
//       num_items++;
//     }
//   }
//   MPI_Win_unlock(rank, rma_win);

//   // process arrived items
//   for (int i = 0; i < num_items; i++) {

//     int *msg = &loc_rma_buf[3 * i]; // pointer into local rma buffer
//     MPI_Aint lb, extent;
//     MPI_Type_get_extent(datatype, &lb, &extent);
//     items[i] = new unsigned char[msg[2] * extent];
//     src_gids[i] = msg[0];
//     dest_gids[i] = msg[1];
//     int src_proc = assign->Gid2Proc(src_gids[i]);
//     MPI_Recv(items[i], msg[2], datatype, src_proc, dest_gids[i], comm, &status);

//   }

//   // clear the RMA buffer, must be locked
//   MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, 0, rma_win);
//   for (int i = 0; i < 3 * DIY_RMA_MAX_ITEMS; i++)
//     rma_buf[i] = -1; // mark as unused
//   rma_buf[DIY_RMA_MAX_ITEMS * 3] = 0; // reset the number of items
//   MPI_Win_unlock(rank, rma_win);

//   // everyone's rma buffers are reset and ready for more
//   MPI_Barrier(comm);

//   return num_items;

// }
//----------------------------------------------------------------------------
