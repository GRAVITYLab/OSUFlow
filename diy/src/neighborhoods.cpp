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

//
// constructs and initializes my neighborhoods
//
// blocking: pointer to blocking class
// assignment: pointer to assignment class
// comm: MPI commnicator
// nhdr: optional number of header counts
//
Neighborhoods::Neighborhoods(Blocking *blocking, Assignment *assignment, 
			     MPI_Comm comm, int nhdr) {

  this->comm = comm;
  this->blocking = blocking;
  this->assign = assignment;
  this->nhdr = nhdr;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // init my blocks list
  for (int i = 0; i < assign->NumBlks(); i++)  {
    bl_t block;
    vector<struct gb_t> gbs;
    block.gid = assign->RoundRobin_lid2gid(i);
    blocking->GetNeighbors(i, gbs);
    for (vector<struct gb_t>::iterator gi = gbs.begin(); gi != gbs.end(); 
	 gi++) {
      ne_t neighbor;
      neighbor.gb = *gi;
      block.neighbors.push_back(neighbor);
    }
    blocks.push_back(block);
  }

  tot_items_sent = 0;
  tag = 0;

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
// enqueues an item for sending to a neighbor
//
// lid: local id of my block
// item: item to be enqueued (char * pointer can point to anything, does not
// need to be chars)
// size: size of item in bytes
// neigh_gid: global id of neighboring destination block
// hdr: if nhdr > 0, pointer tp header counts
// (additional quantities of subitems within the item)
//
void Neighborhoods::EnqueueItem(int lid, char *item, size_t size, 
				int neigh_gid, int *hdr) {

  // todo: provide a non-copy mode for very large items

  int neighbor; // local neighbor number of neigh_gid

  // copy item so that it persists
  char *p = new char[size];
  memcpy(p, item, size);

  // find the neighbor gid in my neighboring blocks
  for (neighbor = 0; neighbor < (int)blocks[lid].neighbors.size(); neighbor++) {
    if (neigh_gid == blocks[lid].neighbors[neighbor].gb.gid)
      break;
  }
  assert(neighbor < (int)blocks[lid].neighbors.size()); // sanity

  blocks[lid].neighbors[neighbor].items.push_back(p);
  if (blocks[lid].neighbors[neighbor].gb.proc != rank)
    tot_items_sent++;

  int j = blocks[lid].neighbors[neighbor].items.size();
  if ((int)blocks[lid].neighbors[neighbor].hdr.size() < j)
    blocks[lid].neighbors[neighbor].hdr.resize(j);
  for (int i = 0; i < nhdr; i++) 
    blocks[lid].neighbors[neighbor].hdr[j].push_back(hdr[i]);

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
// RecvItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message  and
//   creates an MPI datatype for the payloads message
// SendItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and a payloads message and
//   creates an MPI datatype for the payloads message
//
// side effects: allocates vector of vectors to hold items
//
// returns: total number of payload items received
//
int Neighborhoods::ExchangeNeighbors(vector<vector<char *> > &items, float wf,
				     MPI_Datatype* (*RecvItemDtype)(int *),
				     MPI_Datatype* (*SendItemDtype)(int *, 
								    char**)) {

  // total number of neighbor blocks
  nn = 0;
  for (vector<bl_t>::iterator bi = blocks.begin(); bi != blocks.end(); bi++)
    nn += bi->neighbors.size();

  PackMessages();
  PostMessages(SendItemDtype);
  TestMessages(wf, RecvItemDtype);

  // copy received items from list to output vector
  int tot_npr = ListToVector(items);

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
	pp.c.push_back(0);
	pps.push_back(pp); // invalidates iteterator pi
	pi = pps.end() - 1; // reset pi to last item
      }

      // add count information
      pi->c.push_back((ni->gb).gid);
      pi->c.push_back((ni->items).size());
      (pi->c)[0]++;

      // add payload and header to the message
      int j = 0;
      for(vector<char *>::iterator ii = ni->items.begin(); 
	  ii != ni->items.end(); ii++) {
	pi->p.push_back(*ii); // payload
	for (vector<int>::iterator hi = (ni->hdr)[j].begin(); 
	     hi != (ni->hdr)[j].end(); hi++)
	  pi->c.push_back(*hi); // header
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
// SendItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and a payloads message and
//   creates an MPI datatype for the payloads message
//
void Neighborhoods::PostMessages(MPI_Datatype* (*SendItemDtype)(int *, 
								char **)) {

  MPI_Request req;
  pl_t pl; // one payload message
  ct_t ct; // one count message
  int *rcv_ct; // one receive-counts message

  // post counts-sends, points-sends, counts-receives
  for (vector<pp_t>::iterator pi = pps.begin(); pi != pps.end(); pi++) {

    if (!pi->posted) {

      // counts-sends
      MPI_Isend(&(pi->c)[0], (pi->c)[0] * (nhdr + 2) + 1, MPI_INT, pi->proc, 
		tag * 2, comm, &req);
      pi->posted = true;
      ct.req = req;
      ct.proc = pi->proc;
      ct.tag = tag * 2;
      ct.done = false;
      ct.c = NULL; // unused for sending
      send_cts.push_back(ct);

      // payload-sends
      if (pi->p.size() > 0) { // at least one item to send
	MPI_Datatype *itype = SendItemDtype(&(pi->c)[0], &(pi->p)[0]);
	MPI_Datatype *mtype = SendMsgDtype(&(pi->c)[0], &(pi->p)[0], itype);
	MPI_Isend(MPI_BOTTOM, 1, *mtype, pi->proc, tag * 2 + 1, comm, &req);
	pl.req = req;
	pl.proc = pi->proc;
	pl.done = false;
	send_pts.push_back(pl);
	MPI_Type_free(mtype);
	delete mtype;
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
// RecvItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and
//   creates an MPI datatype for the payloads message
//
void Neighborhoods::TestMessages(float wf, 
				 MPI_Datatype* (*RecvItemDtype)(int *)) {

  int npr; // number of received points from each process
  list<ct_t>::iterator ct_it; // request list iterators
  list<pl_t>::iterator pl_it; // request list iterators
  int p; // process number
  char *rcv_p; // one payload-receive
  int i, j, k;
  MPI_Request *reqs; // pending requests
  MPI_Request *arr; // requests that arrived
  MPI_Status *stats; // statuses for arrivals
  int narr; // number of requests that arrived
  MPI_Status stat;
  int tot_narr = 0; // total number counts-receive messages arrived this round
  int nreqs; // number of requests
  if (!assign->GetStaticMode()) // override wf for dynamic repartitioning
    wf = 1.0;
  int min_arr = (int)(wf * pps.size()); // wait for this number of 
                                        // counts-receives
                                        // to arrive in this round

  if (recv_cts.size() > 0) {

    reqs = new MPI_Request[recv_cts.size()];
    arr = new MPI_Request[recv_cts.size()];
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
	    MPI_Datatype *itype = RecvItemDtype(&(ct_it->c)[0]);
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
	    pt.item_size = extent;
	    recv_pts.push_back(pt);
	    MPI_Type_free(mtype);
	    delete mtype;
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
// RecvItemDtype: pointer to user-supplied function
//   that takes a pointer to a counts message and
//   creates an MPI datatype for the payloads message
//
// returns: total number of payload items received
//
int Neighborhoods::FlushNeighbors(vector<vector<char *> > &items, 
				  MPI_Datatype* (*RecvItemDtype)(int *)) {

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
      MPI_Datatype *itype = RecvItemDtype(&(ct_it->c)[0]);
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
      pt.item_size = extent;
      recv_pts.push_back(pt);
      MPI_Type_free(mtype);
      delete mtype;
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
  int tot_npr = ListToVector(items);

  // cleanup
  for (vector<pp_t>::iterator ppi = pps.begin(); ppi != pps.end(); ppi++) {
    for (vector<char *>::iterator pi = ppi->p.begin(); pi != ppi->p.end();
	 pi++)
      delete[] *pi; // item payload, new'ed back when item was enqueued
    ppi->c.clear();
    ppi->p.clear();
  }
  pps.clear();
  for (list<ct_t>::iterator ci = recv_cts.begin(); ci != recv_cts.end(); ci++)
    delete[] ci->c; // array of received counts
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
//
// side effects: allocates vector of vectors to hold items
//
// returns: total number of items received
//
int Neighborhoods::ListToVector(vector<vector<char *> >& items) {

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
	if (ct_it->c[j * 2 + 2] > 0) {
	  gid = (ct_it->c)[2 * j + 1]; // global and local block id of dest.
	  assert((lid = Gid2Lid(gid)) != -1);
	}
	// copy items from recv_pts list to local block of output items vector
	for (int k = 0; k < (ct_it->c)[j * (nhdr + 2) + 2]; k++) {
	  char *pp = &(pl_it->p[n]);
	  items[lid].push_back(pp);
	  tot_npr++;
	  n += pl_it->item_size;
	}
      }
      // delete the messages
      pl_it = recv_pts.erase(pl_it);
      ct_it = recv_cts.erase(ct_it);
    }
    else
      pl_it++;
  }

  return tot_npr;

}
//---------------------------------------------------------------------------
//
// searches local blocks for global block id and returns local block id
// 
// gid: global block id
//
// returns: local block id, -1 if not found
//
int Neighborhoods::Gid2Lid(int gid) {

  int i;

  for (i = 0; i < (int)blocks.size(); i++) {
    if (blocks[i].gid == gid)
      break;
  }

  return(i < (int)blocks.size() ? i : -1);

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

  int size; // datatype size in bytes
  MPI_Type_size(*mtype, &size);
  if (size > 0) // allocate receive buffer
    pts = new char[size];

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

  return mtype;

}
//-----------------------------------------------------------------------
