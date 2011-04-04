//------------------------------------------------------------------------------
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

#ifndef _NEIGHBORHOOD
#define _NEIGHBORHOOD 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <list>
#include <vector>
#include <iterator>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "blocking.hpp"
#include "assignment.hpp"

using namespace std;

// todo: review these limits; hopefully remove them
#define MAX_BLOCK 64
#define MAX_NEIGHBORS 64
#define MAX_MSG_SIZE 256

// one neighboring block
struct ne_t {
  gb_t gb; // global id and process id of the neighbor block
  vector<char *> items; // items to send to the neighbor
  vector<vector <int> > hdr; // optional header info, eg. array of
                             // counts of subitems in each item
};

// my block in the center of a neighborhood
struct bl_t {
  int gid; // global id of this block
  vector<ne_t> neighbors; // neighboring blocks
};

// counts message
struct ct_t {
  MPI_Request req; // MPI request for this message
  int proc; // sending process of this message
  int tag; // tag in for this message
  bool done; // message completed or not
  int *c; // received counts
};

// payload message
struct pl_t {
  MPI_Request req; // MPI request for this message
  int proc; // sending process of this message
  int tag; // tag in for this message
  bool done; // message completed or not
  char *p; // received items (char * for any generic byte)
  int item_size; // size of one item (in bytes, all items same size)
};

// package ready to ship to a single process
struct pp_t {
  int proc; // process id
  vector<int> c; // counts message
  vector<char *> p; // payload message
  bool posted; // a send has been posted for this package
};

class Neighborhoods {

 public: 
  
  Neighborhoods(Blocking *blocking, Assignment *assignment,MPI_Comm comm, 
		int nhdr = 0); 
  ~Neighborhoods(); 
  void EnqueueItem(int lid, char *item, size_t size, int neigh_gid, 
		   int *hdr = NULL);
  int ExchangeNeighbors(vector<vector<char *> > &items, float wf,
			MPI_Datatype* (*RecvItemDtype)(int *),
			MPI_Datatype* (*SendItemDtype)(int *, char**));
  int FlushNeighbors(vector<vector<char *> > &items,
		     MPI_Datatype* (*RecvItemDtype)(int *));
  int TotItemsSent() { return tot_items_sent; }
  int Gid2Lid(int gid);
  int Lid2Gid(int lid) { return blocks[lid].gid; }

 private: 

  void PostMessages(MPI_Datatype* (*SendItemDtype)(int *, char**));
  void TestMessages(float wf, MPI_Datatype* (*RecvItemDtype)(int *));
  void PackMessages();
  int ListToVector(vector<vector<char *> >& items);
  MPI_Datatype* RecvMsgDtype(int *cts, char* &pts, MPI_Datatype *itype);
  MPI_Datatype* SendMsgDtype(int *cts, char** pts, MPI_Datatype *itype);

  vector<bl_t> blocks; // my blocks
  vector<pp_t> pps; // packed messages organized by process
  int nn; // total number of neighbor blocks for this process
  int tag; // tag number for matching count- and point-receives
  int rank; // my process rank in MPI communicator
  int groupsize; // size of MPI communicator
  MPI_Comm comm; // MPI communicator
  Blocking *blocking; // blocking object
  Assignment *assign; // assignment object
  int nhdr; // number of optional header counts
  int tot_items_sent; // total number of items sent

  // message lists
  list<ct_t> send_cts; // MPI requests for sending counts
  list<pl_t> send_pts; // MPI requests for sending payloads
  list<ct_t> recv_cts; // MPI requests for receiving counts
  list<pl_t> recv_pts; // MPI requests for reciving payloads

}; 

#endif 
