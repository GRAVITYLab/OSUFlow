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
#include <stddef.h>
#include <assert.h>
#include <list>
#include <vector>
#include <iterator>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include "blocking.hpp"
#include "assignment.hpp"
#include "diy.h"

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
  unsigned char wrap_dir; // wraparound neighbor, and wrapping direction
};

// my block in the center of a neighborhood
struct bl_t {
  int gid; // global id of this block
  vector<ri_t> rem_data; // remote data used for discovering neighoring blocks
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
  int tag; // tag for this message
  bool done; // message completed or not
  char *p; // received items (char * for any generic byte)
  int item_size; // size of one item (in bytes, all items same size)
};

// payload package ready to ship to a single process
struct pp_t {
  int proc; // process id
  int* c; // counts message
  int cn; // size of counts message array
  vector<char *> p; // payload message
  bool posted; // a send has been posted for this package
};

struct ep_t {
  int proc; // process id
  vector<float> msg; // gids and bounds of blocks that neighbor this process
};

class Neighborhoods {

 public: 
  
  Neighborhoods(int did, Blocking *blocking, Assignment *assignment,
		MPI_Comm comm, bool wrap, int nhdr = 0); 
  Neighborhoods(int did, Blocking *blocking, Assignment *assignment, 
		struct ri_t **rem_ids, int *num_rem_ids, int **vids,
		int *num_vids, gb_t **neighbors, int *num_neighbors,
		MPI_Comm comm, bool wrap, int nhdr = 0);
  ~Neighborhoods(); 
  void WrapNeighbors();
  void EnqueueItem(int lid, char *item, size_t size, int neigh_gid, 
		   int *hdr, void (*TransformItem)(char *, unsigned char));
  void EnqueueItemAll(int lid, char *item, size_t size, int *hdr,
		      void (*TransformItem)(char *, unsigned char),
		      bool symmetrical, bool self);
  void EnqueueItemAllNear(int lid, char *item, size_t size,
			  float *pt, float dist, int *hdr, 
			  void (*TransformItem)(char *, unsigned char),
			  bool symmetrical);
  void EnqueueItemDir(int lid, char *item, size_t size, int *hdr, 
		      void (*TransformItem)
		      (char *, unsigned char),
		      unsigned char neigh_dir);
  int ExchangeNeighbors(vector<vector<char *> > &items, float wf,
			void (*ItemDtype)(MPI_Datatype *),
			bool discovery = false);
  int FlushNeighbors(vector<vector<char *> > &items,
		     void (*ItemDtype)(MPI_Datatype *),
		     bool discovery = false);
  int Pt2NeighGid(int lid, float *pt);
  void BoundsIntersectNeighbors(int lid, bb_t cell_bounds, float t, 
				int *num_intersect, int *gids_intersect);

 private: 

  void PostMessages(void (*ItemDtype)(MPI_Datatype*));
  void TestMessages(float wf, void (*ItemDtype)(MPI_Datatype *));
  void PackMessages();
  int ListToVector(vector<vector<char *> >& items, bool discovery);
  MPI_Datatype* RecvMsgDtype(int *cts, char* &pts, MPI_Datatype *itype);
  MPI_Datatype* SendMsgDtype(int *cts, char** pts, MPI_Datatype *itype);
  void GetNeighborBounds();
  void SetNeighborBounds(int lid, gb_t *gb);
  void GetNeighbors(int **vids, int *num_vids);
  int DiscoverGid2Lid(gb_t *gb, vector<int> &lids);

  int did; // domain id
  vector<bl_t> blocks; // my blocks
  vector<pp_t> pps; // packed messages organized by process
  vector<ep_t> eps; // packed extents organized by process
  int nn; // total number of neighbor blocks for this process
  int tag; // tag number for matching count- and point-receives
  int rank; // my process rank in MPI communicator
  int groupsize; // size of MPI communicator
  MPI_Comm comm; // MPI communicator
  Blocking *blocking; // blocking object
  Assignment *assign; // assignment object
  int nhdr; // number of optional header counts
  ri_t **rem_ids; // remote ids for each local block (for neighbor discovery)
  int *num_rem_ids; // number of remote ids for each local block
  bool discovery; // whether neighbor discovery is being used

  // message lists
  list<ct_t> send_cts; // MPI requests for sending counts
  list<pl_t> send_pts; // MPI requests for sending payloads
  list<ct_t> recv_cts; // MPI requests for receiving counts
  list<pl_t> recv_pts; // MPI requests for reciving payloads

}; 

// callback function not part of the class
static void Nbhds_ItemType(MPI_Datatype *type);

#endif 
