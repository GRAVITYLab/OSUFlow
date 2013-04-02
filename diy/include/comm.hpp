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

#ifndef _COMM
#define _COMM

#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include "diy.h"
#include "mpi.h"
#include "assignment.hpp"

using namespace std;

struct snd_item { // send-item message (includes both header and item)
  int *hdr; // header
  MPI_Request hreq; // header send request
  MPI_Request ireq; // item send request
};

struct rcv_hdr { // receive-header message
  int *hdr; // header
  MPI_Request req; // header receive request
  int tag; // message tag
  int proc; // source process rank
};

class Comm {

 public:

  Comm(MPI_Comm comm);
  ~Comm();
  void SendItem(char *item, int *hdr, int dest_rank, int gid,
		MPI_Datatype *dtype);
  void StartRecvItem(int src_rank, bool use_header);
  int FinishRecvItemsMerge(char **items, int *gids, int *procs,
			   char * (*CreateItem)(int *),
			   void (*RecvItemDtype)(void*, MPI_Datatype*, int *));
  int RecvItemsMerge(char **items, int *gids, int *procs, float wf,
		     char * (*CreateItem)(int *),
		     void (*RecvItemDtype)(void*, MPI_Datatype*, int *));
  int FinishRecvItemsSwap(char **items, int *gids, int *procs, int ne,
			  char * (*CreateItem)(int *, int),
			  void (*RecvItemDtype)(void*, MPI_Datatype*, 
						 int));
  int RecvItemsSwap(char **items, int *gids, int *procs, float wf, int ne,
		    char * (*CreateItem)(int *, int),
		    void (*RecvItemDtype)(void*, MPI_Datatype*, 
					   int));
  void Send(void *item, int count, DIY_Datatype datatype, 
	    int dest_gid, vector <Assignment*> assign);
  int Recv(int my_gid, void** &items, DIY_Datatype datatype, int wait,
	   int *sizes);
  void FlushSendRecv(int barrier);

#ifdef _MPI3

  void RmaSend(void *item, int count, DIY_Datatype datatype, 
	       int my_gid, int dest_gid, vector <Assignment*> assign);
  int RmaRecv(int my_gid, void** &items, DIY_Datatype datatype, 
	      int *src_gids, int wait, Assignment *assign, int *sizes);
  void RmaFlushSendRecv(int barrier);

#endif

 private:

  int FinishRecvItems(char **items, int *gids, int *procs, int ne,
		      char* (*CreateItemMerge)(int *),
		      char* (*CreateItemSwap)(int *, int),
		      void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *),
		      void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int));
  int RecvItems(char **items, int *gids, int *procs, float wf, int ne,
		char* (*CreateItemMerge)(int *),
		char* (*CreateItemSwap)(int *, int),
		void (*RecvItemDtypeMerge)(void*, MPI_Datatype*, int *),
		void (*RecvItemDtypeSwap)(void*, MPI_Datatype*, int));

  bool use_header; // using a header for item communication
  MPI_Comm comm; // communicator
  int rank, groupsize; // MPI usual
  vector<snd_item> snd_items; // posted send-item messages
  vector<rcv_hdr> rcv_hdrs; // posted receive-header messages
  MPI_Datatype *send_dtypes; // datatypes for sending items
  int num_sends; // number of sends in the current round
  int *rma_buf; // RMA buffer
  MPI_Win rma_win; // RMA window
  vector<MPI_Request> rma_reqs; // RMA nonblocking send requests

};

#endif
