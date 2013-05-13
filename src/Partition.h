//------------------------------------------------------------------------------
//
// Partition header
// manages parts
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// All rights reserved. May not be used, modified, or copied
// without permission
//
//--------------------------------------------------------------------------

#ifndef _PARTITION_H_
#define _PARTITION_H_ 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <list>
#include <vector>
#include <iterator>
#include <math.h>
#include <string.h>

#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

// ADD-BY-LEETEN 07/01/2011-BEGIN
// MOD-BY-LEETEN 2011/07/01-FROM:
	// #if !defined( WIN32 )
// TO:
#if defined( WIN32 )
// MOD-BY-LEETEN 2011/07/01-END
typedef long long int64_t;
#endif
// ADD-BY-LEETEN 07/01/2011-END

/* #define MAX_BLOCK 64 */
/* #define MAX_NEIGHBORS 64 */
/* #define MAX_MSG_SIZE 256 */
/* #define WAIT_FACTOR 0.75 // wait for this fraction of pending messages to arrive */

// a global (all processes) block or partition
struct Partition4D {

  // the following arrays are indexed according to neighbor number
  int *NumSendPoints; // number of points ready to send
  int *SizeSendPoints; // size of sending points list (bytes)
  float **SendPoints; // sending points list
  int64_t **SendPointIds; // sending point ids
  int NumNeighbors; // the number of neighbors used
  int AllocNeighbors; // number of neighbors allocated
  int Proc; // process(or) number (mpi rank, core number, node number, etc.)
  int Loaded; // data are loaded

};

class Partition {

 public: 
  
  Partition(int npart, int nproc, int myproc, bool track_ids = false); 
  ~Partition(); 

/*   void SetReq(int gid) { parts[gid].Requed = 1; } */
/*   void ClearReq(int gid) { parts[gid].Requed = 0; } */
/*   int GetReq(int gid) { return parts[gid].Requed;} */
  void SetLoad(int gid) { parts[gid].Loaded = 1; }
  void ClearLoad(int gid) { parts[gid].Loaded = 0; }
  int GetLoad(int gid) { return parts[gid].Loaded;}
/*   void SetComp(int gid, int iter_num) { parts[gid].Comped = iter_num; } */
/*   void ClearComp(int gid) { parts[gid].Comped = -1; } */
/*   int GetComp(int gid, int iter_num) { return(parts[gid].Comped >=  */
/* 						    iter_num); } */
  void PostPoint(int gid, float *p, int neighbor, int64_t seed_id = -1);
/*   void PrintPost(int gid); */
/*   void PrintRecv(int gid); */
/*   void GetRecvPts(int gid, float **ls); */
/*   int GetNumNeighbors(int gid) { return parts[gid].NumNeighbors; } */
  void AddNeighbor(int gid);
/*   void RemoveBlock(int gid); */
/*   void AddBlock(int gid); */

/* #ifdef _MPI */

/*   int OldExchangeNeighbors(int *block_gids, int **neighbor_gids, */
/* 			   int **neighbor_procs, float ***seeds,  */
/* 			   int *alloc_seeds, */
/* 			   int *num_seeds, MPI_Comm comm = MPI_COMM_WORLD, */
/* 			   int64_t **seed_ids = NULL); */
/*   int SyncExchangeNeighbors(int *block_gids, int **neighbor_gids,  */
/* 			    int **neighbor_procs, int *nproc,  */
/* 			    float ***pts, int ***cts, int64_t ***pids = NULL, */
/* 			    MPI_Comm comm = MPI_COMM_WORLD); */

/* #endif */

  Partition4D *parts; // global partition list
  int nb; // number of local blocks

 private: 

/* #ifdef _MPI */

/*   void SyncExchangeNeighborCounts(int *block_gids, int **neighbor_gids,  */
/* 				  int **neighbor_procs, int nn,  */
/* 				  int *nps, int *npr, int ***cts, int *proc,  */
/* 				  int *nproc, MPI_Comm comm = MPI_COMM_WORLD); */
/*   void SyncExchangeNeighborPoints(int *block_gids, int **neighbor_gids, */
/* 				  int **neighbor_procs, int *proc,  */
/* 				  int *nproc, int *nps, int *npr,  */
/* 				  float ***recv_pts,  */
/* 				  int64_t ***recv_pids = NULL, */
/* 				  MPI_Comm comm = MPI_COMM_WORLD); */
/*   void AsyncPostMessages(int *block_gids,  */
/* 			 int **neighbor_gids,  */
/* 			 int **neighbor_procs,  */
/* 			 int *nps, MPI_Comm comm); */
/*   void AsyncTestMessages(int *block_gids,  */
/* 			 int **neighbor_procs, int *nps,  */
/* 			 MPI_Comm comm); */

/* #endif */

/*   void PackCounts(int *block_gids, int **neighbor_gids,  */
/* 		  int **neighbor_procs, int *nps); */
/*   void PackPoints(int *block_gids, int **neighbor_procs, int *nps); */
/*   int ListToArray(int *nproc, float ***pts, int ***cts, int64_t ***pids); */

  int npart; // total number of partitions
  int nproc; // total number of processes
  int myproc; // my process id
  bool track_ids; // tracking point ids
  int nnproc; // number of processes in my neighborhood
  vector<int> proc; // the process ids of my neighboring processes
  int nn; // total number of neighbor blocks for this process
  int **snd_ct; // send-counts messages
  float **snd_pt; // sending points
  int64_t **snd_pd; // sending point ids
  int *mlen; // length of sending messages
  int tag; // tag number for matching count- and point-receives

}; 

#endif 
