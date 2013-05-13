//
// 4D lattice header
//
// Han-Wei Shen
// The Ohio State University
// Columbus, OH
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

#ifndef _LATTICE4D_H_
#define _LATTICE4D_H_ 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <list>
#include <iterator>
#include <math.h>
#include "Partition.h"
#include "Repartition.h"

#ifdef _MPI
#include <mpi.h>
#endif

struct time_bounds_t {
  int tmin;
  int tmax;
};

struct volume_bounds_t {
  int xmin, xmax; 
  int ymin, ymax;
  int zmin, zmax;
  int tmin, tmax;
  int xdim, ydim, zdim, tdim;
};

class  Lattice4D {

 public: 
  
  Lattice4D(int xlen, int ylen, int zlen, int tlen, int nsp, 
	    int *ntp, int nid = 1, int myid = 0, 
	    bool track_ids = false, int ghost  = 0); 
  Lattice4D(char *part_file, int xlen, int ylen, int zlen, int tlen,
	    int nsp, int *ntp, int nid, int myid,
	    bool track_ids = false, int ghost = 0);

  ~Lattice4D(); 
  int GetRank(int i, int j, int k, int l);
  int GetRank(float x, float y, float z, float t);
  int MyGetRank(float x, float y, float z, float t);
  int GetRank(int i);
  int GetRealRank(float x, float y, float z, float t);
  int GetIndices(int rank, int &i, int &j, int&k, int&t);
  int GetIndices(int rank, int &i, int &j, int&k);
  int GetIndices(float x, float y, float z, float t, int& iidx, int &jidx, 
		 int&kidx, int& lidx);
  int GetBounds(int i, int j, int k, int t, volume_bounds_t& vb);
  int GetBounds(int rank, volume_bounds_t& vb);
  int GetRealBounds(int i, int j, int k, int t, volume_bounds_t& vb);
  int GetRealBounds(int rank, volume_bounds_t& vb);
  volume_bounds_t* GetBoundsList(){ return vb_list;}
  void GetLatticeDims(int& i, int&j, int &k, int &l) 
  {i = idim; j=jdim; k=kdim; l = tdim; }
  bool isIn(float, float, float, float, int, int, int, int); 
  bool isIn(float, float, float, int, int, int); 
  bool isInReal(float, float, float, float, int, int, int, int); 
  bool isInReal(float, float, float, int, int, int); 

  int GetNeighbor(int myrank, float x, float y, float z, float t, 
		  int &ei, int &ej, int &ek, int &et); 

  int GetNeighbor(int myrank, float x, float y, float z, float t);
  int CheckNeighbor(int myrank, float x, float y, float z, float t); 
  int CheckNeighbor(int myrank, float x, float y, float z); 
  int myGetNeighbor(int myrank, float x, float y, float z, float t, int *i, 
		    int *j, int *k, int *l);
  int GetProc(int, int, int, int); 
  void RoundRobin_proc(); 
  void Explicit_proc(int *block_procs);
  void GetPartitions(int, int**, int&); 
  void GetPartitions(int proc, int*p_list);
  void NeighborIndices(int n, int i, int j, int k, int l, int &in,
		       int &jn, int &kn, int &ln);
  int GetMyNumPartitions() { return nb; }
  int GetMyNumNeighbors() { return avg_neigh; }
  double GetMyCommTime() { return comm_time; }
  int GetMyTotPtsSend() { return tot_pts_send; }
  void GetVB(int block, float *min_s, float *max_s, int *min_t, int *max_t);
  void GetRealVB(int block, float *min_s, float *max_s, int *min_t, int *max_t);
  void GetTB(int block, int *min_t, int *max_t);
  void GetRealTB(int block, int *min_t, int *max_t);
  void GetTimeGroupBounds(int group, int* min_t, int* max_t);
  void GetRealTimeGroupBounds(int group, int* min_t, int* max_t);
  void GetGlobalVB(int part, float *min_s, float *max_s, 
		   int *min_t, int *max_t);
  void GetExtents(float *min, float *max);
  void swap4(char *n);

  // wrappers around partition methods
/*   void SetReq(int block); */
/*   void ClearReq(int block); */
/*   int GetReq(int block); */
  void SetLoad(int block);
  void ClearLoad(int block);
  int GetLoad(int block);
/*   void SetComp(int block, int iter_num); */
/*   void ClearComp(int block); */
/*   int GetComp(int block, int iter_num); */
  void PostPoint(int block, float *p, int recirc, int64_t seed_id = -1);
/*   void PrintPost(int block); */
/*   int OldExchangeNeighbors(float ***seeds, int *alloc_seeds, int *num_seeds,  */
/* 			int64_t **seed_ids = NULL); */
/*   int SyncExchangeNeighbors(int *nproc, float ***pts, int ***cts, */
/* 			    int64_t ***pids = NULL); */

  int xdim, ydim, zdim, ldim; // total grid size
  int idim, jdim, kdim, tdim; //the lattice range (number of partitions)
  int npart; // total number of partitions in the domain
  int nb; // number of my blocks
  int **neighbor_ranks; // ranks of neighbors for my blocks
  int **neighbor_procs; // procs of neighbors for my blocks
  class Partition *part;

  // benchmarking
  int avg_neigh; // average number of neighbors per block in my process
  double comm_time; // communication time for my process
  int tot_pts_send; // total points sent from my process

  int *alloc_neighbors; // allocated size of neighbor_ranks, neighbor_procs
  int alloc_blocks; // allocated number of blocks
  int *block_ranks; // rank (global partition number) of each of my blocks

 private: 

  int myproc; // my process or thread number
  int nproc; // number of processes or threads
  volume_bounds_t* vbr_list; // real volume bounds w/o ghost
  volume_bounds_t* vb_list;  // volume bounds w/ ghost
  time_bounds_t * tbr_list;  // time bounds for each time block w/o ghost cells
  time_bounds_t * tb_list;   // time bounds for each time block w/ ghost cells
  int* flowMatrix; 
  void GetNeighborRanks(int block);
  void VolumeBounds(float *block_extents, int nblocks, int *block_size,
		    volume_bounds_t *vb_list, int &idim, int &jdim, 
		    int &kdim);
  void ComputePartition(int *data_dim, int ghost, 
				    int nsp, int ntp, int *lat_dim);
  void ApplyGhost(int ghost);
  int GetNumPartitions(int proc);
  void FindTimeBounds();

  // function to sort time bounds
  static int compare_time_bounds(const void* a, const void* b); 

  // overall extents of the entire dataset
  float min_extent[4];
  float max_extent[4];

}; 

// not a member of Lattice4D so that it can be used as a callback function
void AddNeighbor(int myblock, int neighrank, int neighproc, int *block_ranks,
		 int *alloc_neighbors, int ***neighbor_ranks,
		 int ***neighbor_procs, Partition *part);

#endif 
 
