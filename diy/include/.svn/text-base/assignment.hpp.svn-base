//---------------------------------------------------------------------------
//
// assignment class
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

#ifndef _ASSIGNMENT
#define _ASSIGNMENT

#include <vector>
#include <algorithm>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include "diy.h"

using namespace std;

//
// assignment places blocks on processes and gives blocks a unique id numbering
// i.e., the mapping of global block ids (gid) and process ids (pid) to 
//   local blocks (lid)
//
// so far, two assignment modes are supported
//
// 1. round-robin assignment (for postprocessing)
// blocks are uniquely numbered in row-major order and assigned to processes
// in round robin order based on the block gid
//
// 2. process rder assignment (for in situ)
// blocks are already mapped to processes (in situ) and gids are assigned in
// process rank order, ie, process 0's blocks are given thr first gids,
// followed by process 1's blocks, etc.
//
// neither of these modes support block reassignment
// the ZOLTAN-based reassignment code included below is from a previous project
// but has not been integrated into DIY yet. Additional assignment modes
// will support reassignment
//
// Assignment is an abstract base class that applications do not call
// applications use derived classes: RoundRobinAssignment or ProcOrderAssignment
//
class Assignment {

 public:

  ~Assignment() {}

  // local block id to global block id
  virtual int AssignGid(int lid) = 0;
  virtual int Gid2Lid(int gid) = 0;
  virtual int Gid2Proc(int gid)  = 0;

  int NumBlks() { return nb; } // number of local blocks in current domain
  int StartGid() { return start_b; } // starting gid in current domain
  int NumGids() { return tot_b; } // number of global blocks in current domain

  // static partitioning or dynamic repartitioning (static for now)
  bool GetStaticMode() { return true; }

  // MPI rank and groupsize
  int GetRank() { return rank; }
  int GetGroupsize() { return groupsize; }

protected:

  int start_b; // starting gid of this domain (num blocks in prior domains)
  int tot_b; // total number of blocks in this domain
  int nb; // my local number of blocks
  int rank; // my MPI process rank
  int groupsize; // MPI groupsize
  MPI_Comm comm; // MPI communicator

};

class RoundRobinAssignment : public Assignment {

  friend class blocking;

 public:

  RoundRobinAssignment(int start_b, int tot_b, int &nb, int &max_b, 
		       MPI_Comm comm);

  // assigns global block id
  int AssignGid(int lid)  {
    return(start_b + rank + lid * groupsize);
  }

  // global block id to local block id
  int Gid2Lid(int gid)  {
    return((gid - start_b)  / groupsize);
  }

  // global block id to process rank
  int Gid2Proc(int gid) {
    return((gid - start_b) % groupsize);
  }

};

class ProcOrderAssignment : public Assignment {

 public:

  ProcOrderAssignment(int start_b, int tot_b, int &nb, int &max_b, 
		      MPI_Comm comm);

  // assigns global block id
  int AssignGid(int lid) {
    return(start_b + tot_b / groupsize * rank + lid);
  }

  // global block id to local block id
  int Gid2Lid(int gid) {
    return(gid - start_b - tot_b / groupsize * rank);
  }

  // global block id to process rank
  int Gid2Proc(int gid) {
    int proc = (gid - start_b) / (tot_b / groupsize);
    return(proc < groupsize ? proc : groupsize - 1);
  }

};

class ExistingAssignment : public Assignment {

 public:

  ExistingAssignment(int start_b, int nb, int &max_b, int &tot_b, 
		     MPI_Comm comm);

  // local block id to global block id
  int AssignGid(int lid)  {
    lid = lid; // quiet compiler warning
    fprintf(stderr, "ExistingAssingment::Lid2Gid() should not be called. Use DIY_Gid() or Blocking::Lid2Gid() instead.\n");
    assert(0);
    return(-1);
  }

  // global block id to local block id
  int Gid2Lid(int gid)  {
    gid = gid; // quiet compiler warning
    fprintf(stderr, "ExistingAssignment::Gid2lid() is not implemented yet\n");
    assert(0);
    return(-1);
  }

  // global block id to process rank
  int Gid2Proc(int gid) {
    gid = gid; // quiet compiler warning
    fprintf(stderr, "ExistingAssignment::Gid2Proc() is not implemented yet\n");
    assert(0);
    return(-1);
  }

};

//----------------------------------------------------------------------------
//
// reassignment functions based on Zoltan library
// not part of a class, so that various functions can be used as callbacks
// todo: see if some of these can be included in the assignment class
//
//---------------------------------------------------------------------------

// todo: replace Lattice with diy

#ifdef ZOLTAN

#include <stdio.h>
#include <stdlib.h>
#include <zoltan.h>
#include "VectorMatrix.h" 
#include "Lattice4D.h"
#include "LatticeAMR.h"

// todo: this is a hack
#define MAX_NUM_BLOCKS 1000 // maximum number of block per process
#define MAX_NUM_NEIGHBORS 1000 // maximum number of neighbors a block can have
#define MAX_NUM_SEEDS 10000 // maximum number of seeds a block can have

void InitRepartition4D(void *lat, MPI_Comm comm);
void InitRepartitionAMR(void *lat, MPI_Comm comm);
void ChangePartition(int grp, int *nb, int **block_ranks, int ***neighbor_ranks,
		     int ***neighbor_procs, Partition *part, VECTOR4 ***seeds, 
		     int **size_seeds, int **num_seeds, int *avg_neigh,
		     int *alloc_blocks, int **alloc_neighbors, MPI_Comm comm,
		     OSUFlow ***osuflow,
		     void (*add_neighbor)(int, int, int, int *, int *, int ***,
					  int ***, Partition *),
		     int *wgts = NULL);
void RemoveBlock(int myrank, int newproc, int *nb, int **block_ranks,
		 int **alloc_neighbors, int ***neighbor_ranks,
		 int ***neighbor_procs, Partition *part, 
		 VECTOR4 ***seeds, int **size_seeds, int **num_seeds);
void AddBlock(int myrank, int num_neighbors, int *neighbors, int *alloc_blocks,
	      int *nb, int **block_ranks, int **alloc_neighbors, 
	      int ***neighbor_ranks, int ***neighbor_procs, Partition *part,
	      int myproc, VECTOR4 ***seeds, int **size_seeds, int **num_seeds,
	      OSUFlow ***osuflow,
	      void (*add_neighbor)(int, int, int, int *, int *, 
				   int ***, int ***, Partition *));
void ExchangeExports(int nexport, ZOLTAN_ID_PTR export_gids,
		     ZOLTAN_ID_PTR export_lids, int * export_procs,int nb, 
		     int *block_ranks, int **neighbor_ranks, 
		     int **neighbor_procs, Partition *part, MPI_Comm comm);

// zoltan callbacks to be used with Lattice4D
int GetNumberofAssignedObjects4D(void *lat, int *err);
void GetObjectList4D(void *lat, int ngids, int nlids, 
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
		   int wgt_dim, float *obj_wgts, int *err);
int GetObjectSize4D(void *lat, int *err);
void GetObjects4D(void *lat, int ngids, int nlids, int nobjs, 
		ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int ndim,
		double *pts, int *err);


// zoltan callbacks to be used with LatticeAMR
int GetNumberofAssignedObjectsAMR(void *lat, int *err);
void GetObjectListAMR(void *lat, int ngids, int nlids, 
		   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
		   int wgt_dim, float *obj_wgts, int *err);
int GetObjectSizeAMR(void *lat, int *err);
void GetObjectsAMR(void *lat, int ngids, int nlids, int nobjs, 
		ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int ndim,
		double *pts, int *err);

// utility functions
int IsBlockInTimeGroup4D(void *lat, int grp, int blk);
int IsBlockInTimeGroupAMR(void *lat, int grp, int blk);

#endif

#endif
