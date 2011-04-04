//------------------------------------------------------------------------------
//
// Repartitioning header
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

#ifndef _REPARTITION_H_
#define _REPARTITION_H_ 
    
#ifdef ZOLTAN

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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
 
