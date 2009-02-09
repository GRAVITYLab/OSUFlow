
#ifndef _LATTICEAMR_H_
#define _LATTICEAMR_H_ 
    
#include<stdio.h>
#include<stdlib.h>

#include "header.h"
#include "calc_subvolume.h"
#include "VectorMatrix.h" 
#include "OSUFlow.h"

#ifdef MPI
#include <mpi.h>
#endif

struct PartitionAMR4D {
  int Proc; // process(or) number (mpi rank, core number, node number, etc.)
};

class  LatticeAMR {

 public: 
  
  //x/y/zlen are the physical dimensions in the domain, tlen is the total number of 
  //time steps 
  LatticeAMR(float xlen, float ylen, float zlen, float tlen, int total_level); 

  ~LatticeAMR(); 

  //create a new AMR level. x/y/zsize are the physical dimensions of the blocks
  bool CreateLevel(int level, float xsize, float ysize, float zsize, float tsize,
		   float x_min, float x_max, float y_min, float y_max, 
		   float z_min, float z_max, float t_min, float t_max); 

  //Check in a non-empt block in the level 
  bool CheckIn(int level, float x, float y, float z, float t, 
	       float* data); 

  //Go through all levels to collect non-empty blocks 
  void CompleteLevels(); //extract those non-empty blocks 


  int GetIndexinLevel(int level, float x, float y, float z, float t);

  int GetCoordsinLevel(int level, float x, float y, float z, float t, 
			int& i, int &j, int &k, int &t); 

  int GetFinestLevel(float x, float y, float z, float t);

  int GetRank(float x, float y, float z, float t); 

  int GetRank(int i, int j, int k, int t, int level);         // serach by lattice indices

  int GetIndices(int rank, int &i, int &j, int&k, int&t, int&level);  // search by rank 

  int GetIndices(float x, float y, float z, float t, int& iidx, int &jidx, int&kidx, int& tidx, int& level);  // search by physical locations


  bool isIn(float x, float y, float , float t, int i, int j, int k, int l, 
	    int level); 

  int GetBounds(int i, int j, int k, int t, 
		int level, volume_bounds_type& vb);  //return bound 

  int GetBounds(int rank, volume_bounds_type& vb);   // return bound 

  volume_bounds_type* GetBoundsList(int& num){ num = npart; return vb_list;}

  float *GetDataPtr(int rank); 

  void GetLatticeDims(int& i, int&j, int &k, int &l, int level) 
  {i = idim[level]; j=jdim[level]; k=kdim[level]; l = ldim[level]; }

  int CheckNeighbor(int myrank, float x, float y, float z, float t); 

  int GetNeighbor(int myrank, float x, float y, float z, float t, int &i, int &j,
		  int &k, int&l, int&level);

  void InitSeedLists(); 
  void ResetSeedLists(); 
  void ResetSeedLists(int i); 
  bool InsertSeed(int, int, int, int, int, VECTOR4); 
  bool InsertSeed(int rank, VECTOR4); 

  int GetProc(int, int, int, int, int); 
  int GetProc(int); 
  void GetPartitions(int, int**, int&); 
  void GetPartitions(int, int*, int&); 
  void RoundRobin_proc(int n); 
  
  list<VECTOR4> *seedlists; 

 private: 

  int num_levels;  // the AMR refinement level 

  float xdim, ydim, zdim, tdim; //the physical dimensions of the whole domain 

  //the physical ranges, in all levels
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax, *tmin, *tmax; 

  float *xlength, *ylength, *zlength, *tlength; // the physical size of the 
                                        // block in each level

  int *idim, *jdim, *kdim, *ldim; //the lattice dimensions, in all levels 

  bool **has_data;      // whether the elements in this level have data or not 
  int  **finest_level;  // in which level the finest data are defined 
  int  **index_to_rank; 
  float ***data_ptr; // data pointers for every non-empty block, ordered by 
                     // level and index (so three *s here) 


  int* nblocks; // number of blocks in each level regardless of having data or not 

  int npart; 
  volume_bounds_type *vb_list; 
  int *rank_to_index; // from rank to (i,j,k,t, level) indices 

  PartitionAMR4D *parts; 
}; 

#endif 
