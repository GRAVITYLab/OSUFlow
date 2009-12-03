
#ifndef _LATTICEAMR_H_
#define _LATTICEAMR_H_ 
    
#include<stdio.h>
#include<stdlib.h>

#include "header.h"
#include "calc_subvolume.h"
#include "VectorMatrix.h"
#include "OSUFlow.h"
#include "Partition.h"
#include "FlashAMR.h"

#ifdef _MPI
#include <mpi.h>
#endif

class  LatticeAMR {

 public: 
  
  //x/y/zlen are the physical dimensions in the domain, tlen is the total number of 
  //time steps 
  LatticeAMR(char *filename, int tlen, char * vx, char * vy, char *vz,
	     int nid = 1, int myid = 0); 

  ~LatticeAMR(); 

  //create a new AMR level. x/y/zsize are the physical dimensions of the blocks 
  //in this level 
  //we assume a block will contain one time step worth of data 
  bool CreateLevel(int level, float xsize, float ysize, float zsize, 
		   int x_res, int y_res, int z_res, 
		   float x_min, float x_max, float y_min, float y_max, 
		   float z_min, float z_max, int t_min, int t_max); 

  // Check in a non-empty block in the level, along with its data
  bool CheckIn(int level, float x, float y, float z, int t, 
	       float* data); 

  //Go through all levels to collect non-empty blocks 
  void CompleteLevels(); //extract those non-empty blocks 
  void CompleteLevels(int t_interval); 
  void MergeAndCompleteLevels(); 

  // assign contiguous division of blocks to processes and get the process
  // id back based on that division
  int AssignContig(int start_time, int end_time,
		   int *start_block, int *end_block);
  int GetProc(int level, int idx);

  // time-varying AMR object
  TimeVaryingFlashAMR *tamr;

  int GetRank(float x, float y, float z, float t); 
  int GetRank(float x, float y, float z, float t, int level);
  int GetRank(int i, int j, int k, int t, int level);   
  bool isIn(float x, float y, float , float t, int i, int j, int k, int l, 
	    int level); 
  int GetBounds(int i, int j, int k, int t, 
		int level, volume_bounds_type_f& vb);  //return bound 
  int GetBounds(int rank, volume_bounds_type_f& vb);   // return bound 
  volume_bounds_type_f* GetBoundsList(int& num){ num = npart; return vb_list;}
  float**GetDataPtr(int rank); 
  void GetLatticeDims(int& i, int&j, int &k, int &l, int level) 
  {i = idim[level]; j=jdim[level]; k=kdim[level]; l = ldim[level]; }
  int CheckNeighbor(int myrank, float x, float y, float z, float t); 
  int GetNeighbor(int myrank, float x, float y, float z, float t, int &i, 
		  int &j, int &k, int&l, int&level);
  int GetFinestLevel(float x, float y, float z, float t);
  void InitSeedLists(); 
  void ResetSeedLists(); 
  void ResetSeedLists(int i); 
  bool InsertSeed(int, int, int, int, int, VECTOR4); 
  bool InsertSeed(int rank, VECTOR4); 
  int GetProc(int, int, int, int, int); 
  int GetProc(int); 
  void GetPartitions(int, int**, int&); 
  bool Mergeable(int i, int j, int k, int t, int level, int& mergeLevel); 
  void MergeBlocks(); 
  void RoundRobin_proc(); 
  void ContigRange_proc();
  void GetBlockDims(int *dims);
  
  list<VECTOR4> *seedlists; 

  int GetMyNumPartitions(int proc);
  int GetMyNumPartitions();
  void GetMyPartitions(int proc, int* p_list);
  int GetTotalNumPartitions() { return npart; }
  int GetMyNumNeighbors() { return avg_neigh; }
  int GetMyIOBW() { return io_bw; }
  double GetMyIOTime() { return io_time; }
  double GetMyCommTime() { return comm_time; }
  void GetVB(int block, float *min_s, float *max_s, 
	     int *min_t, int *max_t);
  void GetGlobalVB(int part, float *min_s, float *max_s, 
		   int *min_t, int *max_t);
  int GetNeighbor(int myrank, float x, float y, float z, float t);
  void GetNeighborRanks(int block);

  // wrappers around partition methods
  float** GetData(int block);
  void SetReq(int myrank);
  void ClearReq(int myrank);
  int GetReq(int myrank);
  void SetLoad(int myrank);
  void ClearLoad(int myrank);
  int GetLoad(int myrank);
  void SetComp(int myrank, int iter_num);
  void ClearComp(int myrank);
  int GetComp(int myrank, int iter_num);
  void PostPoint(int myrank, VECTOR4 p);
  void PrintPost(int myrank);
  void PrintRecv(int myrank);
  int ExchangeNeighbors(VECTOR4 **seeds, int *num_seeds);

 private: 

  int GetCoords(int rank, int &i, int &j, int&k, int&t, int&level);  // search by rank 

  int GetCoords(float x, float y, float z, float t, int& iidx, int &jidx, int&kidx, int& tidx, int& level);  // search by physical locations

  int GetIndexinLevel(int level, float x, float y, float z, float t);

  int GetCoordsinLevel(int level, float x, float y, float z, float t, 
			int& i, int &j, int &k, int &l); 

  bool MapCells(int fromI, int fromJ, int fromK, int fromT, 
		 int fromlevel, int toLevel, int& toImin, int& toImax, 
		 int& toJmin, int& toJmax, int& toKmin, int& toKmax, 
		 int& toTmin, int& toTmax); 

  int block_dims[3]; // number of voxels per block eg 16x16x16
  int num_levels;  // the number of AMR levels
  float xdim, ydim, zdim; //the lengths of the whole domain (not the data resolution) 
  int tdim;               //number of time steps, always an integer

  // the physical bounds of all levels
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax; 
  int   *tmin, *tmax; 

  float *xlength, *ylength, *zlength; // the physical lengths of blocks 
  int   *tlength;                     // in all levels (all blocks have the 
                                      // same size in each level)

  int *xres, *yres, *zres, *tres; // resolution of data blocks in all levels 

  int *idim, *jdim, *kdim, *ldim; // the lattice resolutions, in all levels 

  bool **has_data;      // whether the block in this level has data or not 
  int ** has_data_from_merger; 
  int  **finest_level;  // in which level the finest data are defined 
  int  **index_to_rank; // mapping index in a level to a partition rank
  int  **index_to_seq; // mapping index in a level to a block sequence number
  float ***data_ptr; // data pointers for every non-empty block, ordered by 
                     // level and index (so three *s here) 

  int* nblocks; // number of blocks in each level (regardless of having data or not) 

  int tot_nblocks; // total number of blocks in a time interval
  int avg_neigh; // average number of neighbors per block in my process
  int npart;   // number of partitions. 
  int io_bw; // I/O bandwidth in MB/s
  double io_time, comm_time; // I/O and communication time for my process
  volume_bounds_type_f *vb_list; // bounds for each partition 
  int *rank_to_index; // from rank to (i,j,k,t, level) indices 

  int *block_ranks; // rank (global partition number) of each of my blocks
  int myproc; // my process or thread number
  int nproc; // number of processes or threads
  class Partition *part; // partition class object
  int **neighbor_ranks; // ranks of neighbors for my blocks
  int nb; // number of my blocks
  void GetCoarseNeighborRanks(int block);
  void GetFineNeighborRanks(int block);
  void AddNeighbor(int myblock, int neighrank);

}; 

#endif 
