
#ifndef _LATTICE4D_H_
#define _LATTICE4D_H_ 
    
#include<stdio.h>
#include<stdlib.h>

#include "header.h"
#include "calc_subvolume.h"
#include "VectorMatrix.h" 
#include "OSUFlow.h"

#ifdef MPI
#include <mpi.h>
#endif

struct Partition4D {
  int NumSendPoints[6]; // number of points ready to send
  int SizeSendPoints[6]; // size of sending points list (bytes)
  float *SendPoints[6]; // sending points list
  int NumRecvPoints[6]; // number of points received
  int SizeRecvPoints[6]; // size of receiving points list (bytes)
  float *RecvPoints[6]; // receiving points list
  int Proc; // process(or) number (mpi rank, core number, node number, etc.)
};

class  Lattice4D {

 public: 
  
  Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp); 
  ~Lattice4D(); 
  int GetRank(int i, int j, int k, int l);         // serach by lattice indices
  int GetRank(float x, float y, float z, float t);  //search by physical location
  void GetNeighborRanks(int myrank, int *neighbor_ranks);   // 80 neighbors
  int GetIndices(int rank, int &i, int &j, int&k, int&t);  // search by rank 
  int GetIndices(float x, float y, float z, float t, int& iidx, int &jidx, int&kidx, int& lidx);  // search by physical locations
  int GetBounds(int i, int j, int k, int t, volume_bounds_type& vb);  //return bound 
  int GetBounds(int rank, volume_bounds_type& vb);   // return bound 
  volume_bounds_type* GetBoundsList(){ return vb_list;}
  void GetLatticeDims(int& i, int&j, int &k, int &l) {i = idim; j=jdim; k=kdim; l = tdim; }
  bool isIn(float, float, float, float, int, int, int, int); 
  int CheckNeighbor(int myrank, float x, float y, float z, float t); 
  int GetNeighbor(int myrank, float x, float y, float z, float t, int &i, int &j, int &k, int&l);
  int GetProc(int, int, int, int); 
  int GetProc(int); 
  void InitSeedLists(); 
  void ResetSeedLists(); 
  void InsertSeed(int, int, int, int, VECTOR4); 
  void RoundRobin_proc(int n); 
  void GetPartitions(int, int**, int&); 
  void GetPartitions(int, int*, int&); 

  list<VECTOR4> *seedlists; 

 private: 

  int idim, jdim, kdim, tdim; //the lattice range
  int xdim, ydim, zdim, ldim; //the whole data range
  int npart; 
  volume_bounds_type *vb_list; 
  Partition4D *parts; // list of partition information

}; 

#endif 
