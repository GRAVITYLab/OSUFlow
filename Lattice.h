
#ifndef _LATTICE_H_
#define _LATTICE_H_ 
    
#include<stdio.h>
#include<stdlib.h>

#include "header.h"
#include "calc_subvolume.h"
#include "VectorMatrix.h" 

class  Lattice {

 public: 
  
  Lattice(int xlen, int ylen, int zlen, int ghost, int np); 
  ~Lattice(); 
  int GetRank(int i, int j, int k);        // serach by lattice indices
  int GetRank(float x, float y, float z);  //search by physical location
  void GetNeighborRanks(int myrank, int *neighbor_ranks);
  int GetIndices(int rank, int &i, int &j, int&k);  // search by rank 
  int GetIndices(float x, float y, float z, int& iidx, int &jidx, int&kidx);  // search by physical locations
  int GetBounds(int i, int j, int k, volume_bounds_type& vb);  //return bound 
  int GetBounds(int rank, volume_bounds_type& vb);   // return bound 
  volume_bounds_type* GetBoundsList(){ return vb_list;}
  void GetLatticeDims(int& i, int&j, int &k) {i = idim; j=jdim; k=kdim; }
  bool isIn(float, float, float, int, int, int); 
  int CheckNeighbor(int myrank, float x, float y, float z); 
  int GetNeighbor(int myrank, float x, float y, float z, int &i, int &j, int &k);
  void InitSeedLists(); 
  void ResetSeedLists(); 
  void InsertSeed(int, int, int, VECTOR3); 
  list<VECTOR3> *seedlists; 
 private: 
  int idim, jdim, kdim; //the lattice range
  int xdim, ydim, zdim; //the data range
  int* lattice; 
  int npart; 
  volume_bounds_type *vb_list; 
}; 

#endif 
