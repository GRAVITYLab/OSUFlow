
#ifndef _LATTICE4D_H_
#define _LATTICE4D_H_ 
    
#include<stdio.h>
#include<stdlib.h>
#include <time.h>

#include "header.h"
#include "calc_subvolume.h"
#include "VectorMatrix.h" 
#include "OSUFlow.h"

#ifdef MPI
#include <mpi.h>
#endif

// maximum number of 4D neighbors, including myself
#define MAX_NEIGHBORS 81

// maximum number of blocks a process can have (todo: allocate dynamically)
#define MAX_BLOCKS 256

// a global (all processes) block or partition
struct Partition4D {
  int NumSendPoints[MAX_NEIGHBORS]; // number of points ready to send
  int SizeSendPoints[MAX_NEIGHBORS]; // size of sending points list (bytes)
  float *SendPoints[MAX_NEIGHBORS]; // sending points list
  int NumRecvPoints[MAX_NEIGHBORS]; // number of points received
  int SizeRecvPoints[MAX_NEIGHBORS]; // size of receiving points list (bytes)
  float *RecvPoints[MAX_NEIGHBORS]; // receiving points list
  int Proc; // process(or) number (mpi rank, core number, node number, etc.)
  int HasData; // data are ready
};

// a local block (for one process)
struct Block {
  int id;               // block ID (number)
  unsigned char status; // bitmap [msb                                  lsb]
                        //        [unused, ..., requested, loaded, computed]
  time_t time;           // last time accessed
};

class  Lattice4D {

 public: 
  
  Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp); 
  ~Lattice4D(); 
  int GetRank(int i, int j, int k, int l);         // serach by lattice indices
  int GetRank(float x, float y, float z, float t);  //search by physical location
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
  void ResetSeedLists(int i); 
  bool InsertSeed(int, int, int, int, VECTOR4); 
  bool InsertSeed(int from_i, int from_j, int from_k, int from_t, 
		  int to_i, int to_j, int to_k, int to_t, VECTOR4); 
  bool InsertSeed(int to_rank, VECTOR4); 
  bool InsertSeed(int from_rank, int to_rank, VECTOR4); 
  void RoundRobin_proc(int n); 
  void GetPartitions(int, int**, int&); 
  void ResetFlowMatrix(); // set all values to zero
  int GetFlowMatrix(int i, int j) {return flowMatrix[i*npart+j];}

  // block status manipulation
  void SetReq(int block_num) { my_blocks[block_num].status |= 0x04; }
  void ClearReq(int block_num) { my_blocks[block_num].status &= 0xf3; }
  int GetReq(int block_num) { return((my_blocks[block_num].status & 0x04) >> 2);}
  void SetLoad(int block_num) { my_blocks[block_num].status |= 0x02; }
  void ClearLoad(int block_num) { my_blocks[block_num].status &= 0xf5; }
  int GetLoad(int block_num) { return((my_blocks[block_num].status & 0x02) >> 1);}
  void SetComp(int block_num) { my_blocks[block_num].status |= 0x01; }
  void ClearComp(int block_num) { my_blocks[block_num].status &= 0xf6; }
  int GetComp(int block_num) { return(my_blocks[block_num].status & 0x01); }
  void SetTime(int block_num) { my_blocks[block_num].time = time(NULL); }
  time_t GetTime(int block_num) { return(my_blocks[block_num].time); }

  list<VECTOR4> *seedlists; 
  Block my_blocks[MAX_BLOCKS]; // list of blocks owned by my process

 private: 

  int idim, jdim, kdim, tdim; //the lattice range
  int xdim, ydim, zdim, ldim; //the whole data range
  int nbhd; // neighborhood size
  int npart; 
  volume_bounds_type *vb_list; 
  Partition4D *parts; // list of partition information
  int* flowMatrix; 

#ifdef MPI

 public:

  Lattice4D(int xlen, int ylen, int zlen, int tlen, int ghost, int nsp, int ntp, int d);
  int GetNumPartitions(int proc);
  void GetPartitions(int proc, int*p_list);
  void GetNeighborRanks(int myrank, int *neighbor_ranks);
  void NeighborIndices(int n, int i, int j, int k, int l, int &in,
		       int &jn, int &kn, int &ln);
  void PostPoint(int myrank, VECTOR4 p, int neighbor);
  void PrintPost(int myrank);
  void PrintRecv(int myrank);
  int GetNumRecv(int myrank);
  void GetRecvPts(int myrank, VECTOR4 *ls);
  void Error(const char *fmt, ...);
  void SendNeighbors(int myrank, MPI_Comm comm);
  int ReceiveNeighbors(int myrank, MPI_Comm comm);
  void SetData(int myrank, int has_data) {parts[myrank].HasData = has_data;}
  int GetData(int myrank) {return parts[myrank].HasData;}

#endif

}; 

#endif 
 
