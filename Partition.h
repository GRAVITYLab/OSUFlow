#ifndef _PARTITION_H_
#define _PARTITION_H_ 

#include<stdio.h>
#include<stdlib.h>
#include "VectorMatrix.h" 

#ifdef _MPI
#include <mpi.h>
#endif

// a global (all processes) block or partition
struct Partition4D {

  // the following arrays are indexed according to neighbor number
  int *NumSendPoints; // number of points ready to send
  int *SizeSendPoints; // size of sending points list (bytes)
  float **SendPoints; // sending points list
  int *NumRecvPoints; // number of points received
  int *SizeRecvPoints; // size of receiving points list (bytes)
  float **RecvPoints; // receiving points list

  int NumNeighbors; // the number of neighbors used
  int AllocNeighbors; // number of neighbors allocated

  // the following array fills up in order, from 0 to number of requests
#ifdef _MPI
  MPI_Request *SendReqs; // send requests
  MPI_Request *RecvReqs; // receive requests
  int NumSendReqs; // number of send requests
  int NumRecvReqs; // number of receive requests
#endif 

  int Proc; // process(or) number (mpi rank, core number, node number, etc.)

  int Requed; // data are requested
  int Loaded; // data are loaded
  int Comped; // data are computed

};

class Partition {

 public: 
  
  Partition(int npart, int nproc); 
  ~Partition(); 

  void SetReq(int myrank) { parts[myrank].Requed = 1; }
  void ClearReq(int myrank) { parts[myrank].Requed = 0; }
  int GetReq(int myrank) { return parts[myrank].Requed;}
  void SetLoad(int myrank) { parts[myrank].Loaded = 1; }
  void ClearLoad(int myrank) { parts[myrank].Loaded = 0; }
  int GetLoad(int myrank) { return parts[myrank].Loaded;}
  void SetComp(int myrank, int iter_num) { parts[myrank].Comped = iter_num; }
  void ClearComp(int myrank) { parts[myrank].Comped = -1; }
  int GetComp(int myrank, int iter_num) { return(parts[myrank].Comped >= 
						    iter_num); }
  void PostPoint(int myrank, VECTOR4 p, int neighbor);
  void PrintPost(int myrank);
  void PrintRecv(int myrank);
  void GetRecvPts(int myrank, VECTOR4 *ls);
  int GetProc(int myrank);
  int GetNumNeighbors(int myrank) { return parts[myrank].NumNeighbors; }
  int GetAllocNeighbors(int myrank) { return parts[myrank].AllocNeighbors; }
  void GrowNeighbors(int myrank);
  void AddNeighbor(int myrank, int neighblock, int neighrank);

#ifdef _MPI
  void ExchangeNeighbors(int **neighbor_ranks, VECTOR4 **seeds, int *size_seeds,
			MPI_Comm comm = MPI_COMM_WORLD);
#else
  int ExchangeNeighbors(int **neighbor_ranks, VECTOR4 **seeds, 
			int *size_seeds){};
#endif

  // global partition list
  Partition4D *parts;

  // data structures indexed by process for faster access when
  // aggregating messages from / to a process
  int **proc_parts; // partition ranks
  int *proc_nparts; // number of paritition ranks
  int ***proc_neighbors; // neighbor blocks
  int **proc_nneighbors; // number of neighbor blocks in proc_neighbors
  int **proc_aneighbors; // number of allocated neighbors in proc_neighbors

 private: 

  int npart; // total number of partitions
  int nproc; // total number of processes

}; 

#endif 
