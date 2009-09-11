#ifndef _PARTITION_H_
#define _PARTITION_H_ 

#include<stdio.h>
#include<stdlib.h>
#include "VectorMatrix.h" 

/* // todo: lattice base class header */
/* #include "Lattice4D.h" */

#ifdef MPI
#include <mpi.h>
#endif

// maximum number of total global partitions
#define MAX_PARTS 10000000

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
#ifdef MPI
  MPI_Request *Reqs; // message requests
  int NumReqs; // number of requests
#endif 

  int Proc; // process(or) number (mpi rank, core number, node number, etc.)

  int Requed; // data are requested
  int Loaded; // data are loaded
  int Comped; // data are computed

};

class Partition {

 public: 
  
  Partition(int nsp, int ntp); 
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
  void SetNumNeighbors(int myrank, int num) {parts[myrank].NumNeighbors = num; }
  int GetNumNeighbors(int myrank) { return parts[myrank].NumNeighbors; }
  int GetAllocNeighbors(int myrank) { return parts[myrank].AllocNeighbors; }
  void GrowNeighbors(int myrank);

#ifdef  MPI
  void SendNeighbors(int myrank, int *ranks, MPI_Comm comm = MPI_COMM_WORLD);
  int ReceiveNeighbors(int *block_ranks, int **neighbor_ranks,
		       int block, int nb, MPI_Comm comm = MPI_COMM_WORLD);
#else
  void SendNeighbors(int myrank, int *ranks){};
  int ReceiveNeighbors(int myrank, int *ranks, int **all_ranks){};
  int ReceiveNeighbors(int *block_ranks, int **neighbor_ranks,
		       int block, int nb){};
#endif

  Partition4D *parts; // list of partition information

  void Check(int myrank);

 private: 

  int npart; 

}; 

#endif 
