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

// maximum number of 4D neighbors, including myself
#define MAX_NEIGHBORS 81

// maximum number of total global partitions
#define MAX_PARTS 100000000

// a global (all processes) block or partition
struct Partition4D {

  // the following arrays are indexed according to neighbor number
  int NumSendPoints[MAX_NEIGHBORS]; // number of points ready to send
  int SizeSendPoints[MAX_NEIGHBORS]; // size of sending points list (bytes)
  float *SendPoints[MAX_NEIGHBORS]; // sending points list
  int NumRecvPoints[MAX_NEIGHBORS]; // number of points received
  int SizeRecvPoints[MAX_NEIGHBORS]; // size of receiving points list (bytes)
  float *RecvPoints[MAX_NEIGHBORS]; // receiving points list

  // the following array fills up in order, from 0 to number of requests
  MPI_Request Reqs[4 * MAX_NEIGHBORS]; // message requests
  int NumReqs; // number of requests

  int Proc; // process(or) number (mpi rank, core number, node number, etc.)

  int Requed; // data are requested
  int Loaded; // data are loaded
  int Comped; // data are computed

};

class Partition {

 public: 
  
  Partition(int nsp, int ntp, int d); 
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
  void SendNeighbors(int myrank, int *ranks, MPI_Comm comm = MPI_COMM_WORLD);
  int ReceiveNeighbors(int myrank, int *ranks, MPI_Comm comm = MPI_COMM_WORLD);

  Partition4D *parts; // list of partition information

 private: 

  int nbhd; // neighborhood size
  int npart; 

/*   // todo: pointer to lattice base class, not one particular lattice */
/*   class Lattice4D *lat; // pointer to lattice */

}; 

#endif 
