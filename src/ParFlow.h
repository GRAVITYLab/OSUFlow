//------------------------------------------------------------------------------
//
// parallel wrapper around OSUFlow
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
#ifndef _PAR_FLOW
#define _PAR_FLOW

#ifdef _MPI
#include <mpi.h>
#include "neighborhoods.hpp"
#include "assignment.hpp"
#include "blocking.hpp"
#endif

#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>
#include "OSUFlow.h"
#include "calc_subvolume.h"
#include "Blocks.h"
#include "Lattice4D.h"
#include "LatticeAMR.h"
#include <vector>

struct Item {
  float pt[4]; // position
  int steps; // number of steps
};
  
struct Particle {
  VECTOR4 pt; // particle position
  int steps; // number of steps this particle traveled so far
};
  
using namespace std;

class ParFlow {

 public:

#ifdef _MPI

  // edited TP 10/12/12
/*   ParFlow(Blocking *blocking, Assignment *assignment, Blocks *blocks,  */
/* 	  OSUFlow **osuflow, list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt, */
/* 	  int **npt, int *tot_ntrace, int nb, int track_seed_id = 0); */

  ParFlow(Blocks *blocks, 
	  OSUFlow **osuflow, list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt,
	  int **npt, int *tot_ntrace, int nb, int track_seed_id = 0);
#endif
  ParFlow(Lattice4D *lat, OSUFlow **osuflow, 
	  list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt, 
	  int **npt, int *tot_ntrace, int nb, int track_seed_id = 0);
  ParFlow(LatticeAMR *lat, OSUFlow **osuflow, 
	  list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt, 
	  int **npt, int *tot_ntrace, int nb, int track_seed_id = 0);
  ~ParFlow();
  void UpdateOSUFlow(OSUFlow **osuflow);
#if 0 // MOD-BY-LEETEN 01/17/2012-FROM:
  void ComputePathlines(vector<Particle> seeds, int block_num, int pf, 
			int end_steps, int *w = NULL);
  void ComputeStreamlines(vector<Particle> seeds, int block_num, int pf, 
			  int end_steps, int *w = NULL);
#else // MOD-BY-LEETEN 01/17/2012-TO:
  void ComputePathlines(const vector<Particle>& seeds, int block_num, int pf, 
			int end_steps, int *w = NULL);
  void ComputeStreamlines(const vector<Particle>& seeds, int block_num, int pf, 
			  int end_steps, int *w = NULL);
#endif // MOD-BY-LEETEN 01/17/2012-END
  void GatherFieldlines(int nblocks, float *size, int tsize);
  void SerialGatherFieldlines(int nblocks, float* size, int tsize);
  int GatherNumPts(int* &ntrace, int all, int nblocks);
  void GatherPts(int *ntrace, int mynpt, int nblocks);
  void PrintPerf(double TotTime, double TotInTime, double TotOutTime, 
		 double TotCompCommTime,
		 int TotParticles, float *size);
  void WriteFieldlines(int *ntrace, int mynpt, char *filename,
		       int nblocks, float *size, int tsize);
  void DistributedWriteFieldlines(char *filename, char *id_filename,
				  int nblocks, float *size, int tsize);
  void InitTraces(vector< vector <Particle> >& Seeds, int tf,
		  int nblocks, int tsize, int tblocks,
		    VECTOR3 *specific_seeds = NULL, 
		    int num_specific_seeds = 0);
  int ExchangeNeighbors(vector< vector<Particle> >& seeds, float wf);
  int SerExchangeNeighbors(vector< vector<Particle> >& seeds);
  int FlushNeighbors(vector< vector<Particle> >& seeds);
  double GetMyCompTime() { return comp_time; }
  void SetSeeds(OSUFlow* osuflow, float* from, float* to, 
		VECTOR3* specific_seeds, int num_specific_seeds, bool* isUsed);

  // seed list management
  void InitSeedLists4D(); 
  void InitSeedLists3D(); 
  void ResetSeedLists4D(); 
  void ResetSeedLists3D(); 
  void ResetSeedLists4D(int i) { seedlists4D[i].clear(); }
  void ResetSeedLists3D(int i) { seedlists3D[i].clear(); }
  bool InsertSeed(int i, int j, int k, int l, VECTOR4 p); 
  bool InsertSeed(int i, int j, int k, VECTOR3 p);
  bool InsertSeed(int from_i, int from_j, int from_k, int from_t, 
		  int to_i, int to_j, int to_k, int to_t, VECTOR4); 
  bool InsertSeed(int from_i, int from_j, int from_k, 
		  int to_i, int to_j, int to_k, VECTOR3); 
  bool InsertSeed(int to_rank, VECTOR4 p); 
  bool InsertSeed(int to_rank, VECTOR3 p); 
  bool InsertSeed(int from_rank, int to_rank, VECTOR4 p); 
  bool InsertSeed(int from_rank, int to_rank, VECTOR3 p); 
  void ResetFlowMatrix();
  int GetFlowMatrix(int i, int j);
  list<VECTOR3> *seedlists3D; 
  list<VECTOR4> *seedlists4D; 

  // integration parameters
  void SetMaxError(float maxError) {this->maxError = maxError;}
  void SetInitialStepSize(float step) {this->initialStepSize = step;}
  void SetMinStepSize(float step) {this->minStepSize = step;}
  void SetMaxStepSize(float step) {this->maxStepSize = step;}
  void SetLowerAngleAccuracy(float angle) {this->lowerAngleAccuracy = angle;}
  void SetUpperAngleAccuracy(float angle) {this->upperAngleAccuracy = angle;}
  void SetIntegrationOrder(INTEG_ORD order) {this->integrationOrder = order;}
  void SetUseAdaptiveStepSize(bool adapt) {this->useAdaptiveStepSize = adapt;}

  void SetIntegrationParams(OSUFlow* osuflow);

  int* flowMatrix; 

#ifdef ZOLTAN
  // wrapper around repartition method
  void Repartition(int grp, int *nblocks, VECTOR4 ***seeds, int type,
		   OSUFlow ***osuflow, Block *block, 
		   int compute_type, MPI_Comm comm, int *wgts = NULL);
#endif

 private:

  /* removed by TP 10/10/12 */
/*   void PointsToSeeds(vector< vector<Particle> >& seeds,  */
/* 		     vector<vector<char *> >points); */

  void PostPoint(int lid, Item *item, int recirc, int end_steps);

  int *block_stats; // block stats
  double *time_stats; // time stats
  int n_block_stats;// number of block stats
  int n_time_stats; // number of time stats
  int TotSeeds; // total number of seeds for all blocks and all rounds
  int TotSteps; // total number of integration steps for all seeds,
                // for all blocks and all rounds in this process
  int TotItemsSent; // total number of items sent
  OSUFlow **osuflow;
  Lattice4D *lat4D;
  LatticeAMR *latAMR;
  list<vtListTimeSeedTrace*> *sl_list;
  VECTOR4 **pt; // points in everyone's traces
  int **npt; // everyone's number of points in their traces
  int *tot_ntrace; // total number of everyone's traces
  int track_seed_id;
  int nb;

  // removed by TP 9/12/12
/* #ifdef _MPI */
/*   Blocking *blocking; // blocking object */
/*   Assignment *assign; // assignment object */
/*   Neighborhoods *nbhds; // neighborhoods object */
/* #endif */
// end TP 9/12/12

  Blocks *blocks; // blocks object
  double comp_time; // computation time for my process
  double comm_time1, comm_time2, comm_time3;

  // integration parameters
  float initialStepSize;
  float minStepSize;
  float maxStepSize;
  float maxError;
  float lowerAngleAccuracy;
  float upperAngleAccuracy;
  INTEG_ORD integrationOrder;
  bool useAdaptiveStepSize;

};

#ifdef _MPI
// the following are outside of the clasee so they can be used as callbacks
void CreateDtype(DIY_Datatype *dtype);
#endif

#endif
