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

// ADD-BY-LEETEN 04/09/2011-BEGIN
#ifdef _MPI 
// ADD-BY-LEETEN 04/09/2011-END
#include <mpi.h>
#include "diy.h"
// ADD-BY-LEETEN 04/09/2011-BEGIN
#endif	// #ifdef _MPI 
// ADD-BY-LEETEN 04/09/2011-END
#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>
#include "OSUFlow.h"
#include "ParFlow.h"
#include "Partition.h"

// todo: replace 2D vector of seeds with pointer to it?
// or pass by reference (need these functions to modify seeds)

//-----------------------------------------------------------------------

#ifdef _MPI

//-----------------------------------------------------------------------
//
// parallel version
//

// edited TP 10/12/12
// ParFlow::ParFlow(Blocking *blocking, Assignment *assignment, Blocks *blocks, 
// 		 OSUFlow **osuflow, list<vtListTimeSeedTrace*> *sl_list, 
// 		 VECTOR4 **pt, int **npt, int *tot_ntrace, int nb, 
// 		 int track_seed_id) {

ParFlow::ParFlow(Blocks *blocks, 
		 OSUFlow **osuflow, list<vtListTimeSeedTrace*> *sl_list, 
		 VECTOR4 **pt, int **npt, int *tot_ntrace, int nb, 
		 int track_seed_id) {

  this->osuflow = osuflow;
  this->sl_list = sl_list;
  this->pt = pt;
  this->npt = npt;
  this->tot_ntrace = tot_ntrace;
  this->track_seed_id = track_seed_id;
  this->nb = nb;
  this->blocks = blocks;
  this->comm = MPI_COMM_WORLD; // Jimmy added for default comm.

  // deleted by TP 9/12/12  
//   this->blocking = blocking;
//   this->assign = assignment;
  // nbhds = new Neighborhoods(blocking, assign, MPI_COMM_WORLD);
  // end TP 9/12/12

  TotSeeds = 0;
  TotSteps = 0;
  TotItemsSent = 0;
  numStrms = 0;	    // added by Zhanping Liu on 05/30/2013 ZPL
  curvlinr = 0;     // added by Zhanping Liu on 07/11/2013 ZPL

  // integration parameters
  initialStepSize = 1.0;
  minStepSize = 0.01;
  maxStepSize = 5;
  maxError = 0.01;
  lowerAngleAccuracy = 3.0;
  upperAngleAccuracy = 15.0;
  integrationOrder = RK45;
  useAdaptiveStepSize = true;

  // performance stats
  n_block_stats = 5;
  n_time_stats = 4;
  block_stats = (int *)malloc(n_block_stats * sizeof(int));
  assert(block_stats != NULL);
  time_stats = (double *)malloc(n_time_stats * sizeof(double));
  assert(time_stats != NULL);

}
//----------------------------------------------------------------------------

#endif

//----------------------------------------------------------------------------
//
// serial lattice4D version
//
ParFlow::ParFlow(Lattice4D *lat, OSUFlow **osuflow, 
		 list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt, 
		 int **npt, int *tot_ntrace, int nb,
		 int track_seed_id) {

  lat4D = lat;
  latAMR = NULL;
  this->osuflow = osuflow;
  this->sl_list = sl_list;
  this->pt = pt;
  this->npt = npt;
  this->tot_ntrace = tot_ntrace;
  this->track_seed_id = track_seed_id;
  this->nb = nb;

  // deleted TP 9/12/12  
//   this->nbhds = NULL;

  TotSeeds = 0;
  TotSteps = 0;

  // performance stats
  n_block_stats = 5;
  n_time_stats = 6;
  block_stats = (int *)malloc(n_block_stats * sizeof(int));
  assert(block_stats != NULL);
  time_stats = (double *)malloc(n_time_stats * sizeof(double));
  assert(time_stats != NULL);

}
//----------------------------------------------------------------------------
//
// serial latticeAMR version
//
ParFlow::ParFlow(LatticeAMR *lat, OSUFlow **osuflow, 
		 list<vtListTimeSeedTrace*> *sl_list, VECTOR4 **pt, 
		 int **npt, int *tot_ntrace, int nb,
		 int track_seed_id) {

  latAMR = lat;
  lat4D = NULL;
  this->osuflow = osuflow;
  this->sl_list = sl_list;
  this->pt = pt;
  this->npt = npt;
  this->tot_ntrace = tot_ntrace;
  this->track_seed_id = track_seed_id;
  this->nb = nb;
  integrationDir = FORWARD_DIR;

  TotSeeds = 0;
  TotSteps = 0;

  // performance stats
  n_block_stats = 5;
  n_time_stats = 6;
  block_stats = (int *)malloc(n_block_stats * sizeof(int));
  assert(block_stats != NULL);
  time_stats = (double *)malloc(n_time_stats * sizeof(double));
  assert(time_stats != NULL);

}
//----------------------------------------------------------------------------
//
void ParFlow::UpdateOSUFlow(OSUFlow **osuflow) {

  this->osuflow = osuflow;

}
//-----------------------------------------------------------------------
//
// destructor
//
ParFlow::~ParFlow() {

  if (block_stats != NULL) {
    free(block_stats);
  }
  if (time_stats != NULL) {
    free(time_stats);
  }

  // deleted TP 9/12/12
//   if (nbhds != NULL) {
//     delete nbhds;
//   }
  // end TP 9/12/12

}
//-----------------------------------------------------------------------
//
// initializes traces and seeds them
// can randomly seed blocks or set specified seeds to blocks
//
// Seeds: generated seeds (output)
// tf: number of traces per block
// nblocks: local number of blocks
// tsize: number of timesteps
// tblocks: number of time blocks
// specific seeds: list of specified seeds
// num_specific_seeds: number of specified seeds
//
void ParFlow::InitTraces(vector< vector<Particle> >& Seeds, int tf,
			 int nblocks, int tsize, int tblocks, 
			 VECTOR3* specific_seeds,
			 int num_specific_seeds) {

  int i, j;
  float from[3], to[3]; // block spatial extent
  int min_t, max_t; // block temporal extent
  int64_t starts[4], sizes[4]; // block starts and sizes
  Particle seed; // one seed
  VECTOR3 *seeds; // 3D seeds for one block
  int nseeds; // number of seeds

  // if specific seeds are used, an array indicating if that seed has been
  // placed in a block or not. true if it has been placed, false otherwise.
  bool* isUsed = NULL;
  if(num_specific_seeds > 0)
  {
    isUsed = new bool[num_specific_seeds];
    for(i=0; i<num_specific_seeds; i++)
    {
      isUsed[i] = false;
    }
  }

  // init all blocks
  for (i = 0; i < nblocks; i++) {

#ifdef _MPI // parallel version
    blocks->ClearLoad(i);

    // edited TP 9/12/12
//     int64_t from64[4];
//     int64_t to64[4];
    float from[3];
    float to[3];
    bb_t bb;
//     blocks->GetRealBlockBounds(i, from, to);
    DIY_No_ghost_block_bounds(0, i, &bb);
//     from[0] = from64[0];
//     from[1] = from64[1];
//     from[2] = from64[2];
//     to[0] = to64[0];
//     to[1] = to64[1];
//     to[2] = to64[2];
//     int min_t = from[3];
//     int max_t = from[3];
    int min_t = bb.min[3];
    int max_t = bb.max[3]; // todo: why was this from instead of to above?

//     int gid = blocking->assign->RoundRobin_lid2gid(i);
    int gid = DIY_Gid(0, i);
    // end TP 9/12/12

#else // serial version
    if (lat4D) {
      lat4D->ClearLoad(i);
      lat4D->GetRealVB(i, from, to, &min_t, &max_t);
    }
    else {
      latAMR->ClearLoad(i);
      // TODO: get the bounds without ghost cells, this currently includes
      // ghost cells
      latAMR->GetVB(i, from, to, &min_t, &max_t);
    }
#endif

    // init seeds for blocks in first time group
    if (tsize == 1 || tblocks == 1 || min_t == 0) {

      if(num_specific_seeds > 0)
	// edited TP 10/12/12
// 	SetSeeds(osuflow[i], from, to, specific_seeds, num_specific_seeds, isUsed);
	#ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
	SetSeeds(osuflow[i], bb.min, bb.max, specific_seeds, num_specific_seeds, isUsed);
	// ADD-BY-LEETEN 10/29/2012-BEGIN
	#else // #ifdef _MPI
 	SetSeeds(osuflow[i], from, to, specific_seeds, num_specific_seeds, isUsed);
	#endif // #ifdef _MPI
	// ADD-BY-LEETEN 10/29/2012-END
      else
	// edited TP 10/12/12
// 	osuflow[i]->SetRandomSeedPoints(from, to, tf); 
	#ifdef _MPI	// ADD-BY-LEETEN 10/29/2012
	osuflow[i]->SetRandomSeedPoints(bb.min, bb.max, tf);
	// ADD-BY-LEETEN 10/29/2012-BEGIN
	#else // #ifdef _MPI
 	osuflow[i]->SetRandomSeedPoints(from, to, tf); 
	#endif // #ifdef _MPI
	// ADD-BY-LEETEN 10/29/2012-END
      seeds = osuflow[i]->GetSeeds(nseeds); 

      for (j = 0; j < nseeds; j++) {
	seed.pt.Set(seeds[j][0], seeds[j][1], seeds[j][2], min_t);
	seed.steps = 0;
	Seeds[i].push_back(seed);
      }

    }

  }

  if(isUsed != NULL)
  {
    delete[] isUsed;
  }

}
//-----------------------------------------------------------------------
//
// SetSeeds
//
// given a list of specific seeds to use, find the subset that is inside a
// given block, and set those seeds as the seeds for that block
//
//
void ParFlow::SetSeeds(OSUFlow* osuflow, float* from, float* to, 
                       VECTOR3* specific_seeds, int num_specific_seeds, 
                       bool* isUsed) {

  // first find the indices of the seeds which are in this block, then copy
  // those seeds over
  std::vector<int> indices;
  float x, y, z;
  for(int i=0; i<num_specific_seeds; i++)
  {
    if(!isUsed[i])
    {
      x = specific_seeds[i][0];
      y = specific_seeds[i][1];
      z = specific_seeds[i][2];
      if(x >= from[0] && x <= to[0] &&
	 y >= from[1] && y <= to[1] &&
	 z >= from[2] && z <= to[2])
      {
	indices.push_back(i);
	isUsed[i] = true;
      }
    }
  }

  int num_seeds = indices.size();
  VECTOR3* block_seeds = new VECTOR3[num_seeds]; // throws exception if fail
  for(int i=0; i<num_seeds; i++)
  {
    int j = indices[i];
    x = specific_seeds[j][0];
    y = specific_seeds[j][1];
    z = specific_seeds[j][2];
    block_seeds[i].Set(x, y, z);
  }

  osuflow->SetSeedPoints(block_seeds, num_seeds);
  delete [] block_seeds;

}
//-----------------------------------------------------------------------
//
// set all the integration parameters for an osuflow object
//
// osuflow: pointer to osuflow object
//
void ParFlow::SetIntegrationParams(OSUFlow* osuflow) {
  osuflow->SetMaxError(maxError);
  osuflow->SetInitialStepSize(initialStepSize);
  osuflow->SetMinStepSize(minStepSize);
  osuflow->SetMaxStepSize(maxStepSize);
  osuflow->SetLowerAngleAccuracy(lowerAngleAccuracy);
  osuflow->SetUpperAngleAccuracy(upperAngleAccuracy);
  osuflow->SetIntegrationOrder(integrationOrder);
  osuflow->SetUseAdaptiveStepSize(useAdaptiveStepSize);
}
//-----------------------------------------------------------------------
//
// computes streamlines for a block
//
// seeds: seeds for this block
// block_num: local block number (0 to nblocks-1)
// pf: point factor (points per trace)
// end_steps: number of steps a particle must travel before stopping
// w: weight of this block (output) (optional)
//
void ParFlow::ComputeStreamlines(const vector<Particle>& seeds, int block_num,
				 int pf, int end_steps, int *w) {

  list<vtListSeedTrace*> list3; // 3D list of traces
  std::list<VECTOR3*>::iterator pt_iter3; // 3D iterator over pts in one trace
  std::list<vtListSeedTrace*>::iterator trace_iter3; // 3D iter. over traces
  VECTOR3 p3; // 3D current point
  VECTOR4 *p = NULL; // 4D current point
  vtListTimeSeedTrace *trace; // 4D single trace
  int nseeds = seeds.size();

#ifdef _MPI
  comp_time = MPI_Wtime();
#endif

  if (nseeds > 0) {
	

    TotSeeds += nseeds;

    // make VECTOR3s (temporary)
    VECTOR3* temp_seeds = new VECTOR3[nseeds];
    for (int i = 0; i < nseeds; i++) {
      temp_seeds[i][0]= seeds[i].pt[0];
      temp_seeds[i][1]= seeds[i].pt[1];
      temp_seeds[i][2]= seeds[i].pt[2];
    }
	  
    // perform the integration
    SetIntegrationParams(osuflow[block_num]);
    osuflow[block_num]->GenStreamLines(temp_seeds, this->integrationDir, nseeds, pf, list3);

    // copy each 3D trace to a 4D trace and then to the streamline list
    // post end point of each trace to the send list
    int n = 0;
    for (trace_iter3 = list3.begin(); trace_iter3 != list3.end(); 
	 trace_iter3++)
    {

      TotSteps += (*trace_iter3)->size();
      if (w != NULL)
	*w += (*trace_iter3)->size(); // number of steps accrues to block weight

      // any  degenerate  'streamline' (with only one sample,  i.e.,  the  seed point) ZPL begin
      // MUST / SHOULD be discarded,  otherwise the  'current'  block A  would  manage
      // to forward its 'ending point' (actually the seed point, possibly still within
      // the REAL / non-ghost boundary of the 'current' block A) to the 'next' block B
      // which would not be able to locate the 'ending point' within its (even  ghost)
      // boundary and then would produce and forward a new 'ending point' (i.e., the
      // same as the seed point, due to the bug with the line above) back to block A
      //
      // flow line integration would be stuck by such endless end-point forwarding 
      // between neighborng blocks OVER CURVILINEAR GRIDS
      // 
      // the fix below discards degenerate streamlines such that only REAL end points
      // (different from the seed points) are forwarded to the 'next' block
      //
      // the fix also prevents single points (seed points) from being treated & hence
      // exported as streamlines, cleaning up the flow line integration result
      //
      // added by Zhanping Liu on 07/03/2013
      //
      if (    int(  ( *trace_iter3 )->size()  )    <    2    )    continue; 	    // ZPL end


      trace = new vtListTimeSeedTrace;
      for (pt_iter3 = (*trace_iter3)->begin(); pt_iter3 != 
	     (*trace_iter3)->end(); pt_iter3++) {
	p3 = **pt_iter3;
	p = new VECTOR4;
	p->Set(p3[0], p3[1], p3[2], 0.0f);
	trace->push_back(p);
      }

    	// find the matching seed for the trace
    	// (need the number of steps from the seed
    	p3 = **(*trace_iter3)->begin();
    	while ((seeds[n].pt[0] != p3[0] ||
    		  seeds[n].pt[1] != p3[1] ||
    		  seeds[n].pt[2] != p3[2])
    		  && n<nseeds ) // Jimmy modified: n stops when exceeding nseeds
    	  	  n++;
    	if (n == nseeds)
  		  	 fprintf(stderr, "Error: cannot find a match between seeds and list. "
  		  			 "This should not happen.\n");

    	// enqueue last point in the trace
    	Item item;
    	item.pt[0] = (*p)[0];
    	item.pt[1] = (*p)[1];
    	item.pt[2] = (*p)[2];
    	item.pt[3] = (*p)[3];
    	item.steps = seeds[n].steps + (*trace_iter3)->size();
    	PostPoint(block_num, &item, 0, end_steps);
    	sl_list[block_num].push_back(trace); // for later rendering
    	n++;

    }

    // cleanup
    delete[] temp_seeds;
    for(trace_iter3 = list3.begin(); trace_iter3 != list3.end();trace_iter3++)
    {
      for(pt_iter3 = (*trace_iter3)->begin(); pt_iter3 != (*trace_iter3)->end();
	  pt_iter3++) {
	delete *pt_iter3;
      }
      (*trace_iter3)->clear();
      delete *trace_iter3;
    }
    list3.clear();


  }
  
#ifdef _MPI
  comp_time = MPI_Wtime() - comp_time;
#endif

}
//-----------------------------------------------------------------------
//
// computes pathlines
//
// seeds: seeds for this block
// block_num: local block number (0 to nblocks-1)
// not global partition number
// pf: point factor (points per trace)
// end_steps: number of steps a particle must travel before stopping
// w: weight of this block (output) (optional)
//
void ParFlow::ComputePathlines(const vector<Particle>& seeds, int block_num,
			       int pf, int end_steps, int *w) {

  list<vtListTimeSeedTrace*> list; // list of traces
  std::list<VECTOR4*>::iterator pt_iter; // iterator over pts in one trace
  std::list<vtListTimeSeedTrace*>::iterator trace_iter; // iter. over traces
  VECTOR4 p; // current point
  int nseeds = seeds.size();
  int loaded; // whether block is loaded

#ifdef _MPI
  comp_time = MPI_Wtime();
  loaded = blocks->GetLoad(block_num);
#else
  if (lat4D)
    loaded = lat4D->GetLoad(block_num);
  else
    loaded = latAMR->GetLoad(block_num);
#endif

  if (loaded && nseeds > 0) {

    TotSeeds += nseeds;

    // make a temporary list of seeds
    VECTOR4 *temp_seeds = new VECTOR4[nseeds];
    for (int i = 0; i < nseeds; i++) {
      temp_seeds[i][0]= seeds[i].pt[0];
      temp_seeds[i][1]= seeds[i].pt[1];
      temp_seeds[i][2]= seeds[i].pt[2];
      temp_seeds[i][3]= seeds[i].pt[3];
    }

    // perform the integration
    SetIntegrationParams(osuflow[block_num]);
    osuflow[block_num]->GenPathLines(temp_seeds, list, FORWARD, nseeds, pf); 

    // copy each trace to the streamline list for later rendering
    // post end point of each trace to the send list
    int n = 0;
    for (trace_iter = list.begin(); trace_iter != list.end(); trace_iter++) {

      TotSteps += (*trace_iter)->size();
      if (w != NULL)
	*w += (*trace_iter)->size(); // number of steps accrues to block weight
      if (!(*trace_iter)->size())
	continue;

      // find the matching seed for the trace 
      // (need the number of steps from the seed)
      pt_iter = (*trace_iter)->begin();
      p = **pt_iter; // start point of computed trace
      while (n < nseeds && (seeds[n].pt[0] != p[0] || 
			    seeds[n].pt[1] != p[1] ||
			    seeds[n].pt[2] != p[2] || 
			    seeds[n].pt[3] != p[3]))
	n++;

      if (n == nseeds)
	fprintf(stderr, "Error: cannot find a match between seeds and list. "
		"This should not happen.\n");

      // get the end point
      pt_iter = (*trace_iter)->end();
      pt_iter--;
      p = **pt_iter;

      // enqueue last point in trace
      Item item;
      item.pt[0] = p[0];
      item.pt[1] = p[1];
      item.pt[2] = p[2];
      item.pt[3] = p[3];
      item.steps = seeds[n].steps + (*trace_iter)->size();
      PostPoint(block_num, &item, 0, end_steps);
      sl_list[block_num].push_back(*trace_iter); // for later rendering
      n++;

    }

    delete[] temp_seeds;

  }

  // recirculate seeds to myself if the block is not loaded yet
  if (!loaded) {

    for (int i = 0; i < nseeds; i++) {
      Item item;
      item.pt[0] = seeds[i].pt[0];
      item.pt[1] = seeds[i].pt[1];
      item.pt[2] = seeds[i].pt[2];
      item.pt[3] = seeds[i].pt[3];
      item.steps = seeds[i].steps;
      PostPoint(block_num, &item, 1, end_steps);  
    }
  }

#ifdef _MPI
  comp_time = MPI_Wtime() - comp_time;
#endif


}
//-----------------------------------------------------------------------
//
// gathers all fieldlines at the root for rendering
//
// nblocks: local number of blocks
// size: spatial data size
// tsize: temporal data size
//
void ParFlow::GatherFieldlines(int nblocks, float *size, int tsize) {
  
  static int *ntrace = NULL; // number of traces for each proc
  int n; // total number of my points

#ifdef GRAPHICS
  
  // gather number of points in each trace at the root
  n = GatherNumPts(ntrace, 0, nblocks);
  
  // gather the actual points in each trace at the root
  GatherPts(ntrace, n, nblocks);

  // write a file too
  WriteFieldlines(ntrace, n, (char *)"field_lines.out", nblocks, size, tsize);

#elif defined TRACK_SEED_ID
//#ifdef TRACK_SEED_ID
  
  DistributedWriteFieldlines((char *)"field_lines.out",
			     (char *)"field_line_ids.out", nblocks, size,
			     tsize);
#else
  // gather number of points in each trace to everyone
  n = GatherNumPts(ntrace, 1, nblocks);
   
  // write the traces collectively
  WriteFieldlines(ntrace, n, (char *)"field_lines.out", nblocks, size, tsize);
#endif

}
//-----------------------------------------------------------------------
//
// gathers all fieldlines for rendering, serial version
//
// nblocks: local number of blocks
// size: spatial data size
// tsize: temporal data size
//
void ParFlow::SerialGatherFieldlines(int nblocks, float* size, int tsize) {

  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  int i, j, k;
  int tot_npts = 0;

  // compute number of traces and points
  for (i = 0; i < nblocks; i++) {
    *tot_ntrace += sl_list[i].size();
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++)
      tot_npts += (*trace_iter)->size();
  }

  // allocate rendering data
  *npt = new int[*tot_ntrace]; // throws exception if fail
  *pt = new VECTOR4[tot_npts]; // points in everyones traces

  // compute number of points in each trace and collect the points
  j = 0;
  k = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
           pt_iter++)
	(*pt)[k++] = **pt_iter;
      (*npt)[j++] = (*trace_iter)->size();
    }
  }

  WriteFieldlines(tot_ntrace, tot_npts, (char *)"field_lines.out", nblocks,
		  size, tsize);
}

//-----------------------------------------------------------------------
//
// gathers number of points in each trace to the root
//
// ntrace: number of traces in each process (passed by reference)
// all: all = 0 gather to root, all = 1 gather to all
// nblocks: local number of blocks
//
// returns: total number of points in my process
//
int ParFlow::GatherNumPts(int* &ntrace, int all, int nblocks) {

  int myntrace = 0; // my number of traces
  static int *ofst = NULL; // offsets into ntrace
  int *mynpt; // number of points in each of my traces
  int tot_mynpt = 0; // total number of my points
  int rank, nproc; // MPI usual
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  int i, j;

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Comm_rank(this->comm, &rank);
  MPI_Comm_size(this->comm, &nproc);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END


  // allocate memory
  if (ntrace == NULL)
    ntrace = new int[nproc];
  if (ofst == NULL)
    ofst = new int[nproc];

  // compute number of my traces
  for (i = 0; i < nblocks; i++)
    myntrace += sl_list[i].size();

  numStrms = myntrace;	// added by Zhanping Liu on 05/30/2013 ZPL

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  // gather number of traces
  MPI_Allgather(&myntrace, 1, MPI_INT, ntrace, 1, MPI_INT, this->comm);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END


  // compute number of points in each of my traces
  mynpt = new int[myntrace];
  j = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      assert(j < myntrace); // sanity
      mynpt[j] = (*trace_iter)->size();
      tot_mynpt += mynpt[j++];
    }
  }

  // gather number of points in each trace
  *tot_ntrace = 0;
  for (i = 0; i < nproc; i++) {
    ofst[i] = (i == 0) ? 0 : ofst[i - 1] + ntrace[i - 1];
    *tot_ntrace += ntrace[i];
  }
  *npt = new int[*tot_ntrace];
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Allgatherv(mynpt, myntrace, MPI_INT, *npt, ntrace, ofst, MPI_INT,
		  this->comm);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END

  delete[] mynpt;

  return tot_mynpt;

}
//-----------------------------------------------------------------------
//
// gathers the points in each trace at the root
//
// ntrace: number of traces in each process
// mynpt: total number of points in my process
// nblocks: local number of blocks
//
void ParFlow::GatherPts(int *ntrace, int mynpt, int nblocks) {
  
  static int *nflt = NULL; // number of floats in points from each proc
  static int *ofst = NULL; // offsets into pt
  VECTOR4 *mypt; // points in my traces
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  int rank, nproc; // MPI usual
  int tot_npt = 0; // total number of points in all traces from everyone
  int i, j, k;

  // init
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Comm_rank(this->comm, &rank);
  MPI_Comm_size(this->comm, &nproc);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  if (nflt == NULL)
    nflt = new int[nproc];
  if (ofst == NULL)
    ofst = new int[nproc];
  mypt = new VECTOR4[mynpt];

  // collect my own points
  j = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
         pt_iter++)
	mypt[j++] = **pt_iter;
    }
  }

  // gather the points at the root
  if (rank == 0) 
  {
    	k = 0;
    	for (i = 0; i < nproc; i++) 
   	{
      		nflt[i] = 0;

      		for (j = 0; j < ntrace[i]; j++)
		nflt[i] += ((*npt)[k++] * 4);

      		ofst[i] = (i == 0) ? 0 : ofst[i - 1] + nflt[i - 1];
    	}
  }

  for(i = 0; i < *tot_ntrace; i++)
    tot_npt += (*npt)[i];

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  *pt = new VECTOR4[tot_npt];
  assert((*pt) != NULL);
  MPI_Gatherv(mypt, mynpt * 4, MPI_FLOAT, *pt, nflt, ofst,
	      MPI_FLOAT, 0, this->comm);

  		delete[] mypt;
		// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
}

// ADD-BY-LEETEN 04/09/2011-BEGIN
#ifdef _MPI 
// ADD-BY-LEETEN 04/09/2011-END
//-----------------------------------------------------------------------
//
// writes field lines keeping the header distributed
//
// filename: output file name
// id_filename: output id file name
// nblocks: local number of blocks
// size: spatial data size
// tsize: temporal data size
//
void ParFlow::DistributedWriteFieldlines(char *filename, char *id_filename,
				       int nblocks, float *size, int tsize) {

  MPI_File fd;
  MPI_Status status;
  int myproc;
  float *mypt; // points in my traces
  VECTOR4 temp; // temporary point
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  std::list<int64_t>::iterator trace_id_iter; // iterator over trace ids
  float min[4], max[4]; // extents
  int delim = -1; // delimits numbers of points from the points in the file
  int i, j, n;
    
  MPI_Comm_rank(this->comm, &myproc);
  
  int myntrace = 0; // my number of traces
  // compute number of my traces
  for (i = 0; i < nblocks; i++)    
    myntrace += sl_list[i].size();	// write numbers of points in each trace

  int totntrace = 0;
  MPI_Allreduce(&myntrace, &totntrace, 1, MPI_INT, MPI_SUM, this->comm);

  char filename2[1024];
  sprintf(filename2, "fieldlines.%d.out", totntrace);

  int stat = MPI_File_open(this->comm, filename2, MPI_MODE_CREATE |
		       MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
  assert(stat == MPI_SUCCESS);
  MPI_File_set_size(fd, 0); // start with an empty file every time

  if (myproc == 0) {

    // write extents
    min[0] = min[1] = min[2] = min[3] = 0.0;
    max[0] = size[0]; max[1] = size[1]; max[2] = size[2]; max[3] = tsize - 1;
    MPI_File_write(fd, min, 4, MPI_FLOAT, &status);
    assert(status.count == 4 * sizeof(float)); // bytes
    MPI_File_write(fd, max, 4, MPI_FLOAT, &status);
    assert(status.count == 4 * sizeof(float)); // bytes
	   
  }

  int *mynpt; // number of points in each of my traces
  int tot_mynpt = 0; // total number of my points

  // compute number of points in each of my traces  
  mynpt = (int *)malloc(sizeof(int) * myntrace);
  assert(mynpt != NULL);
  j = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end();
         trace_iter++) {
      assert(j < myntrace); // sanity
      mynpt[j] = (*trace_iter)->size();
      tot_mynpt += mynpt[j++];
    }
  } 

  // find my offset into the header
  int64_t myntrace64 = myntrace;
  int64_t header_offset = 0;
  MPI_Scan(&myntrace64, &header_offset, 1, MPI_LONG_LONG, 
	   MPI_SUM, this->comm);
  header_offset -= myntrace64;
  header_offset *= sizeof(int);
  header_offset += sizeof(float) * 8;

  stat = MPI_File_write_at_all(fd, header_offset, mynpt, myntrace, MPI_INT,
			       &status);
  assert(stat == MPI_SUCCESS);
  assert(status.count == myntrace * (int)(sizeof(int)));

  int64_t total_header_size;
  MPI_Allreduce(&myntrace64, &total_header_size, 1, MPI_LONG_LONG, MPI_SUM,
		  this->comm);

  total_header_size *= sizeof(int);
  total_header_size += sizeof(float) * 8;
  if (myproc == 0) {// write delimiter
    stat = MPI_File_write_at(fd, total_header_size, &delim, 1,
			     MPI_INT, &status);
    assert(stat == MPI_SUCCESS);
  }
  total_header_size += sizeof(int);
	
  // find point offset
  int64_t my_point_offset = 0;
  int64_t my_num_points = tot_mynpt;
  MPI_Scan(&my_num_points, &my_point_offset, 1, MPI_LONG_LONG, 
	   MPI_SUM, this->comm);
  my_point_offset -= my_num_points;
  my_point_offset *= sizeof(float) * 4;
  my_point_offset += total_header_size;

  // collect my points in a buffer
  mypt = (float *)malloc(tot_mynpt * 4 * sizeof(float));
  assert(mypt != NULL);
  n = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
	   pt_iter++) {
	temp = **pt_iter;
	mypt[4 * n]     = temp[0];
	mypt[4 * n + 1] = temp[1];
	mypt[4 * n + 2] = temp[2];
	mypt[4 * n + 3] = temp[3];
	n++;
      }
    }
  }

  // write my points
  stat = MPI_File_write_at_all(fd, my_point_offset, mypt, tot_mynpt * 4,
			       MPI_FLOAT, &status);
  assert(stat == MPI_SUCCESS);
  assert(status.count == tot_mynpt * 4 * (int)(sizeof(float))); // in bytes

  free(mypt);
  MPI_File_close(&fd);

#ifdef TRACK_SEED_ID
  if (id_filename != NULL) {

    char id_filename2[1024];
    sprintf(id_filename2, "fieldline_ids.%d.out", totntrace);
    stat = MPI_File_open(this->comm, id_filename2, MPI_MODE_CREATE |
			 MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
    assert(stat == MPI_SUCCESS);
    int64_t *my_trace_ids = (int64_t *)malloc(sizeof(int64_t) * myntrace);
    assert(my_trace_ids != NULL);
    int32_t which_id = 0;
    for (i = 0; i < nblocks; i++)
      {
	for (trace_id_iter = sl_id_list[i].begin(); 
	     trace_id_iter != sl_id_list[i].end(); trace_id_iter++, which_id++)
	  {
	    assert(which_id < myntrace);
	    my_trace_ids[which_id] = *trace_id_iter;
	  }
      }
    assert(which_id == myntrace);
    header_offset -= sizeof(float) * 8;
    header_offset /= sizeof(int);
    header_offset *= sizeof(int64_t);
    stat = MPI_File_write_at_all(fd, header_offset, my_trace_ids, myntrace,
				 MPI_LONG_LONG, &status);
    assert(stat == MPI_SUCCESS);
    MPI_File_close(&fd);
    free(my_trace_ids);
  }
#endif

}
// ADD-BY-LEETEN 04/09/2011-BEGIN
#endif	// #ifdef _MPI 
// ADD-BY-LEETEN 04/09/2011-END

//-----------------------------------------------------------------------
//
// writes field lines
//
// ntrace: number of traces in each process
// mynpt: total number of points in my process
// filename: output file name
// nblocks: local number of blocks
// size: spatial data size
// tsize: temporal data size
//
void ParFlow::WriteFieldlines(int *ntrace, int mynpt, char *filename,
	  		      int nblocks, float *size, int tsize) {

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_File fd;
  MPI_Status status;
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#else	// #ifdef _MPI 
  FILE *fd;
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  int myproc = 0;
  int nproc = 1;
  int64_t ofst; // offset into the file (bytes)
  int64_t pts_ofst = 0; // number of points before mine
  float *mypt; // points in my traces
  VECTOR4 temp; // temporary point
  std::list<vtListTimeSeedTrace *>::iterator trace_iter; // iterator over traces
  std::list<VECTOR4 *>::iterator pt_iter; // iterator over points in one trace
  float min[4], max[4]; // extents
  int delim = -1; // delimits numbers of points from the points in the file
  int i, j, n;
  int totntrace = 0;

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Comm_rank(this->comm, &myproc);
  MPI_Comm_size(this->comm, &nproc);

  int stat = MPI_File_open(this->comm, filename, MPI_MODE_CREATE |
		       MPI_MODE_WRONLY, MPI_INFO_NULL, &fd);
  assert(stat == MPI_SUCCESS);
  MPI_File_set_size(fd, 0); // start with an empty file every time
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#else	// #ifdef _MPI 
	fd = fopen(filename, "wb");
	assert(fd);
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
   if (myproc == 0) {

    // write extents (4D case)
    min[0] = min[1] = min[2] = min[3] = 0.0;
    max[0] = size[0]; max[1] = size[1]; max[2] = size[2]; max[3] = tsize - 1;

      // todo: diy amr case
#if 0
      latAMR->GetExtents(min, max);
#endif

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
    stat = MPI_File_write(fd, min, 4, MPI_FLOAT, &status);
    assert(stat == MPI_SUCCESS);
    stat = MPI_File_write(fd, max, 4, MPI_FLOAT, &status);
    assert(stat == MPI_SUCCESS);
	   
    // write numbers of points in each trace
    stat = MPI_File_write(fd, *npt, *tot_ntrace, MPI_INT,
			  &status);
    assert(stat == MPI_SUCCESS);
    assert(status.count == *tot_ntrace * (int)(sizeof(int)));

    // write delimiter
    stat = MPI_File_write(fd, &delim, 1, MPI_INT, &status);
    assert(stat == MPI_SUCCESS);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#else	// #ifdef _MPI 
    fwrite(min, sizeof(min[0]), 4, fd);
    fwrite(max, sizeof(max[0]), 4, fd);
    fwrite(*npt, 	sizeof(int), *tot_ntrace, fd);
    fwrite(&delim, 	sizeof(int), 1, fd);
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
   }

  // find total num traces
  for (i = 0; i < nproc; i++)
    totntrace += ntrace[i];

  // set file pointer to start of my points
  ofst = 8 * sizeof(float) + (*tot_ntrace + 1) * sizeof(int);
  n = 0;
  for (i = 0; i < myproc; i++) {
    for (j = 0; j < ntrace[i]; j++)
      pts_ofst += (*npt)[n++];
  }
  ofst += pts_ofst * 4 * sizeof(float); // pts before mine

  // collect my points in a buffer
  mypt = (float *)malloc(mynpt * 4 * sizeof(float));
  assert(mypt != NULL);
  n = 0;
  for (i = 0; i < nblocks; i++) {
    for (trace_iter = sl_list[i].begin(); trace_iter != sl_list[i].end(); 
         trace_iter++) {
      for (pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
	   pt_iter++) {
	temp = **pt_iter;
	mypt[4 * n]     = temp[0];
	mypt[4 * n + 1] = temp[1];
	mypt[4 * n + 2] = temp[2];
	mypt[4 * n + 3] = temp[3];
	n++;
      }
    }
  }

  // write my points
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_File_set_view(fd, ofst, MPI_FLOAT, MPI_FLOAT, (char *)"native", 
		    MPI_INFO_NULL);
  stat = MPI_File_write_all(fd, mypt, mynpt * 4, MPI_FLOAT, &status);
  assert(stat == MPI_SUCCESS);
  assert(status.count == mynpt * 4 * (int)(sizeof(float))); // in bytes
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#else	// #ifdef _MPI 
  fwrite(mypt, sizeof(float)*4, mynpt, fd);
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
 
  free(mypt);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_File_close(&fd);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#else	// #ifdef _MPI 
	fclose(fd);
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
}

//-----------------------------------------------------------------------
//
// todo: not using these for now
// need to diyitize them still; are they even being used?
//
//-----------------------------------------------------------------------

#if 0

void ParFlow::InitSeedLists4D() {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  seedlists4D = new list<VECTOR4>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists4D[i].clear(); 

}
//--------------------------------------------------------------------------
void ParFlow::InitSeedLists3D() {

  int npart = lat4D->npart; // total number of partitions

  seedlists3D = new list<VECTOR3>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists3D[i].clear(); 

}
//--------------------------------------------------------------------------
//
void ParFlow::ResetSeedLists4D() {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  for (int i = 0; i < npart; i++)
    seedlists4D[i].clear(); 

}
//--------------------------------------------------------------------------
//
void ParFlow::ResetSeedLists3D() {

  int npart = lat4D->npart; // total number of partitions

  for (int i = 0; i < npart; i++)
    seedlists3D[i].clear(); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int i, int j, int k, int t, VECTOR4 p) {

  int rank;
  if (lat4D != NULL)
    rank = lat4D->GetRank(i, j, k, t);
  else
    rank = latAMR->GetRank(i, j, k,t);

  if (rank ==-1) return(false); 

  seedlists4D[rank].push_back(p); 
  return(true); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int i, int j, int k, VECTOR3 p) {

  int rank = lat4D->GetRank(i, j, k, 0); // not using t

  if (rank ==-1) return(false); 

  seedlists3D[rank].push_back(p); 
  return(true); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int from_i, int from_j, int from_k, int from_t, 
			   int i, int j, int k, int t, VECTOR4 p) {

  int npart; // total number of partitions
  int from_rank;
  int to_rank;

  if (lat4D != NULL) {
    npart = lat4D->npart;
    from_rank = lat4D->GetRank(from_i, from_j, from_k, from_t); 
    to_rank = lat4D->GetRank(i,j,k,t); 
  } else {
    npart = latAMR->npart;
    from_rank = latAMR->GetRank(from_i, from_j, from_k, from_t); 
    to_rank = latAMR->GetRank(i,j,k, t); 
  }

  if (to_rank ==-1 || from_rank==-1) return(false); 

  seedlists4D[to_rank].push_back(p); 
  flowMatrix[from_rank*npart+to_rank]++; 
  return(true); 


}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int from_i, int from_j, int from_k,
			   int i, int j, int k, VECTOR3 p) {

  int npart; // total number of partitions
  int from_rank;
  int to_rank;

  npart = lat4D->npart;
  from_rank = lat4D->GetRank(from_i, from_j, from_k, 0); // not using t
  to_rank = lat4D->GetRank(i,j,k,0); // not using t

  if (to_rank ==-1 || from_rank==-1) return(false); 

  seedlists3D[to_rank].push_back(p); 
  flowMatrix[from_rank*npart+to_rank]++; 
  return(true); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int i, VECTOR4 p) {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  if (i>=npart) return(false); 

  seedlists4D[i].push_back(p); 
  return(true); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int i, VECTOR3 p) {

  int npart = lat4D->npart; // total number of partitions

  if (i>=npart) return(false); 

  seedlists3D[i].push_back(p); 
  return(true); 

}
//--------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int from_rank, int to_rank, VECTOR4 p) {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  if (from_rank >=npart || to_rank>=npart) return(false); 

  flowMatrix[from_rank*npart+to_rank]++; 
  seedlists4D[to_rank].push_back(p); 
  return(true); 

}
//---------------------------------------------------------------------------
//
bool ParFlow::InsertSeed(int from_rank, int to_rank, VECTOR3 p) {

  int npart = lat4D->npart; // total number of partitions

  if (from_rank >=npart || to_rank>=npart) return(false); 

  flowMatrix[from_rank*npart+to_rank]++; 
  seedlists3D[to_rank].push_back(p); 
  return(true); 

}
//---------------------------------------------------------------------------
//
void ParFlow::ResetFlowMatrix() {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  if (flowMatrix !=NULL) 
    memset(flowMatrix, '\0', npart*npart*sizeof(int)); 

}
//--------------------------------------------------------------------------
//
int ParFlow::GetFlowMatrix(int i, int j) {

  int npart; // total number of partitions

  if (lat4D != NULL)
    npart = lat4D->npart;
  else
    npart = latAMR->npart;

  return flowMatrix[i*npart+j];

}
//--------------------------------------------------------------------------

#endif

// MPI version of communication

#ifdef _MPI

//---------------------------------------------------------------------------
//
// exhanges points with all neighbors
//
// seeds: locations to store received points, indexed by local block number
// wf: wait factor (wait for this fraction of pending messages in each round)
//
// returns: total number of points received by this process
//
int ParFlow::ExchangeNeighbors(vector< vector<Particle> >& seeds, float wf) {

  // exchange points
  // edited TP 9/12/12
  void ***items = new void**[nb];
  int *num_items = new int[nb];
  int npr = 0;
//   int npr = nbhds->ExchangeNeighbors(pts, wf, &RecvItemDtype, &SendItemDtype);
  DIY_Exchange_neighbors(0, items, num_items, wf, &CreateDtype);

  // copy received points to seeds
  Particle seed; // one 4D seed
  for (int i = 0; i < nb; i++) { // for each block
    npr += num_items[i];
    seeds[i].clear();
    for (int j = 0; j < num_items[i]; j++) { // for each point in block
      seed.pt.Set(((Item *)items[i][j])->pt[0],
		  ((Item *)items[i][j])->pt[1],
		  ((Item *)items[i][j])->pt[2],
		  ((Item *)items[i][j])->pt[3]);
      seed.steps = ((Item *)items[i][j])->steps;
      seeds[i].push_back(seed);
    }
  }
  // end TP 9/12/12

  return npr;

}
//---------------------------------------------------------------------------
//
// completes neighbor exhange
//
// seeds: locations to store received points, indexed by local block number
//
// returns: total number of points received
//
int ParFlow::FlushNeighbors(vector< vector<Particle> >& seeds) {

  // edited TP 9/12/12
  void ***items = new void**[nb];
  int *num_items = new int[nb];
  int npr = 0;
//   int npr = nbhds->FlushNeighbors(pts, &RecvItemDtype);
  DIY_Flush_neighbors(0, items, num_items, &CreateDtype);

  // copy received points to seeds
  Particle seed; // one 4D seed
  for (int i = 0; i < nb; i++) { // for each block
    npr += num_items[i];
    // note that seeds are not cleared on the flush, so they carry over
    // to the next time group, ie, no seeds[i].clear() here
    for (int j = 0; j < num_items[i]; j++) { // for each point in block
      seed.pt.Set(((Item *)items[i][j])->pt[0],
		  ((Item *)items[i][j])->pt[1],
		  ((Item *)items[i][j])->pt[2],
		  ((Item *)items[i][j])->pt[3]);
      seed.steps = ((Item *)items[i][j])->steps;
      seeds[i].push_back(seed);
    }
  }
  // end TP 9/12/12

  return npr;

}
//------------------------------------------------------------------------
// //
// // DEPRECATED, removed by TP 10/10/12
// //
// //
// // copies exhanged points to seeds
// //
// // seeds: locations to store received points, indexed by local block number
// // pts: points for each block
// //
// void ParFlow::PointsToSeeds(vector< vector<Particle> >& seeds, 
// 			  vector< vector< char *> > points) {

//   Particle seed; // one 4D seed
//   for (int b = 0; b < points.size(); b++) { // each block
//     for (int p = 0; p < points[b].size(); p++) { // each point in the block
//       seed.pt.Set(((Item *)points[b][p])->pt[0],
// 		  ((Item *)points[b][p])->pt[1],
// 		  ((Item *)points[b][p])->pt[2],
// 		  ((Item *)points[b][p])->pt[3]);
//       seed.steps = ((Item *)points[b][p])->steps;
//       seeds[b].push_back(seed);
//     }
//   }

// }
//---------------------------------------------------------------------------

#endif

// serial version of communication


//---------------------------------------------------------------------------
//
// exhanges points with all neighbors
//
// seeds: locations to store received points, indexed by block number
//
// returns: total number of points exchanged
//
int ParFlow::SerExchangeNeighbors(vector< vector<Particle> >& seeds) {

  int r; // destination (global) block
  int np = 0; // total number of points
  int npart; // total number of partitions
  Partition4D *parts; // pointer to partition data structure
  int **neighbor_ranks; // ranks of neighbors of my blocks
  Particle seed; // one 4d seed

  // clear old seeds
  for(int i=0; i<seeds.size(); i++) {
    seeds[i].clear();
  }

  if (lat4D != NULL) {
    neighbor_ranks = lat4D->neighbor_ranks;
    npart = lat4D->npart;
    parts = lat4D->part->parts;
  }
  else {
    neighbor_ranks = latAMR->neighbor_ranks;
    npart = latAMR->npart;
    parts = latAMR->part->parts;
  }

  // for all global blocks
  for (int i = 0; i < npart; i++) {

    // j is (local) neighbor number in i's neighbor list
    for (int j = 0; j < parts[i].NumNeighbors; j++) {

      r = neighbor_ranks[i][j]; // global block id of neighbor

      // copy points to seeds
      for (int k = 0; k < parts[i].NumSendPoints[j]; k++) {
	seed.pt.Set(parts[i].SendPoints[j][4 * k],
		    parts[i].SendPoints[j][4 * k + 1],
		    parts[i].SendPoints[j][4 * k + 2],
		    parts[i].SendPoints[j][4 * k + 3]);
	seed.steps = 0; // not used for serial
	seeds[r].push_back(seed);
	np++;
      }

      parts[i].NumSendPoints[j] = 0;

    } // neighbor of partition rank

  } // global partition rank

  return np;

}
//------------------------------------------------------------------------

#ifdef ZOLTAN

//---------------------------------------------------------------------------

// todo: vectorize seeds

//
// wrapper for changing the partition
// assumes some existing paritition is in place already
//
// grp: time group
// nblocks: current number of my blocks (will be updated by Repartition)
// seeds: locations to store received points, indexed by local block number
// type: type of repartitioning algorithm used
//       0 = geometry-based (RCB)
//       1 = graph-based (not yet implemented)
// osuflow: pointer to array of osuflow objects for this process
// block: pointer to block class object
// comm: MPI communicator
// wgts: weights of blocks
//
void ParFlow::Repartition(int grp, int *nblocks, vector< vector<VECTOR4> >& seeds, 
			int type, 
			OSUFlow ***osuflow, Block *block, int compute_type,
			MPI_Comm comm, int *wgts) {

//   Partition *part; // pointer to partition data structure
//   int **block_ranks; // ranks of my blocks
//   int ***neighbor_ranks; // ranks of neighbors of my blocks
//   int ***neighbor_procs; // procs of neighbors of my blocks
//   int *avg_neigh; // average number of neighbors
//   int *alloc_blocks; // allocated number of blocks
//   int **alloc_neighbors; // allocated size of neighbor_ranks, neighbor_procs

//   if (lat4D != NULL) {
//     block_ranks = &lat4D->block_ranks;
//     neighbor_ranks = &lat4D->neighbor_ranks;
//     neighbor_procs = &lat4D->neighbor_procs;
//     avg_neigh = &lat4D->avg_neigh;
//     alloc_blocks = &lat4D->alloc_blocks;
//     alloc_neighbors = &lat4D->alloc_neighbors;
//     part = lat4D->part;
//   }
//   else {
//     block_ranks = &latAMR->block_ranks;
//     neighbor_ranks = &latAMR->neighbor_ranks;
//     neighbor_procs = &latAMR->neighbor_procs;
//     avg_neigh = &latAMR->avg_neigh;
//     alloc_blocks = &latAMR->alloc_blocks;
//     alloc_neighbors = &latAMR->alloc_neighbors;
//     part = latAMR->part;
//   }

//   // only geometry-based so far
//   assert(type == 0);

//   ChangePartition(grp, nblocks, block_ranks, neighbor_ranks, neighbor_procs, 
// 		  part, seeds, size_seeds, num_seeds, avg_neigh, alloc_blocks, 
// 		  alloc_neighbors, this->comm, osuflow, &AddNeighbor, wgts);

//   // update number of blocks
//   nb = *nblocks; //ParFlow's version
//   if (lat4D != NULL) // Lattice's version
//     lat4D->nb = *nblocks;
//   else
//     latAMR->nb = *nblocks;

//   // update ParFlow's version of seed info
//   this->Seeds = *seeds;
//   this->NumSeeds = *num_seeds;

//   // update ParFlow's and Block's version of osuflow info
//   UpdateOSUFlow(*osuflow);
//   block->UpdateCompute((void *)*osuflow, compute_type);

}
//---------------------------------------------------------------------------

#endif

//---------------------------------------------------------------------------
//
// prints performance stats
//
void ParFlow::PrintPerf(double TotTime, double TotInTime, double TotOutTime, 
		      double TotCompCommTime, int TotParticles, float *size) {

  int nproc;  // mpi groupsize
  int rank; // mpi rank
  int *all_block_stats; // gathered block stats
  double *all_time_stats; // gathered time stats
  float TotCells; // total number of spatial data cells in billions
  float TotDataSize; // total size of a time step in GB
  int i;

  int tot_npart = 0; // number of partitions per proc
  int min_npart = 0;
  int max_npart = 0;
  int mean_npart;
  float var_npart = 0.0;
  float std_npart;

  int tot_nneigh = 0; // number of neighbors 
  int min_nneigh = 0;
  int max_nneigh = 0;
  int mean_nneigh;
  float var_nneigh = 0.0;
  float std_nneigh;

  int tot_nseed = 0; // number of seeds
  int min_nseed = 0;
  int max_nseed = 0;
  int p_min_nseed = 0;
  int p_max_nseed = 0;
  int mean_nseed;
  float var_nseed = 0.0;
  float std_nseed;

  int tot_nstep = 0; // number of steps
  int min_nstep = 0;
  int max_nstep = 0;
  int p_min_nstep = 0;
  int p_max_nstep = 0;
  int mean_nstep;
  float var_nstep = 0.0;
  float std_nstep;

  int tot_nptsend = 0; // total number of points sent
  int min_nptsend = 0;
  int max_nptsend = 0;
  int p_min_nptsend = 0;
  int p_max_nptsend = 0;
  int mean_nptsend;
  float var_nptsend = 0.0;
  float std_nptsend;

  double tot_iotime = 0.0; // I/O time
  double min_iotime = 0;
  double max_iotime = 0;
  int p_min_iotime = 0;
  int p_max_iotime = 0;
  double mean_iotime;
  double var_iotime = 0.0;
  double std_iotime;

  double tot_commtime = 0.0; // communication time
  double min_commtime = 0;
  double max_commtime = 0;
  int p_min_commtime = 0;
  int p_max_commtime = 0;
  double mean_commtime;
  double var_commtime = 0.0;
  double std_commtime;

  double tot_comptime = 0.0; // computation time
  double min_comptime = 0;
  double max_comptime = 0;
  int p_min_comptime = 0;
  int p_max_comptime = 0;
  double mean_comptime;
  double var_comptime = 0.0;
  double std_comptime;

  double tot_outtime = 0.0; // output time
  double min_outtime = 0;
  double max_outtime = 0;
  int p_min_outtime = 0;
  int p_max_outtime = 0;
  double mean_outtime;
  double var_outtime = 0.0;
  double std_outtime;

	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Comm_rank(this->comm, &rank);
  MPI_Comm_size(this->comm, &nproc);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END

  // get stats of my process
#ifdef _MPI
  // edited TP 9/12/12
  block_stats[0] = nb; // number of blocks
  block_stats[3] = TotItemsSent; // total points sent
  // end TP 9/12/12
#else
  block_stats[0] = 0; // number of blocks  todo: fill in later
  block_stats[3] = 0; // total points sent todo: fill in later
#endif
  block_stats[1] = 0; // avg number of neighbors unused
  block_stats[2] = TotSeeds; // total number of particles advected
  block_stats[4] = TotSteps; // total number of steps advected
  time_stats[0] = TotInTime; // I/O time
  time_stats[1] = 0.0; // communication time (unused)
  time_stats[2] = 0.0; // computation time (unused)
  time_stats[3] = TotOutTime; // output time

  // alloc space and gather the stats
  all_block_stats = (int *)malloc(n_block_stats * nproc *
					  sizeof(int));
  assert(all_block_stats != NULL);
  all_time_stats = (double *)malloc(n_time_stats * nproc *
  					  sizeof(double));
  assert(all_time_stats != NULL);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END
  MPI_Gather(block_stats, n_block_stats, MPI_INT, all_block_stats, 
	     n_block_stats, MPI_INT, 0, this->comm);
  MPI_Gather(time_stats, n_time_stats, MPI_DOUBLE, all_time_stats, 
	     n_time_stats, MPI_DOUBLE, 0, this->comm);
	// ADD-BY-LEETEN 04/09/2011-BEGIN
	#endif	// #ifdef _MPI 
	// ADD-BY-LEETEN 04/09/2011-END

  // print the stats
  if (rank == 0) {

    // totals, mins, maxs
    for (i = 0; i < nproc; i++) {

      // total
      tot_npart    += all_block_stats[n_block_stats * i];
      tot_nneigh   += all_block_stats[n_block_stats * i + 1];
      tot_nseed    += all_block_stats[n_block_stats * i + 2];
      tot_nptsend  += all_block_stats[n_block_stats * i + 3];
      tot_nstep    += all_block_stats[n_block_stats * i + 4];
      tot_iotime   += all_time_stats[n_time_stats * i];
      tot_commtime += all_time_stats[n_time_stats * i + 1];
      tot_comptime += all_time_stats[n_time_stats * i + 2];
      tot_outtime  += all_time_stats[n_time_stats * i + 3];

      // min, max
      if (i == 0) {
	min_npart = max_npart       = all_block_stats[n_block_stats * i];
	min_nneigh = max_nneigh     = all_block_stats[n_block_stats * i + 1];
	min_nseed = max_nseed       = all_block_stats[n_block_stats * i + 2];
	min_nptsend = max_nptsend   = all_block_stats[n_block_stats * i + 3];
	min_nstep = max_nstep       = all_block_stats[n_block_stats * i + 4];
	min_iotime = max_iotime     = all_time_stats[n_time_stats * i];
	min_commtime = max_commtime = all_time_stats[n_time_stats * i + 1];
	min_comptime = max_comptime = all_time_stats[n_time_stats * i + 2];
	min_outtime = max_outtime   = all_time_stats[n_time_stats * i + 3];
	p_min_nseed = p_max_nseed = 0;
	p_min_nptsend = p_max_nptsend = 0;
	p_min_nstep = p_max_nstep = 0;
	p_min_iotime = p_max_iotime = 0;
	p_min_commtime = p_max_commtime = 0;
	p_min_comptime = p_max_comptime = 0;
	p_min_outtime = p_max_outtime = 0;
      }
      else {

	if (all_block_stats[n_block_stats * i] < min_npart)
	  min_npart = all_block_stats[n_block_stats * i];
	if (all_block_stats[n_block_stats * i] > max_npart)
	  max_npart = all_block_stats[n_block_stats * i];

	if (all_block_stats[n_block_stats * i + 1] < min_nneigh)
	  min_nneigh = all_block_stats[n_block_stats * i + 1];
	if (all_block_stats[n_block_stats * i + 1] > max_nneigh)
	  max_nneigh = all_block_stats[n_block_stats * i + 1];

	if (all_block_stats[n_block_stats * i + 2] < min_nseed) {
	  min_nseed = all_block_stats[n_block_stats * i + 2];
	  p_min_nseed = i;
	}
	if (all_block_stats[n_block_stats * i + 2] > max_nseed) {
	  max_nseed = all_block_stats[n_block_stats * i + 2];
	  p_max_nseed = i;
	}

	if (all_block_stats[n_block_stats * i + 3] < min_nptsend) {
	  min_nptsend = all_block_stats[n_block_stats * i + 3];
	  p_min_nptsend = i;
	}
	if (all_block_stats[n_block_stats * i + 3] > max_nptsend) {
	  max_nptsend = all_block_stats[n_block_stats * i + 3];
	  p_max_nptsend = i;
	}

	if (all_block_stats[n_block_stats * i + 4] < min_nstep) {
	  min_nstep = all_block_stats[n_block_stats * i + 4];
	  p_min_nstep = i;
	}
	if (all_block_stats[n_block_stats * i + 4] > max_nstep) {
	  max_nstep = all_block_stats[n_block_stats * i + 4];
	  p_max_nstep = i;
	}

	if (all_time_stats[n_time_stats * i] < min_iotime) {
	  min_iotime = all_time_stats[n_time_stats * i];
	  p_min_iotime = i;
	}
	if (all_time_stats[n_time_stats * i] > max_iotime) {
	  max_iotime = all_time_stats[n_time_stats * i];
	  p_max_iotime = i;
	}

	if (all_time_stats[n_time_stats * i + 1] < min_commtime) {
	  min_commtime = all_time_stats[n_time_stats * i + 1];
	  p_min_commtime = i;
	}
	if (all_time_stats[n_time_stats * i + 1] > max_commtime) {
	  max_commtime = all_time_stats[n_time_stats * i + 1];
	  p_max_commtime = i;
	}

	if (all_time_stats[n_time_stats * i + 2] < min_comptime) {
	  min_comptime = all_time_stats[n_time_stats * i + 2];
	  p_min_comptime = i;
	}
	if (all_time_stats[n_time_stats * i + 2] > max_comptime) {
	  max_comptime = all_time_stats[n_time_stats * i + 2];
	  p_max_comptime = i;
	}

	if (all_time_stats[n_time_stats * i + 3] < min_outtime) {
	  min_outtime = all_time_stats[n_time_stats * i + 3];
	  p_min_outtime = i;
	}
	if (all_time_stats[n_time_stats * i + 3] > max_outtime) {
	  max_outtime = all_time_stats[n_time_stats * i + 3];
	  p_max_outtime = i;
	}

      }
    }

    // means
    mean_npart = tot_npart / nproc;
    mean_nneigh = tot_nneigh / nproc;
    mean_nseed = tot_nseed / nproc;
    mean_nptsend = tot_nptsend / nproc;
    mean_nstep = tot_nstep / nproc;
    mean_iotime = tot_iotime / nproc;
    mean_commtime = tot_commtime / nproc;
    mean_comptime = tot_comptime / nproc;
    mean_outtime = tot_outtime / nproc;

    // variances
    for (i = 0; i < nproc; i++) {
      var_npart += (all_block_stats[n_block_stats * i] - mean_npart) *
	(all_block_stats[n_block_stats * i] - mean_npart);
      var_nneigh += (all_block_stats[n_block_stats * i + 1] - mean_nneigh) *
	(all_block_stats[n_block_stats * i + 1] - mean_nneigh);
      var_nseed += (all_block_stats[n_block_stats * i + 2] - mean_nseed) *
	(all_block_stats[n_block_stats * i + 2] - mean_nseed);
      var_nptsend += (all_block_stats[n_block_stats * i + 3] - mean_nptsend) *
	(all_block_stats[n_block_stats * i + 3] - mean_nptsend);
      var_nstep += (all_block_stats[n_block_stats * i + 4] - mean_nstep) *
	(all_block_stats[n_block_stats * i + 4] - mean_nstep);
      var_iotime += (all_time_stats[n_time_stats * i] - mean_iotime) *
	(all_time_stats[n_time_stats * i] - mean_iotime);
      var_commtime += (all_time_stats[n_time_stats * i + 1] - mean_commtime) *
	(all_time_stats[n_time_stats * i + 1] - mean_commtime);
      var_comptime += (all_time_stats[n_time_stats * i + 2] - mean_comptime) *
	(all_time_stats[n_time_stats * i + 2] - mean_comptime);
      var_outtime += (all_time_stats[n_time_stats * i + 3] - mean_outtime) *
	(all_time_stats[n_time_stats * i + 3] - mean_outtime);

    }
    var_npart /= nproc;
    var_nneigh /= nproc;
    var_nseed /= nproc;
    var_nptsend /= nproc;
    var_nstep /= nproc;
    var_iotime /= nproc;
    var_commtime /= nproc;
    var_comptime /= nproc;
    var_outtime /= nproc;

    // standard deviations
    std_npart = sqrt(var_npart);
    std_nneigh = sqrt(var_nneigh);
    std_nseed = sqrt(var_nseed);
    std_nptsend = sqrt(var_nptsend);
    std_nstep = sqrt(var_nstep);
    std_iotime = sqrt(var_iotime);
    std_commtime = sqrt(var_commtime);
    std_comptime = sqrt(var_comptime);
    std_outtime = sqrt(var_outtime);

    // misc: data size and aggregate bandwidth
    TotCells = size[0] * size[1] * size[2] / 1.0e6;
    TotDataSize = size[0] * size[1] * size[2] * 12 / 1048576;

    // print results
    fprintf(stderr, "------------------------- Performance Summary -------------------------\n");
    fprintf(stderr, "   -- Aggregate Values --\n");
    fprintf(stderr, "Number of procs = %d\n", nproc);
    fprintf(stderr, "Total time = %.2lf s\n", TotTime);
    fprintf(stderr, "Input time = %.2lf s\n", TotInTime);
    fprintf(stderr, "Compute + Communicate time = %.2lf s\n", TotCompCommTime);
    fprintf(stderr, "Output time = %.2lf s\n", TotOutTime);
    fprintf(stderr, "Total data size = %.3f million cells = %.2f MB\n", TotCells, TotDataSize);
    fprintf(stderr, "Total particles = %.3f thousand\n", TotParticles / 1.0e3);
    fprintf(stderr, "Total advection steps computed = %.3f million\n", tot_nstep / 1.0e6);
    fprintf(stderr, "Aggregate input I/O bandwidth = %.0lf MB/s\n", TotDataSize / TotInTime);
//     fprintf(stderr, "   -- Component Values For All Time Blocks --\n");
//     fprintf(stderr, "Comp time / proc (s) %7s min [%5d] = %-8.2lf max [%5d] = %-8.2lf avg = %-8.2lf std = %-8.2lf\n", "", p_min_comptime, min_comptime, p_max_comptime, max_comptime, mean_comptime, std_comptime);
//     fprintf(stderr, "Comm time / proc (s) %7s min [%5d] = %-8.2lf max [%5d] = %-8.2lf avg = %-8.2lf std = %-8.2lf\n", "", p_min_commtime, min_commtime, p_max_commtime, max_commtime, mean_commtime, std_commtime);
//     fprintf(stderr, "Output time / proc (s) %5s min [%5d] = %-8.2lf max [%5d] = %-8.2lf avg = %-8.2lf std = %-8.2lf\n", "", p_min_outtime, min_outtime, p_max_outtime, max_outtime, mean_outtime, std_outtime);
    fprintf(stderr, "Total pts comp / proc %6s min [%5d] = %-8d max [%5d] = %-8d avg = %-8d std = %-8.0f\n", "", p_min_nseed, min_nseed, p_max_nseed, max_nseed, mean_nseed, std_nseed);
    fprintf(stderr, "Total steps comp / proc %4s min [%5d] = %-8d max [%5d] = %-8d avg = %-8d std = %-8.0f\n", "", p_min_nstep, min_nstep, p_max_nstep, max_nstep, mean_nstep, std_nstep);
//     fprintf(stderr, "Total pts sent / proc %6s min [%5d] = %-8d max [%5d] = %-8d avg = %-8d std = %-8.0f\n", "", p_min_nptsend, min_nptsend, p_max_nptsend, max_nptsend, mean_nptsend, std_nptsend);
//     fprintf(stderr, "   -- Component Values For Final Time Block --\n");
//     fprintf(stderr, "Blocks / proc %14s min = %-8d max = %-8d avg = %-8d std = %-8.0f\n", "", min_npart, max_npart, mean_npart, std_npart);
//     fprintf(stderr, "Neighbors / block %10s min = %-8d max = %-8d avg = %-8d std = %-8.0f\n", "", min_nneigh, max_nneigh, mean_nneigh, std_nneigh);
//     fprintf(stderr, "-----------------------------------------------------------------------\n");

    free(all_block_stats);
    free(all_time_stats);
  } // rank = 0

}
//-----------------------------------------------------------------------
//
// posts a point for sending
//
// lid: local block id of the current block
// item: point to be sent
// recirc: whether to recirculate a point that does not leave the current block
// end_steps: number of steps a particle must go before it should stop
//
void ParFlow::PostPoint(int lid, Item *item, int recirc, int end_steps) {

#ifdef _MPI


  // NOTE: due to the current  'back-off'  problem with  OSUFlow,  the ending ZPL begin
  //       point  of  a flow line (still within the 'current' block) fails to
  //       be   forwarded   to   the  'next'  block  FOR  CURVILINEAR  GRIDS:
  //
  //       ( i == 4 ) always holds!!!
  //
  //       it is what I call a 'block-bound' point:  a point bound to a block
  //
  //       once the 'back-off' problem is fixed,  the ending point  will be a
  //       little bit  beyond the REAL (non-ghost) boundary  of the 'current'
  //       block --- into the ghost cell (that is shared by the 'next' block)
  //       such that the 'next' block  will  successfully  catch  and  accept
  //       such an ending point (as we expect very much)
  //
  //       above commented by Zhanping Liu on 07/03/2013 for CURVILINEAR GRID
  //
  // the if-statement code segment below is executed ONLY IF  NOT circulation
  // and NOT curvilinear
  // 
  // this  combined  logical  condition  supports curvilinear grids (with the
  // 'back-off' problem  ALLEVIATED)  without  affecting  the  original  code 
  // (NON-curvilinear grids) AND without any performance penalty AT ALL
  //
  // curvlinr: an  instance variable of the class  initialized  to  0  (MUST)
  //           while  it can be explicitly specified (by the user of ParFlow)
  //           to 1 through a public member function: SetGrid2Curvilinear(  )
  //
  // combined logical condition added by Zhanping Liu on 07/11/2013
  //
  int  toExecut = !( recirc + curvlinr );                                  // ZPL end


  // start edited TP 9/12/12
  if ( toExecut )  // modified (from a single condition to the current combined one) by Zhanping Liu on 07/11/2013  ZPL
  {   	
	// the 3 lines below were moved from outside in, by Zhanping Liu on 07/11/2013 ZPL begin
	// without any influence on the flow of execution
        //
	int gid = DIY_Gid(0, lid);                   
  	bb_t bounds;
  	DIY_No_ghost_block_bounds(0, lid, &bounds);                                 // ZPL end

  	// check if point is in lid and not recirculating
      	int i;
      	for (i = 0; i < 4; i++) 
      	if (item->pt[i] < bounds.min[i] || item->pt[i] > bounds.max[i]) break;

      	if (i == 4) return;
  }
  
  if (item->steps < end_steps)
    DIY_Enqueue_item_points(0, lid, item, NULL, sizeof(Item), item->pt, 
			    1, NULL);
  // end TP

#else

  lat4D ? lat4D->PostPoint(lid, item->pt, recirc) : 
    latAMR->PostPoint(lid, item->pt, recirc);

#endif

}
//-----------------------------------------------------------------------

#ifdef _MPI

//-----------------------------------------------------------------------
//
// creates a DIY datatype for an item
//
void CreateDtype(DIY_Datatype *dtype) {

  MPI_Datatype types[2];
  int lengths[2];
  MPI_Aint displs[2];

  types[0] = MPI_FLOAT;
  displs[0] = offsetof(struct Item, pt);
  lengths[0] = 4;

  types[1] = MPI_INT;
  displs[1] = offsetof(struct Item, steps);
  lengths[1] = 1;

  MPI_Type_create_struct(2, lengths, displs, types, dtype);

}
//-----------------------------------------------------------------------

#endif
