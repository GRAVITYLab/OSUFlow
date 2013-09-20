//------------------------------------------------------------------------------
//
// block class
// initializes, seeds, loads and unloads data blocks
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

#ifdef _MPI
#include <mpi.h>
#include "diy.h"
#include "bil.h"
#endif

// ADD-BY-LEETEN 12/16/2011-BEGIN
#if WITH_UNISTD
#include <unistd.h>
#endif // #if WITH_UNISTD
// ADD-BY-LEETEN 12/16/2011-END

#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>
#include "Blocks.h"

// pointers to computation objects
// need to be outside of class when ifdef'd
// valgrind complains otherwise-not sure why
#ifdef _OSUFLOW
  OSUFlow **block_osuflow;
#endif

//----------------------------------------------------------------------------

#ifdef _MPI

//----------------------------------------------------------------------------
//
// parallel version
//
// changed by TP 10/12/12
// Blocks::Blocks(Blocking *blocking, Assignment *assignment, void *compute, 
// 	       int compute_type, char **dataset_files, int num_dataset_files,
// 	       DataMode data_mode) {
Blocks::Blocks(int nblocks, void *compute, int compute_type, 
	       char **dataset_files, int num_dataset_files, 
	       DataMode data_mode) {

  // deleted by TP 10/12/12
//   this->ghost = ghost;

  //strncpy(this->filename, filename, sizeof(this->filename));
  this->data_mode = data_mode;
  this->dataset_files = dataset_files;
  this->num_dataset_files = num_dataset_files;
  this->compute_type = compute_type;

  // editedTP 10/12/12
//   this->blocking = blocking;
//   assign = assignment;
//   int nb = assign->NumBlks();
  int nb = nblocks;
  // end TP

  blocks.reserve(nb);
  for (int i = 0; i < nb; i++) {
    lb_t lb;
    // changed by TP 10/9/12
//     lb.gid = assign->RoundRobin_lid2gid(i);
    lb.gid = DIY_Gid(0, i);
    lb.loaded = 0;
    blocks.push_back(lb);
  }

  switch (compute_type) {
#ifdef _OSUFLOW
  case OSUFLOW:
    block_osuflow = (OSUFlow **)compute;
    break;
#endif
  default:
    break;
  }

}
//----------------------------------------------------------------------------

#endif

//----------------------------------------------------------------------------
//
// serial 4D version
//
Blocks::Blocks(Lattice4D *lat, void *compute, int compute_type, 
	       char **dataset_files, int num_dataset_files,
	       DataMode data_mode) {

  lat4D = lat;
  latAMR = NULL;

  // deleted TP 10/12/12
//   ghost = 0;

  this->data_mode = data_mode;
  this->dataset_files = dataset_files;
  this->num_dataset_files = num_dataset_files;
  this->compute_type = compute_type;

  switch (compute_type) {
#ifdef _OSUFLOW
  case OSUFLOW:
    block_osuflow = (OSUFlow **)compute;
    break;
#endif
  default:
    break;
  }

}
//----------------------------------------------------------------------------
//
// serial AMR version
//
Blocks::Blocks(LatticeAMR *lat, void *compute, int compute_type, 
	       char **dataset_files, int num_dataset_files,
	       DataMode data_mode) {

  latAMR = lat;
  lat4D = NULL;

  // deleted 10/12/12
//   ghost = 0;

  this->data_mode = data_mode;
  this->dataset_files = dataset_files;
  this->num_dataset_files = num_dataset_files;
  this->compute_type = compute_type;

  switch (compute_type) {
#ifdef _OSUFLOW
  case OSUFLOW:
    block_osuflow = (OSUFlow **)compute;
    break;
#endif
  default:
    break;
  }

}
//----------------------------------------------------------------------------
//
void Blocks::UpdateCompute(void *compute, int compute_type) {

  switch (compute_type) {
#ifdef _OSUFLOW
  case OSUFLOW:
    block_osuflow = (OSUFlow **)compute;
    break;
#endif
  default:
    break;
  }

}
//----------------------------------------------------------------------------
//
Blocks::~Blocks() {

}
//--------------------------------------------------------------------------

#ifdef _MPI
#ifdef USE_BIL

//--------------------------------------------------------------------------
//
// Only used if USE_BIL is defined. This function calls the bil library
// on groups of blocks from timesteps and then returns the data in a
// float array
//
float ***Blocks::BilLoadTimeGroupBlocks(int t_group, int nblocks,
					float *size, int tsize, int tb) 
{
  	MPI_Datatype 	bil_datatype;
  	MPI_Type_contiguous(3, MPI_FLOAT, &bil_datatype );
  	float     ***   data = new float ** [nblocks];

  	// edited TP 10/12/12
	// int64_t block_min[4], block_size[4]; // block extents
  	int block_min[4], block_size[4]; // block extents

  	// account for a header if necessary
  	// DIY does not provide this function, use BIL directly for it
  	if(data_mode == RAW_HEADER)
    	BIL_Set_io_header_size(12);  // header is 3 ints

  	// load blocks for this time block
  	for (int i = 0; i < nblocks; i++) 
	{ 	
    		// edited TP 10/12/12
    		int time_block;
    		DIY_In_time_block(0, i, &time_block);
		//     if (blocking->InTimeBlock(t_group, i, tsize, tb)) {
    		if (time_block == t_group) 
		{
			// blocking->BlockStartsSizes(i, block_min, block_size);
      			DIY_Block_starts_sizes(0, i, block_min, block_size);
      			// end TP

      			data[i] = new float*[block_size[3]];

      			// edited TP 10/12/12
			//       // convert from int64_t to int, and also reverse order for BIL
			//       int bil_data_size[3] = { size[2], size[1], size[0] };
			//       int bil_min[3] = { block_min[2], block_min[1], block_min[0] };
			//       int bil_size[3] = { block_size[2], block_size[1], block_size[0] };
			//       // post a BIL read for each time step in this block
			//       for (int j = 0; j < block_size[3]; j++) { // for all timesteps
			// 	data[i][j] = new float[block_size[0] * block_size[1] * 
			// 			       block_size[2] * 3];
			// 	BIL_Add_block_raw(3, bil_data_size, bil_min, bil_size, 
			// 			  dataset_files[block_min[3] + j], bil_datatype, 
			// 			  (void **)&data[i][j]);

      			// post a read for each time step in this block
      			for (int j = 0; j < block_size[3]; j++) 
			{ 	// for all timesteps
				data[i][j] = new float[block_size[0] * block_size[1] * block_size[2] * 3];
				DIY_Add_data_raw(  block_min,  block_size, 
			  	                   dataset_files[ block_min[3] + j ], bil_datatype, (void **)&data[i][j]);
      			} // for all timesteps
    		} // in time block
  	} // for all blocks

  	// edited TP 10/12/12
	//   BIL_Read();
  	DIY_Read_data_all();

  	// swap bytes
	#ifdef BYTE_SWAP
  	for (int i = 0; i < nblocks; i++) 
	{ 	// for all my blocks
    		// edited TP 10/12/12
    		DIY_In_time_block(i, &time_block)
		//     if (blocking->InTimeBlock(t_group, i, tsize, tb)) {
    		if (time_block == t_group) 
		{
			// blocking->BlockStartsSizes(i, block_min, block_size);
      			DIY_Block_starts_sizes(i, block_min, block_size);
      			// end TP

      			for (int j = 0; j < block_size[3]; j++) 
			{
				for (int k = 0; k < block_size[0] * block_size[1] * block_size[2] * 3; k++)
	  			swap4((char *)&data[i][j][k]);
      			}
    		}
  	}
	#endif

  	MPI_Type_free(&bil_datatype);

  	return data;
}
//-----------------------------------------------------------------------

#endif
#endif

//-----------------------------------------------------------------------
//
// loads blocks in this time group
//
// grp: current group
// time (output): the time it took to load
// nblocks: number of local blocks
// size: spatial domain size
// tsize: number of time steps in dataset
// tb: number of global time blocks
// data: data that has already been read in (defaults to NULL)
//
// returns: number of bytes actually loaded
//
int Blocks::LoadBlocks4D(int grp, double *time, int nblocks,
			 float *size, int tsize, int tb, float ***data) {

  int s = 0; // data size (bytes)

  // edited TP 10/12/12
//   int64_t starts[4], sizes[4]; // block starts and sizes
  int starts[4], sizes[4]; // block starts and sizes

  int min_t, max_t; // block temporal extent
  int loaded; // whether a block is loaded
  double t0;
  int i;

  *time = 0.0;

  for (i = 0; i < nblocks; i++) {

#ifdef _MPI
    loaded = GetLoad(i);
#else
    loaded = lat4D->GetLoad(i);
#endif

    #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
    // edited by TP 10/10/20
    int time_block;
    DIY_In_time_block(0, i, &time_block);
//     if (!loaded && IsBlockInTimeGroup(grp, i, tsize, tb)) {
    if (!loaded && time_block == grp) {
      // end TP
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

      // compute block extents
#ifdef _MPI // parallel version

      // edited TP 10/12/12
//       blocking->BlockStartsSizes(i, starts, sizes);
      DIY_Block_starts_sizes(0, i, starts, sizes);

      float from[3] = {starts[0], starts[1], starts[2]};
      float to[3] = {starts[0] + sizes[0] - 1, starts[1] + sizes[1] - 1,
		     starts[2] + sizes[2] - 1};
      min_t = starts[3];
      max_t = starts[3] + sizes[3] - 1;

      // get the real bounds (not including ghost cells)
      int real_from[4];
      int real_to[4];

      // deleted by TP 10/10/12 (gid not used later)
      // int gid = blocking->assign->RoundRobin_lid2gid(i); // global id of block

//       for(int j=0; j<4; j++)
//       {
// 	int64_t real_from64[4];
// 	int64_t real_to64[4];
// 	blocking->GetRealBlockBounds(i, real_from64, real_to64);
// 	real_from[0] = real_from64[0];
// 	real_from[1] = real_from64[1];
// 	real_from[2] = real_from64[2];
// 	real_from[3] = real_from64[3];
// 	real_to[0] = real_to64[0];
// 	real_to[1] = real_to64[1];
// 	real_to[2] = real_to64[2];
// 	real_to[3] = real_to64[3];
//       }

      struct bb_t bb;
      DIY_No_ghost_block_bounds(0, i, &bb);
      real_from[0] = bb.min[0];
      real_from[1] = bb.min[1];
      real_from[2] = bb.min[2];
      real_from[3] = bb.min[3];
      real_to[0] = bb.max[0];
      real_to[1] = bb.max[1];
      real_to[2] = bb.max[2];
      real_to[3] = bb.max[3];

      // end TP

      t0 = MPI_Wtime();
#else // serial version
      float from[3], to[3]; // bounds with ghost
      float r_from[3], r_to[3];  // bounds w/o ghost, in float format
      int real_from[4];  // bounds w/o ghost, int format
      int real_to[4];    // bounds w/o ghost, int format

      lat4D->GetVB(i, from, to, &min_t, &max_t);
      lat4D->GetRealVB(i, r_from, r_to, &real_from[3], &real_to[3]);

      // convert from float to int
      real_from[0] = r_from[0];
      real_from[1] = r_from[1];
      real_from[2] = r_from[2];
      real_to[0] = r_to[0];
      real_to[1] = r_to[1];
      real_to[2] = r_to[2];
#endif

      // load the block
      switch (compute_type) {
      case OSUFLOW:
	if(data == NULL)
	  block_osuflow[i]->LoadData(dataset_files, num_dataset_files, from, 
				     to, real_from, real_to, size, min_t,
				     max_t, data_mode);
	else
	  block_osuflow[i]->LoadData(dataset_files, num_dataset_files, from, 
				     to, real_from, real_to, size, min_t,
				     max_t, data_mode, data[i]); 	
	break;
      default:
	break;
      }

#ifdef _MPI
      *time += (MPI_Wtime() - t0);
      SetLoad(i);
#else
      lat4D->SetLoad(i);
#endif

      s += ((to[0] - from[0]) * (to[1] - from[1]) *
	    (to[2] - from[2]) * 3 * sizeof(float));

    #ifdef _MPI			// ADD-BY-LEETEN 10/29/2012
    }
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

  }

  return s;

}

#ifdef _MPI // compatible for builds without MPI
// ======================================================================= // ZPL begin
//                                                                         //
//                         below are the functions                         //
//                                                                         //
//                              in support of                              //
//                                                                         //
//                    integrating  VTK I/O with OSUFlow                    //
//                                                                         //
//             added by Zhanping Liu (05/23/2013 ~ 07/08/2013)             //
//                                                                         //
// ======================================================================= //


// -----------------------------------------------------------------------
//
// attach DIY data blocks (distributed to THIS process from process #0) to
// the OSUFlow engine
//
// time_grp: index of the time group (or time block)
// timeUsed: time used for this attachment process
// numBlcks: number of spatial-temporal blocks of THIS process
//           (NOTE: in the time-varying case, NOT all these blocks are ACTIVE)
// volSizes: size of the WHOLE volume (NOT that of any individual block)
// dataBlks: [numBlcks][timeGroupSize][blockSize]
// 
// added by Zhanping Liu on 05/26/2013 and last updated on 06/04/2013
//
int 	Blocks::AttachDataBlocks4D
	( int time_grp, double * timeUsed, 
	  int numBlcks, float  * volSizes, float *** dataBlks )
{	
	cout << endl << "==== attaching DIY data blocks (time group #" 
	             << time_grp << ") to the OSUFlow engine ====" << endl;

	int  i, tmGrpIdx;
	int  blockMin[4], blockSiz[4];


	#ifdef BYTE_SWAP
  	for ( i = 0; i < numBlcks; i ++ )
	{ 
    		DIY_In_time_block( 0, i, &tmGrpIdx );
    		if ( tmGrpIdx != time_grp ) continue;

      		DIY_Block_starts_sizes( 0, i, blockMin, blockSiz );
      		for ( int j = 0; j < blockSiz[3]; j ++ ) 
		for ( int k = 0; k < blockSiz[0] * blockSiz[1] * blockSiz[2] * 3; k ++ )
	  	swap4(  ( char * ) &dataBlks[i][j][k]  );
  	}
	#endif


 	int 	dataSize = 0;
  	int 	beLoaded;
	int	minTmStp;
	int	maxTmStp;
	int     real_min[4];
      	int     real_max[4];
	float   spaceMin[3];
	float   spaceMax[3];
  	double 	start_tm;
	struct  bb_t real_bb4;

  	*timeUsed = 0.0;
  	for ( i = 0; i < numBlcks; i ++ ) 
	{
    		beLoaded = GetLoad( i );	// get the flag
    		DIY_In_time_block ( 0, i, &tmGrpIdx );

    		if (  ( !beLoaded )  &&  ( tmGrpIdx == time_grp )  ) 
		{
			// DIY block info
      			DIY_Block_starts_sizes( 0,  i,  blockMin,   blockSiz );
      			spaceMin[0] = blockMin[0];    	spaceMax[0] = blockMin[0] + blockSiz[0] - 1;
			spaceMin[1] = blockMin[1];	spaceMax[1] = blockMin[1] + blockSiz[1] - 1;
			spaceMin[2] = blockMin[2];	spaceMax[2] = blockMin[2] + blockSiz[2] - 1;
      			minTmStp    = blockMin[3];	maxTmStp    = blockMin[3] + blockSiz[3] - 1;
 
			// real (non-ghost) block info
      			DIY_No_ghost_block_bounds( 0, i, &real_bb4 );
      			real_min[0] = real_bb4.min[0];	real_max[0] = real_bb4.max[0];
      			real_min[1] = real_bb4.min[1];	real_max[1] = real_bb4.max[1];
      			real_min[2] = real_bb4.min[2];	real_max[2] = real_bb4.max[2];
      			real_min[3] = real_bb4.min[3];	real_max[3] = real_bb4.max[3];
      			
			// attach the data blocks to the OSUFlow engine
      			start_tm = MPI_Wtime();
      			block_osuflow[i]->LoadData(  dataset_files, num_dataset_files, 
						     spaceMin,  spaceMax, 
						     real_min,  real_max,  volSizes,
						     minTmStp,  maxTmStp,  data_mode,  
						     dataBlks[i]  );
      			*timeUsed += ( MPI_Wtime() - start_tm );
		
			// update the flag
      			SetLoad( i );
      			dataSize += ( spaceMax[0] - spaceMin[0] ) * 
				    ( spaceMax[1] - spaceMin[1] ) * 
				    ( spaceMax[2] - spaceMin[2] ) * 3 * sizeof(float);
    		}
  	}

  	return  dataSize;
}


// -----------------------------------------------------------------------
//
// attach DIY curvilinear blocks (GRID blocks and DATA blocks, distributed
// to THIS process from process #0) to the OSUFlow engine
//
// the caller MUST destroy gridBlks, but must NOT destroy dataBlks
//
// time_grp:    index of the time group (or time block)
// timeUsed:    time used for this attachment process
// numBlcks:    number of spatial-temporal blocks of THIS process
//              (NOTE: in the time-varying case, NOT all these blocks are ACTIVE)
// sBlckIds[i]: ID of the SPATIAL block corresponding to spatial-temporal block #i
// gridBlks[i]: array of grid coordinates of the i-th (out of nsBlocks) SPATIAL 
//              block on THIS process
// dataBlks:    [numBlcks][timeGroupSize][blockSize]
//
// added by Zhanping Liu on 06/12/2013 and last updated on 06/19/2013
//
int 	Blocks::AttachCurvilinearGridDataBlocks4D
	( int   time_grp, double * timeUsed, int       numBlcks, 
	  int * sBlckIds, float ** gridBlks, float *** dataBlks )
{	
	cout << endl << "==== attaching DIY curvilinear blocks (time group #" 
	             << time_grp << ") to the OSUFlow engine ====" << endl;
	
	int  i, tmGrpIdx;
	int  blockMin[4], blockSiz[4];


	#ifdef BYTE_SWAP
  	for ( i = 0; i < numBlcks; i ++ )
	{ 
    		DIY_In_time_block( 0, i, &tmGrpIdx );
    		if ( tmGrpIdx != time_grp ) continue;

      		DIY_Block_starts_sizes( 0, i, blockMin, blockSiz );
		// solution data
      		for ( int j = 0; j < blockSiz[3]; j ++ ) 
		for ( int k = 0; k < blockSiz[0] * blockSiz[1] * blockSiz[2] * 3; k ++ )
	  	swap4(  ( char * ) &dataBlks[i][j][k]  );
  	}
	#endif


 	int 	dataSize = 0;
  	int 	beLoaded = 0;
  	double 	start_tm;

  	*timeUsed = 0.0;
  	for ( i = 0; i < numBlcks; i ++ ) 
	{
    		beLoaded = GetLoad( i );	// get the flag
    		DIY_In_time_block ( 0, i, &tmGrpIdx );

    		if (  ( !beLoaded )  &&  ( tmGrpIdx == time_grp )  ) 
		{
			// DIY block info
      			DIY_Block_starts_sizes( 0, i, blockMin, blockSiz );

			// attach the curvilinear blocks (GRID and DATA) to the OSUFlow engine
      			start_tm = MPI_Wtime();
			block_osuflow[i]->AttachCurvilinearGridData
              				  ( blockMin[3], blockSiz[3], blockSiz,    
		                            gridBlks[ sBlckIds[i] ],  dataBlks[i] );
      			*timeUsed += ( MPI_Wtime() - start_tm );
		
			// update the flag
      			SetLoad( i );
      			dataSize += blockSiz[0] * blockSiz[1] 
			          * blockSiz[2] * 3 * sizeof( float );
    		}
  	}

  	return  dataSize;
}


// ======================================================================= //
//                                                                         //
//                         above are the functions                         //
//                                                                         //
//                              in support of                              //
//                                                                         //
//                    integrating  VTK I/O with OSUFlow                    //
//                                                                         //
//             added by Zhanping Liu (05/23/2013 ~ 07/08/2013)             //
//                                                                         //
// ======================================================================= // ZPL end
#endif

//----------------------------------------------------------------------------
//
// DEPRECATED, removed by TP 10/10/12
//
// //
// // tests whether block b is in time block g
// //
// // g: current time block
// // blk: local block id
// // tsize: total number of timesteps
// // tb: total number of global time blocks
// //
// int Blocks::IsBlockInTimeGroup(int g, int b, int tsize, int tb) {

//   int64_t from[4]; // block starts
//   int64_t to[4]; // block starts
//   int min_t, max_t; // time range of the block

// #ifdef _MPI
//   int64_t time_group_start;
//   int64_t time_group_end;
//   blocking->GetRealBlockBounds(b, from, to);
//   min_t = from[3];
//   max_t = to[3];
//   blocking->GetRealTimeBounds(g, &time_group_start, &time_group_end);
// #else
//   int time_group_start;
//   int time_group_end;
//   if (lat4D) {
//     lat4D->GetRealTB(b, &min_t, &max_t);
//     lat4D->GetRealTimeGroupBounds(g, &time_group_start, &time_group_end);
//   }
//   else {
//     latAMR->GetTB(b, &min_t, &max_t);
//     time_group_start = g * tsize / tb;
//   }
// #endif

//   if (tsize == 1 || tb == 1 || (min_t == time_group_start && 
// 	                        max_t == time_group_end))
//     return 1;

//   return 0;

// }
//-----------------------------------------------------------------------
//
// evicts blocks in time group just prior to the current group
//
// grp: current group
// tsize: number of time steps
// tb: total number of global time blocks
// nblocks: number of local blocks
// ntblocks: total number of global time blocks
//
void Blocks::DeleteBlocks(int grp, int tsize, int tb, int nblocks) {

  int i;
  int loaded; // whehter a block is loaded

  if (grp == 0)
    return;

  for (i = 0; i < nblocks; i++) 
  {

#ifdef _MPI
    loaded = GetLoad(i);
#else
    loaded = lat4D ? lat4D->GetLoad(i) : latAMR->GetLoad(i);
#endif

    #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
    // edited by TP 10/10/20
    int time_block;
    DIY_In_time_block(0, i, &time_block);
//     if (loaded && IsBlockInTimeGroup(grp - 1, i, tsize, tb)) {
    if (loaded && time_block == grp - 1) 
    {
      // end TP
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

#ifdef _MPI
      ClearLoad(i);
#else
      lat4D ? lat4D->ClearLoad(i) : latAMR->ClearLoad(i);
#endif

      switch (compute_type) 
      {
        case OSUFLOW:
	  block_osuflow[i]->DeleteData();
	  break;
        default:
	  break;
      }
    #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
    }
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

  }

}

//-----------------------------------------------------------------------
//
// loads a new block
//
// grp: current group
// blk: desired block
// time (output): the time it took to load
// size: spatial domain size
// tsize: number of timesteps in dataset
// tb: total number of global time blocks
// nblocks: number of local blocks
// data: data that has already been read in (defaults to NULL)
//
// returns: number of bytes actually loaded
// 0 if the block is already in memory
// or its temporal extent does not match the group
//
int Blocks::LoadBlock4D(int grp, int blk, double *time, float *size,
			int tsize, int tb, int nblocks, float **data) {

  int s = 0; // data size (bytes)

  // edited TP 10/12/12
//   int64_t starts[4], sizes[4]; // block starts and sizes
  int starts[4], sizes[4]; // block starts and sizes

  int min_t, max_t; // block temporal extent
  int loaded; // whether a block is loaded

#ifdef _MPI
    loaded = GetLoad(blk);
#else
    loaded = lat4D ? lat4D->GetLoad(blk) : latAMR->GetLoad(blk);
#endif

    #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
    // edited by TP 10/10/20
    int time_block;
    DIY_In_time_block(0, blk, &time_block);
    //   if (!loaded && IsBlockInTimeGroup(grp, blk, tsize, tb)) {
    if (!loaded && time_block == grp) {
      // end TP
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

      // compute block extents
#ifdef _MPI // parallel version

      // edited TP 10/12/12
//       blocking->BlockStartsSizes(blk, starts, sizes);
      DIY_Block_starts_sizes(0, blk, starts, sizes);

      float from[3] = {starts[0], starts[1], starts[2]};
      float to[3] = {starts[0] + sizes[0] - 1, starts[1] + sizes[1] - 1,
		     starts[2] + sizes[2] - 1};
      int min_t = starts[3];
      int max_t = starts[3] + sizes[3] - 1;
      *time = MPI_Wtime();
#else // serial version
      float from[3], to[3]; // block spatial extent
      if (lat4D)
	lat4D->GetVB(blk, from, to, &min_t, &max_t);
      else 
	latAMR->GetVB(blk, from, to, &min_t, &max_t);
#endif

    // load the block
    switch (compute_type) {
    case OSUFLOW:
      block_osuflow[blk]->LoadData(dataset_files, num_dataset_files, from, to, 
				   size, min_t, max_t, data_mode, data); 
      break;
    default:
      break;
    }

#ifdef _MPI
    *time = (MPI_Wtime() - *time);
#else
      lat4D ? lat4D->SetLoad(blk) : latAMR->SetLoad(blk);
#endif

    s = ((to[0] - from[0]) * (to[1] - from[1]) *
	 (to[2] - from[2]) * 3 * sizeof(float));

  #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
  }
  #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012
  return s;

}
//-----------------------------------------------------------------------
//
// todo: diy AMR
//
// loads blocks in this time group
//
// grp: current group
// time (output): the time it took to load
// dm: data mode
//
// returns: number of bytes actually loaded
//
int Blocks::LoadBlocksAMR(int grp, double *time, DataMode dm) {

  int s = 0; // data size (bytes)
  float from[3], to[3]; // block spatial extent
  int min_t, max_t; // time group extent
  int i;
  float **data;
  int ntpart; // number of time partitions
  int tsize; // total number of timesteps
  int dims[3]; // block size (eg. 16x16x16)

  assert(dm == HDF_FLOAT || dm == HDF_DOUBLE);
  assert(latAMR != NULL);
  ntpart = latAMR->ntpart;
  tsize = latAMR->tdim;
  latAMR->GetBlockDims(dims);

  min_t = grp * tsize / ntpart;
  max_t = (grp == ntpart - 1 ? tsize - 1 : (grp + 1) * tsize / ntpart);

#ifdef _MPI
  *time = MPI_Wtime();
#endif

  for (i = min_t; i <= max_t; i++) // all time steps in the time group
    latAMR->LoadData(i, dm);

#ifdef _MPI
  *time  = MPI_Wtime() - *time;
#endif

  for (i = 0; i < latAMR->nb; i++ ) {

    #ifdef _MPI		// ADD-BY-LEETEN 10/29/2012
    // edited by TP 10/10/20
    // DIY does not support AMR right now, so this won't work
    int time_block;
    DIY_In_time_block(0, i, &time_block);
//     if (!latAMR->GetLoad(i) && IsBlockInTimeGroup(grp, i, tsize, ntpart)) {
    if (!latAMR->GetLoad(i) && time_block == grp) {
      // end TP
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012


      // create time varying flow field for this block
      latAMR->GetVB(i, from, to, &min_t, &max_t);
      data = latAMR->GetData(i);

      switch (compute_type) {
      case OSUFLOW:
	block_osuflow[i]->CreateTimeVaryingFlowField(data, dims[0], dims[1], 
						     dims[2], from, to, 
						     min_t, max_t); 
	break;
      default:
	break;
      }

      s += ((to[0] - from[0]) * (to[1] - from[1]) *
	    (to[2] - from[2]) * 3 * sizeof(float));
      latAMR->SetLoad(i);

    #ifdef _MPI			// ADD-BY-LEETEN 10/29/2012
    }
    #endif // #ifdef _MPI	// ADD-BY-LEETEN 10/29/2012

  }

  return s;

}
//-----------------------------------------------------------------------
//
// return the resident memory of the process in MB
//
// courtesy of Matt Sutter, UIUC
//
int mem_size(double *vsizeMB, double *residentMemoryMB) {


#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)

  // todo: work this out for mac
  *vsizeMB = (double) 0.0;
  *residentMemoryMB = (double) 0.0;

#endif

#ifdef LINUX

  // linux does not have getrusage, but has the /proc filesystem

  FILE *selfFile;
  
  // get resident use

  unsigned long int vsize;
  long int rss;

  int ret;
  char s[100];
  long pagesize = sysconf(_SC_PAGESIZE);

  selfFile = fopen("/proc/self/stat", "r");
  if (!selfFile) {
    perror("io_memory_usage could not open /proc/self/stat file (io_memory_usage only works on linux, right now)");
    *residentMemoryMB = -1.0;
    *vsizeMB = -1.0;
  } else {

    ret = fscanf(selfFile,"%*d %*s %*c %*d %*d %*d %*d %*d %*lu %*lu %*lu %*lu %*lu %*lu %*lu %*ld %*ld %*ld %*ld %*ld %*ld %*lu %lu %ld %*lu ", &vsize, &rss);
    if ((ret == EOF) || (ret < 2)) {
      fprintf(stderr, "[io_memory_usage] fscanf(/proc/self/stat) returned %d!\n", ret);
    }
    
    ret = fclose(selfFile);
    if (ret != 0) {
      fprintf(stderr, "[io_memory_usage] fclose(/proc/self/stat) returned %d!\n", ret);
    }

    rss = pagesize * rss;
    /* divide by 1024 * 1024 to convert to MB */
    *residentMemoryMB = ( (double) rss ) * 9.765625e-4 * 9.765625e-4;
    *vsizeMB = ( (double) vsize ) * 9.765625e-4 * 9.765625e-4;
  }

#endif

#ifdef BGP

  // getrusage is supposed to work on ibm

  struct rusage r_usage;

  getrusage ( RUSAGE_SELF, &r_usage);

  *residentMemoryMB = ( (double) r_usage.ru_maxrss ) * 9.765625e-4;
  *vsizeMB = (double) 0.0;

#endif

  return 0;

}
//-----------------------------------------------------------------------
