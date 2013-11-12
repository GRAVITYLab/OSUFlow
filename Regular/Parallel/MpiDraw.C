//------------------------------------------------------------------------------
//
// mpi test draw
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>
#include "OSUFlow.h"
#include "Draw.h"
#include "Blocks.h"
#include "ParFlow.h"
#include "diy.h"


// -------------------
//
// VTK I/O integration ---- added by Zhanping Liu on 05/23/2013 and last updated on 07/08/2013 ZPL begin
//


// NOTE: _USE_VTK_RENDERING_ is essentially a functionality added for _VTKIO_INTEGRATION_
//       while  it  is  exposed  here  (by NOT being under _VTKIO_INTEGRATION_  such that
//       the  original  code  (prior to _VTKIO_INTEGRATION_)  can  exploit  VTK rendering
//
#define  _USE_VTK_RENDERING_	      // to use the VTK renderer or the GCB renderer
#ifdef   _USE_VTK_RENDERING_          // search 'RenderFlowLines' below for appropriate tube radi
   //#define  _USE_SELF_ROTATION_     // to let  objects (flow lines herein) rotate automatically 
#endif			              // or allow interaction (through vtkRenderWindowInteractor)


#define  _VTKIO_INTEGRATION_	      // to integrate  VTK I/O with OSUFlow (ON) or not (OFF)


// VTK data formats
#ifdef   _VTKIO_INTEGRATION_
   #define _CURVILINEAR_GRIDS_        // load a curvilinear grid  dataset (either VTK or RAW)
   #define _LOAD_VTK_DATA_SET_	      // load VTK data (ON) or the counterpart RAW data (OFF)

   #define _DUMP_VTK_FILE_INF_	      // dump the info about a VTK data file
   //#define _GENERATE_RAW_DATA_      // when loading a VTK dataset,  create a RAW data file and its LIST file (for test only)
   //#define _TEST_VTKIO_RESULT_      // test if  the (serial data loading + block distribution) scheme  produces  exactly the
                                      // same result as that created by the  ORIGINAL BIL-based parallel data loading approach
   //#define _RECT_TO_CURV_DATA_      // save a rectilinear grid dataset as a curvilinear grid RAW dataset (for test purposes)
   //#define _CURVILINEAR_2_VTK_      // save an in-memory curvilinear grid dataset to a VTK file (of type  vtkStructuredGrid)
#endif			      


// VTK data access
#include "vtkCellData.h"	
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkStructuredGrid.h"	      // curvilinear grids
#include "vtkRectilinearGrid.h"
#include "vtkStructuredGridReader.h"  // curvilinear grids
#include "vtkStructuredGridWriter.h"  // curvilinear grids
#include "vtkRectilinearGridReader.h"
#include "vtkRectilinearGridWriter.h"


// VTK output and rendering
#include "vtkActor.h"
#include "vtkPoints.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTubeFilter.h"
#include "vtkRenderWindow.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataWriter.h"
#include "vtkRenderWindowInteractor.h"


// 
// ------------------- ZPL end


#ifdef MPE
#include "mpe_log.h"
#endif

// define this if you want to repartition between time groups
// #define REPARTITION

// WES - define this if you want to assign unique ids to each seed
// and then track these ids. the ids will be written to a separate file
// along with the normal output file
// #define TRACK_SEED_ID

// todo: this is a temporary hack
#ifdef REPARTITION
#define MAX_BLK 512 // maximum blocks per process
int wgts[MAX_BLOCK]; // weights of blocks
void AdvanceWeights(int g);
#endif

// drawing data at the root process
VECTOR4 *pt; // points in everyone's traces
int *npt; // everyone's number of points in their traces
int tot_ntrace; // total number of everyone's traces

// performance stats
double TotTime = 0.0; // total time
double TotInTime = 0.0; // total input time
double TotOutTime = 0.0; // total output time
double TotCompCommTime = 0.0; // comp + comm time
int TotParticles; // total number of particles in the system


// -------------------
//
// VTK I/O integration ---- added by Zhanping Liu on 05/23/2013 and last updated on 07/08/2013 ZPL begin
//


int           maxBlcks = 1;	// max number of spatial-temporal blocks on a single process
int           max_tblk = 1;	// max number of time steps per 'time block' (or time group)
int	      nsBlocks = 1;	// number of SPATIAL blocks on THIS / EACH process
int	    * sBlckIds = NULL;	// array: SPATIAL index of each of 'nblocks' spatial-temporal blocks on THIS process
int         * stBlkIds = NULL;  // array: index of the FIRST spatial-temporal block corresponding to each of 'nsBlocks' SPATIAL blocks on THIS process
int	    * STBLKIDS = NULL;	// array: index of the FIRST spatial-temporal block corresponding to each SPATIAL block on ALL processes (maintained by process #0)
int	    * tBlkMins = NULL;  // array: index  of the beginning time step  of each time group
int         * tBlkSizs = NULL;  // array: number of the           time steps of each time group
int         * tmGrpIds = NULL;  // array: time group index of each spatial-temporal block on a single process
int	    * beActive = NULL;  // array: status (active / inactive) of each spatial-temporal block on a single process
int         * BEACTIVE = NULL;  // array: status (active / inactive) of ALL  spatial-temporal blocks on ALL processes (maintained by process #0)
int         * bb_infor = NULL;  // array: ALL (ghost) bounding-box info on THIS process
int 	    * BB_INFOR = NULL;	// array: ALL (ghost) bounding-box info of ALL processes (collected and maintained by process #0)
int         * numBlcks = NULL;	// array: number of (SPATIAL-TEMPORAL) blocks on each process
float	    * grd_cord = NULL;	// array: the grid data (3D coordinates) of the whole volume (loaded and maintained by process #0 only)
float      ** gridBlks = NULL;  // array: gridBlks[i] refers to the coordinates of the i-th (out of nsBlocks) SPATIAL block on THIS process
float	   ** vec_data = NULL;	// array: vec_data[timeGroupSize][VolumeSize] for loading and maintaining multiple time steps of data
float     *** dataBlks = NULL;	// array: dataBlks[numBlocks][timeGroupSize][subVolumeSize] for blocks (i.e., sub-volumes) maintained by EACH process
vtkPolyData * polyData = NULL;  // a vtkPolyData object for VTK output and rendering of flow lines


//
// ------------------- ZPL end


// globals
static char filename[256]; // dataset file name
static char **dataset_files = NULL; // files in the dataset
int num_dataset_files = 0;
char part_file[256]; // partition file name
float size[3]; // spatial domain size
static int tsize; // temporal domain size
vector< vector<Particle> > Seeds; // list of seeds lists
OSUFlow **osuflow; // one flow object for each block
list<vtListTimeSeedTrace*> *sl_list; // fieldlines list
int nspart; // global total number of spatial blocks
int ntpart; // global total number of temporal blocks
int nblocks; // my number of blocks
int tf; // max number of traces per block
int pf = 1000; // max number of points per trace in each round
const int max_rounds = 100; // max number of rounds

// deleted TP 10/23/12
// const int check_rounds = 5; // how often to flush / check seeds

int end_steps; // final number of points each particle should travel
DataMode data_mode; // data format
Blocks *blocks; // block class object
ParFlow *parflow; // parallel flow class object

// deleted TP 10/12/12
// Blocking *blocking; // blocking class object
// Assignment *assign; // assignment class object
// Neighborhoods *nbhds; // neighborhoods class object
// end TP

int compute_begin, compute_end; // jumpshot states
int rank; // mpi rank
float vec_scale; // vector scaling factor
char seed_file[256]; // seed file name
int seed_file_num = 0; // number of seeds read from the seed file
VECTOR3* seed_file_seeds = NULL; // seeds read from the seed file
float wf = 0.1; // wait factor for nonblocking communication
                // wait for this portion of messages to arrive each round

// deleted TP 10/23/12
// const int ghost = 1;  // a ghost layer of at least 1 is required for correct
// 		      // advection between blocks

// integration parameters
const float maxError = 0.001;
const float initialStepSize = 1.0;
const float minStepSize = 0.01;
const float maxStepSize = 5.0;
const float lowerAngleAccuracy = 3.0;
const float upperAngleAccuracy = 15.0;

const INTEG_ORD integrationOrder = RK45;
//const INTEG_ORD integrationOrder = FOURTH;
const bool useAdaptiveStepSize = true;

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void Cleanup();
void Run(MPI_Comm comm);
void Header(char *filename, float *size, int *tsize, float *vec_scale);
void LoadSeedsFromFile();
int isSeedInTimeGroup(int g);
bool isSeedInTimeGroupTotal(int g);
int getNumSeedsInTimeGroup(int g);


// -------------------
//
// VTK I/O integration ---- added by Zhanping Liu on 05/23/2013 and last updated on 07/08/2013 ZPL begin
//


// VTK data generation
void	VEC2VTK( char * vec_name, char * vtk_name );                                  // convert a  single  RAW  VEC data  file
void	TornadoDataFiles2VTK( char * thePrefx, int fileIdx0 = 1, int fileIdx1 = 50 ); // convert time-varying RAW tornado files
void	Rectilinear2CurvilinearRaw      // curvilinear grid data but in RAW binary format
	( int   * gridSizs, float * gridVecs, char * pathName, char * mainName );
void	CurvilinearGridData2VTK		// in-memory curvilinear data (grid + vector) to a VTK file (of type vtkStructuredGrid)
	( int   * dimensns, float * curv_grd, 
	  float * curv_vec, char  * fileName, int    fileType = VTK_BINARY );


// serial data loading
void 	LoadCartesianRawData  ( char * raw_file, float *   vec_buff ); // for only one time step
int  	LoadRectilinearVTKdata( char * rectFile, float *   vec_buff ); // for only one time step, used as CARTESIAN
void  	LoadCurvilinearRawGrid( char * gridFile, float * & grd_buff ); // once for ALL time steps: only GRID (without VECTOR DATA)
void  	LoadCurvilinearRawData( char * vec_file, float *   vec_buff ); // one for each time step:  only VECTOR DATA (without GRID)
void	LoadCurvilinearVtkGrid( char * curvFile, float * & grd_buff ); // once for ALL time steps: only GRID (without VECTOR DATA)
int	LoadCurvilinearVtkData( char * curvFile, float *   vec_buff ); // one for each time step:  only VECTOR DATA (without GRID)
void	LoadCurvilinearVtkData_DME	// with  INTERMEDIATE attributes (Density, 3 Menmentum components, and Energy)
	                      ( char * curvFile, float *   vec_buff ); // from  which  3D velocity vectors are derived


// blocks distribution
void 	FillBlock( int * volSizes, float * vol_data, // use the whole volume (vector data  or  grid coordinates,  the latter for
		   int * blckMins, int   * blkSizes, float * blckData ); // NON-Carteisan grids) to fill in a block (sub-volume)
void 	InitBlocksDistributor(); // init the context for  the distribution of blocks (vector data and / or grid coordinates)
void 	DistributeGridBlocks();	 // distribute DIY GRID blocks (coordinates)                       (for NON-Cartesian grids)
void	AttachPhysicalBounds();  // attach (ghost and real) physical bounds to each DIY GRID block (for NON-Cartesian grids)
void 	DistributeDataBlocks( int tmGrpIdx );  // distribute DIY data blocks (vectors): this is common to ALL kinds of grids


// placement of seeds for each block
int	InitTraces4CurvilinearBlocks( char * seedFile );  // import seeds from a file (for fixing a bug with CurvilinearGrid.C)
void	InitTraces4CurvilinearBlocks();  // since  a random seed placement scheme might not work (well)  for  curvilinear grids


// VTK output and rendering
void	AddFlowLine2Polydata
	( int         numSamps, VECTOR4      * strmSmps,
          vtkPoints * vtk_pnts, vtkCellArray * vtkLines,
	  int         procIndx, vtkIntArray  * procIdxs,
          int         strmIndx, vtkIntArray  * strmIdxs, vtkIntArray * strmLens );
void	FlowLines2PolyData
	( vtkPolyData * polyData,     int   pIdxInfo,
	  int           sIdxInfo,     int   sLenInfo,
	  int           numProcs = 0, int * numStrms = NULL, int actvScal = 0 );
void	FlowLines2PolyData
	( int pIdxInfo,     int    sIdxInfo,        int sLenInfo,
	  int actvScal = 0, char * vtk_file = NULL, int fileType = VTK_BINARY );
void	RenderFlowLines( int be_tubed = 1, float tube_rad = 0.1, int tube_res = 20 );


// memory clean-up
void 	ReleaseMemory();


//
// ------------------- ZPL end


//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  int nproc; // mpi groupsize

  // init
  MPI_Init(&argc, &argv);
#ifdef MPE
  MPE_Log_get_state_eventIDs(&compute_begin, &compute_end);
  MPE_Describe_state(compute_begin, compute_end, "Compute", "red");
#endif


  GetArgs(argc, argv);
  

  // deleted TP 10/12/12
// #ifdef USE_BIL
//   BIL_Init(MPI_COMM_WORLD);
// #endif

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  Init();

  // run
  MPI_Barrier(MPI_COMM_WORLD);
  TotTime = MPI_Wtime();
  Run(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  TotTime = MPI_Wtime() - TotTime;

  // print the performance stats
  parflow->PrintPerf(TotTime, TotInTime, TotOutTime,
		     TotCompCommTime, TotParticles, size);

#ifdef GRAPHICS
  
  
  #ifdef _USE_VTK_RENDERING_	// added by Zhanping Liu on 05/28/2013                        ZPL begin
  FlowLines2PolyData( 1, 1, 1, 0, "flowlines.vtk" );
  RenderFlowLines( 1, 0.02 ); 	// 1 / 0: use tubes or not; tube radius: 0.1 for vtkDataRaw & RectGrid2 & curvGridRaw & curvGridVTK, 
  #else				// 0.35 for test & Tornado, and 0.02 for office & Blunt Fin   ZPL end
  if (rank == 0)
  {
    VECTOR3 min, max;
    min = VECTOR3(0.0f, 0.0f, 0.0f);
    max = VECTOR3((float)(size[0] - 1), (float)(size[1] - 1),
		  (float)(size[2] - 1));
    DrawInit(pt, npt, tot_ntrace, argc, argv, min, max, 0 );
  }
  #endif                        // ZPL


#endif

  Cleanup();
  MPI_Barrier(MPI_COMM_WORLD);

  // edited TP 10/12/12
// #ifdef USE_BIL
//   BIL_Finalize();
// #endif
  DIY_Finalize();
  // end TP

  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// Run
//
void Run(MPI_Comm comm) {

  int i, j;
  int g; // current group
  int g_io; // the last group in which blocks were loaded
  double time; // time to load a block
  double t0;

#ifdef REPARTITION
  assert(nblocks <= MAX_BLK);
  for (i = 0; i < nblocks; i++)
    wgts[i] = 0;
#endif


  #ifdef _CURVILINEAR_GRIDS_  // added by Zhanping Liu on 06/14/2013 ZPL begin
    InitTraces4CurvilinearBlocks();
    //InitTraces4CurvilinearBlocks( "/home/zhanping/Integration/Data/curvSeeds_RectGrid2.dat" );
  #else                       // ZPL end
    parflow->InitTraces(Seeds, tf, nblocks, tsize, ntpart, seed_file_seeds, seed_file_num);
  #endif                      // ZPL


#ifdef USE_BIL
  float ***bil_data = NULL;
#endif
  
  // for all groups
  for (g = 0; g < ntpart; g++) 
  {	
    	// check if there is any work to be done for this time group
    	if( !isSeedInTimeGroupTotal(g) ) continue;  // go to next time group
	
    	// synchronize before starting I/O
    	MPI_Barrier(comm);
    	t0 = MPI_Wtime();

    	// delete blocks from previous time group
    	if (g > 0) 
    	{ 
      		blocks->DeleteBlocks(g_io+1, tsize, ntpart, nblocks);
    	}
	
    	// todo: change seeds to vector in repartition
    	// #ifdef REPARTITION
    	//     parflow->Repartition(g, &nblocks, &Seeds, 0, &osuflow,
    	// 		       block, OSUFLOW, MPI_COMM_WORLD, wgts);
    	// #endif


        // data blocks distribution (changes made to Blocks.h/C  accordingly) ZPL begin
	//
	// added by Zhanping Liu on 05/26/2013 and last updated on 07/08/2013
	//
    	#ifdef _VTKIO_INTEGRATION_
     

  	int     thisProc; 													
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );

	if ( thisProc == 0 )
	cout << "==== loading data for time group #"
	     << g << " (time steps #" << tBlkMins[g] << " ~ #" 
	     << tBlkMins[g] + tBlkSizs[g] - 1 << ") ====" << endl << endl;


	for ( i = 0; i < tBlkSizs[g]; i ++ )
	{
		#ifdef _LOAD_VTK_DATA_SET_
			#ifdef  _CURVILINEAR_GRIDS_
				LoadCurvilinearVtkData( dataset_files[ tBlkMins[g] + i ], vec_data[i] ); 	  // VTK curvilinear
						
				// the line below is used to load the ORIGINAL blunt-fin dataset (curvGrid_bluntFin_DME.vtk coupled
				// with curvGrid_bluntFin_DME.list) from which the velocity vectors have to be derived from the DME
				// (Density, Momentum, and Energy) attribtes
				//
				// since   this  original  dataset   has  been  converted  to  curvGrid_bluntFin.vtk  (coupled  with 
				// curvGrid_bluntFin.list) where the velocity vectors are readily available,  the line above  should
				// be used instead
				//
				// the line below remains here to maintain the capability of reading (other) DME curvilinear data
				//
				// LoadCurvilinearVtkData_DME( dataset_files[ tBlkMins[g] + i ], vec_data[i] );   // VTK curvilinear
			#else	
				// NOTE: though treated / used as a CARTESIAN grid		
				LoadRectilinearVTKdata( dataset_files[ tBlkMins[g] + i ], vec_data[i] ); 	  // VTK rectilinear
			#endif
  		#else
			#ifdef  _CURVILINEAR_GRIDS_
				LoadCurvilinearRawData( dataset_files[ tBlkMins[g] + i ], vec_data[i] ); 	  // raw curvilinear
			#else
				LoadCartesianRawData  ( dataset_files[ tBlkMins[g] + i ], vec_data[i] ); 	  // raw Cartesian
			#endif
  		#endif
	}
			
  	DistributeDataBlocks( g );  // dataBlks[][][] is ready AFTER this call --- distribution of data blocks (vectors)

        // attach the data blocks to the underlying OSUFlow engine
	#ifdef _CURVILINEAR_GRIDS_
	blocks->AttachCurvilinearGridDataBlocks4D( g, &time, nblocks, sBlckIds, gridBlks, dataBlks );
	#else
    	blocks->AttachDataBlocks4D( g, &time, nblocks, size, dataBlks );
	#endif


    	#else	// _VTKIO_INTEGRATION_ is OFF ZPL end


    		// load blocks for this time group
		#ifdef USE_BIL
    			bil_data = blocks->BilLoadTimeGroupBlocks(g, nblocks, size, tsize, ntpart);
    			blocks->LoadBlocks4D(g, &time, nblocks, size, tsize, ntpart, bil_data);
    			delete[] bil_data;
		#else
    			blocks->LoadBlocks4D(g, &time, nblocks, size, tsize, ntpart);
		#endif


        #endif  // _VTKIO_INTEGRATION_        ZPL


    	g_io = g;

    	// synchronize after I/O
    	MPI_Barrier(comm);
    	TotInTime += (MPI_Wtime() - t0);
    	t0 = MPI_Wtime();

    	// scale blocks to improve visibility
    	for (i = 0; i < nblocks; i++) 
	{
      		if (blocks->GetLoad(i))
		osuflow[i]->ScaleField(vec_scale);
    	}

	#ifdef REPARTITION
    	assert(nblocks <= MAX_BLK);
    	for (i = 0; i < nblocks; i++) wgts[i] = 0;
	#endif

    	// for all rounds
    	for (j = 0; j < max_rounds; j++) 
	{
      		// for all blocks
      		for (i = 0; i < nblocks; i++) 
		{
			#ifdef MPE
			MPE_Log_event(compute_begin, 0, NULL);
			#endif

			// compute fieldlines
			if (tsize > 1) 
			{
				#ifdef REPARTITION
	  			assert(i < MAX_BLK);
	  			parflow->ComputePathlines(Seeds[i], i, pf, end_steps, &wgts[i]);
				#else
	  			parflow->ComputePathlines(Seeds[i], i, pf, end_steps);
				#endif
			}
			else 
			{
				#ifdef REPARTITION
	  			assert(i < MAX_BLK);
	  			parflow->ComputeStreamlines(Seeds[i], i, pf, end_steps, &wgts[i]);
				#else
	  			parflow->ComputeStreamlines(Seeds[i], i, pf, end_steps);
				#endif
			}

			#ifdef MPE
			MPE_Log_event(compute_end, 0, NULL);
			#endif

      		} // for all blocks

      		// exchange neighbors
      		parflow->ExchangeNeighbors(Seeds, wf);

      		// deleted TP 10/12/12
		//       if(j % check_rounds == 0)
		//       {
		// 	// check if there is any more work to do in this time group
		// 	parflow->FlushNeighbors(Seeds);

		// 	if(!isSeedInTimeGroupTotal(g))
		// 	{
		// 	  break;  // break out of loop going through every round
		// 	}
		//       }
		// end TP

    	} // for all rounds

    	// flush any remaining messages
    	parflow->FlushNeighbors(Seeds);

	#ifdef REPARTITION
    	AdvanceWeights(g);
	#endif

    	// end time group synchronized to get accurate timing
    	MPI_Barrier(comm);
    	TotCompCommTime += (MPI_Wtime() - t0);

  } // for all groups

  // synchronize prior to gathering
  MPI_Barrier(comm);
  TotOutTime = MPI_Wtime();

  // gather fieldlines for rendering
  parflow->GatherFieldlines(nblocks, size, tsize);
  MPI_Barrier(comm);
  TotOutTime = MPI_Wtime() - TotOutTime;

  //   // debug
  //   int rank;
  //   double vsize, rsize; // virtual and resident memory size in MB
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //   mem_size(&vsize, &rsize);
  //   fprintf(stderr, "rank = %d vsize = %.3lf MB, rsize = %.3lf MB\n", 
  // 	  rank, vsize, rsize);

}
//-----------------------------------------------------------------------

#ifdef REPARTITION

//-----------------------------------------------------------------------
//
// advances the weights of spatial blocks in the current time block to the 
// to the same spatial block in the next time block
//
// g: current time block
//
void AdvanceWeights(int g) {

  int nwts; // my number of weights (blocks) in this time block
  int tot_nwts; // total number of weights from all processes
  static int *recv_counts; // number of weights from each proess
  static int *recv_displs; // displacements for weights from each process
  int *all_wts; // all the weights from all processes
  int wts[2 * MAX_BLK]; // my block ranks and weights to send to everyone
  static int groupsize;
  static int first = 1;
  int i, j, n;

  if (first) {
    MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
    assert((recv_counts = new int[groupsize]) != NULL);
    assert((recv_displs = new int[groupsize]) != NULL);
    first = 0;
  }

  // count the number of blocks in my time group and exchange with everyone
  assert(nblocks <= MAX_BLK);
  n= 0;
  nwts = 0; // number of my weights (blocks) in this time block
  for (i = 0; i < nblocks; i++) {
    if (block->IsBlockInTimeGroup(g, i, tsize, ntpart)) {
      wts[n++] = nbhds->Lid2Gid(i);
      wts[n++] = wgts[i];
      nwts++;
    }
  }
  MPI_Allgather(&nwts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

  // allocate space for receiving all the weights
  tot_nwts = 0;
  for (i = 0; i < groupsize; i++)
    tot_nwts += recv_counts[i];
  assert((all_wts = new int[2 * tot_nwts]) != NULL);

  // exchange weights
  recv_displs[0] = 0;
  recv_counts[0] *= 2;
  for (i = 1; i < groupsize; i++) {
    recv_counts[i] *= 2;
    recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
  }
  MPI_Allgatherv(wts, 2 * nwts, MPI_INT, 
		all_wts, recv_counts, recv_displs, MPI_INT, MPI_COMM_WORLD);

  // for all my blocks in the next time group, find the weight in the
  // current time group from all_wts
  for (i = 0; i < nblocks; i++) {
    if (block->IsBlockInTimeGroup4D(g + 1, i)) {
      for (j = 0; j < tot_nwts; j++) {
	if (all_wts[2 * j] == nbhds->Lid2Gid(i) - nspart)
	  break;
      }
      if (j < tot_nwts) // sanity check, but always true
	wgts[i] = all_wts[2 * j + 1]; // advance the weight
    }
  }

  delete[] all_wts;

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// GetArgs
//
// gets command line args
//
void GetArgs(int argc, char *argv[]) {

  int groupsize;

  MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
  assert(argc >= 9);
		
  strncpy( filename, argv[1], sizeof(filename) );
  Header ( filename, size, &tsize, &vec_scale );
  
  nspart = groupsize * atoi(argv[2]); // total space partitions
  ntpart = (tsize == 1 ? 1 : atoi(argv[3])); // total time partitions
  tf = atoi(argv[4]) / nspart; // traces per block
  end_steps = atoi(argv[5]); // desired ending number of steps per field line
  strncpy(part_file, argv[6], sizeof(part_file));
  switch(atoi(argv[7])) {
  case 0:
    data_mode = RAW;
    break;
  case 1:
    data_mode = RAW_HEADER;
    break;
  case 2:
    data_mode = NETCDF;
    break;
  case 3:
    data_mode = HDF_FLOAT;
    break;
  case 4:
    data_mode = HDF_DOUBLE;
    break;
  default:
    break;
  }
  strncpy(seed_file, argv[8], sizeof(seed_file));

  pf = end_steps;

}
//-----------------------------------------------------------------------
//
// SetArgs
//
// sets parameters from given values
// SWIFT_TEST: args = $data $bf $tb $tp $st $pf $dm $sf
//
void SetArgs(MPI_Comm comm, const char *data,
             int bf, int tb, int tp, int st, const char* pfile,
             DataMode dm, const char *sf) {

  int groupsize;

  MPI_Comm_size(comm, &groupsize);
  strncpy(filename, data, sizeof(filename));
  Header(filename, size, &tsize, &vec_scale);
  
  nspart = groupsize * bf; // total space partitions
  ntpart = (tsize == 1 ? 1 : tb); // total time partitions
  tf = tp / nspart; // traces per block
  end_steps = st; // desired ending number of steps per field line
  strncpy(part_file, pfile, sizeof(part_file));
  switch(dm) {
  case 0:
    data_mode = RAW;
    break;
  case 1:
    data_mode = RAW_HEADER;
    break;
  case 2:
    data_mode = NETCDF;
    break;
  case 3:
    data_mode = HDF_FLOAT;
    break;
  case 4:
    data_mode = HDF_DOUBLE;
    break;
  default:
    break;
  }
  strncpy(seed_file, sf, sizeof(seed_file));

  pf = end_steps;

}


//-----------------------------------------------------------------------
//
// Init
//
// inits the app
//
void Init() {

  bool track_seed_ids = false;
#ifdef TRACK_SEED_ID
  track_seed_ids = true;
#endif

  int myproc, nproc; // usual MPI
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);


  // this code segment is used to convert the time-varying tornado ZPL begin
  // data set  (with 50 time steps)  from  VEC  files to VTK files
  //
  // added by Zhanping Liu on 06/02/2013
  //
  // if ( myproc == 0 ) 
  // TornadoDataFiles2VTK( "/home/zhanping/Tornado/", 1, 50 );  // ZPL end


  assert(nspart * ntpart >= nproc);
  assert(ntpart <= tsize);
  if(seed_file[0] == '!')
    assert(tf > 0);

  // partition domain
  // todo: don't create partition if there is a part file?
  int data_size[4] = {size[0], size[1], size[2], tsize};
  int given[4] = {0, 0, 0, ntpart}; // constraints in x, y, z, t
  int ghost[8] = {1, 1, 1, 1, 1, 1, 0, 0}; // -x, +x, -y, +y, -z, +z, -t, +t
  if (tsize > 1) 
    ghost[7] = 1;
  DIY_Init(4, data_size, 1, MPI_COMM_WORLD);
  DIY_Decompose(ROUND_ROBIN_ORDER, nspart * ntpart, &nblocks, 1, ghost, given);


  #ifdef _VTKIO_INTEGRATION_	    // added by Zhanping Liu ZPL begin

  InitBlocksDistributor();          // added on 06/02/2013

  	#ifdef _CURVILINEAR_GRIDS_
  	DistributeGridBlocks();	    // added on 06/19/2013
	AttachPhysicalBounds();     // added on 07/02/2013
  	#endif

  #endif			    // ZPL end
	

  // create osuflow object for each block
  // todo: switch to vectors and get rid of memory management
  assert((osuflow = (OSUFlow**)malloc(nblocks * sizeof(OSUFlow))) != NULL);
  for (i = 0; i < nblocks; i++)
    osuflow[i] = new OSUFlow;

  // Seeds and fieldline list
  Seeds.resize(nblocks);
  sl_list = new list<vtListTimeSeedTrace*>[nblocks];

  // create remaining classes
  // todo: rename

  // edited TP 10/12/12
//   blocks = new Blocks(blocking, assign, (void *)osuflow, OSUFLOW, 
// 		      dataset_files, num_dataset_files, data_mode, ghost);
//   parflow = new ParFlow(blocking, assign, blocks, osuflow, sl_list, 
// 			&pt, &npt, &tot_ntrace, nblocks, 0);
  blocks = new Blocks(nblocks, (void *)osuflow, OSUFLOW, 
		      dataset_files, num_dataset_files, data_mode);

  parflow = new ParFlow(blocks, osuflow, sl_list, 
			&pt, &npt, &tot_ntrace, nblocks, 0);
  // end TP

  
  // this code segment supports the (current) curvilinear grid case  ZPL begin
  // without affecting the  original code (NON-curvilinear grids) &
  // without any performance penalty AT ALL
  //
  // added by Zhanping Liu on 07/11/2013
  //
  #ifdef _CURVILINEAR_GRIDS_
  parflow->SetGrid2Curvilinear();
  #endif                                                          // ZPL end


  parflow->SetMaxError(maxError);
  parflow->SetInitialStepSize(initialStepSize);
  parflow->SetMinStepSize(minStepSize);
  parflow->SetMaxStepSize(maxStepSize);
  parflow->SetLowerAngleAccuracy(lowerAngleAccuracy);
  parflow->SetUpperAngleAccuracy(upperAngleAccuracy);
  parflow->SetIntegrationOrder(integrationOrder);
  parflow->SetUseAdaptiveStepSize(useAdaptiveStepSize);

  TotParticles = nspart * tf;

  if(seed_file[0] != '!') {
    LoadSeedsFromFile();
    TotParticles = seed_file_num;
    assert(seed_file_num > 0);
  }

}
//-----------------------------------------------------------------------
//
// Cleanup
//
// frees memory and such
//
void Cleanup() {

  int i;

  delete [] pt;
  delete [] npt;

  
  // ------------------- ZPL begin
  //

  #ifdef _VTKIO_INTEGRATION_
  ReleaseMemory();	// added by Zhanping Liu on 05/28/2013
  #endif


  // NOTE: _USE_VTK_RENDERING_ is essentially a functionality  added for 
  //       _VTKIO_INTEGRATION_, whereas it is exposed here (by NOT being
  //       under _VTKIO_INTEGRATION_) to the original code for rendering
  //
  // added by Zhanping Liu on 07/08/2013
  // 
  #ifdef _USE_VTK_RENDERING_
  if ( polyData ) polyData->Delete();	polyData = NULL;  // for process #0 only
  #endif

  //
  // ------------------- ZPL end

  
  for(i=0; i<nblocks; i++)
  {
    list<vtListTimeSeedTrace*>::iterator trace_iter;
    for(trace_iter=sl_list[i].begin();trace_iter!=sl_list[i].end();trace_iter++)
    {
      vtListTimeSeedTrace::iterator pt_iter;
      for(pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end(); 
	  pt_iter++)
      {
	delete *pt_iter;
      }
      (*trace_iter)->clear();
      delete *trace_iter;
    }
    sl_list[i].clear();
  }
  delete [] sl_list;

  for (i = 0; i < nblocks; i++)
  {
    if (osuflow[i] != NULL)
    {
      delete osuflow[i];
    }
  }
  for (i = 0; i < Seeds.size(); i++)
    Seeds[i].clear();
  Seeds.clear();

  delete blocks;
  delete parflow;

  // deleted TP 10/12/12
//   delete blocking;

  free(osuflow);

  for (i = 0; i < num_dataset_files; i++)
    free(dataset_files[i]);
  free(dataset_files);

}
//-----------------------------------------------------------------------
//
// reads header file
//
// filename: filename (input)
// size: data spatial size (output)
// tsize: number of timesteps (output)
// vec_scale: vector scaling factor (output)
//
void Header(char *filename, float *size, int *tsize, float *vec_scale) {
  
  FILE *fp;
  char line[1024];
  char *token;
  const char delims[] = " \t\n";
  int n = 0; // number of tokens parsed so far

  fp = fopen(filename, "r");
  assert(fp != NULL);

	// ADD-BY-LEETEN 11/19/2011-BEGIN
  char szPath[1024];
  strcpy(szPath, filename);

  char *szSeparator;
  szSeparator = strrchr(szPath, '/');
        #ifdef WIN32
  if( !szSeparator )
    szSeparator = strrchr(szPath, '\\');
        #endif
  if( NULL == szSeparator )
    szSeparator = &szPath[0];

  *szSeparator = '\0';
	// ADD-BY-LEETEN 11/19/2011-END

  while (fgets(line, sizeof(line), fp) != NULL) {

    token = strtok(line, delims);

    while (*line != '\0' && *line != '\n' && *token != '#') {

      if (n < 3)
	size[n++] = atof(token);

      else if (n == 3) {
	*tsize = atoi(token);
	num_dataset_files = *tsize;
	assert((dataset_files = (char **)malloc(sizeof(char *) * 
						num_dataset_files)) != NULL);
	n++;
      }

      else if (n == 4) {
	*vec_scale = atof(token);
	n++;
      }

      else if (n >= 5) {
	if (n - 5 >= *tsize)
	  break;
	// MOD-BY-LEETEN 11/19/2011-FROM:
		// dataset_files[n - 5] = strdup(token);
	// TO:
	char *szFilename = strdup(token);
	dataset_files[n - 5] = (char*)calloc(strlen(szPath) + strlen(szFilename) + 128, 1); // 128: size of a temp buffer
	if( szFilename[0] != '/' )
	  sprintf(dataset_files[n - 5], "%s/%s", szPath, szFilename);
	else
	  strcpy(dataset_files[n - 5], szFilename);
	free(szFilename);
	// MOD-BY-LEETEN 11/19/2011-END
	assert(dataset_files[n - 5] != NULL);
	n++;
      }

      token = strtok(NULL, delims);
      if (token == NULL)
	break;
    }

  }

  fclose(fp);

}
//-----------------------------------------------------------------------
//
// LoadSeedsFromFile
//
// loads seed locations from a file
//
void LoadSeedsFromFile()
{
  FILE* fileptr = fopen(seed_file, "r");
  assert(fileptr != NULL);

  // first line is the number of seeds
  fscanf(fileptr, "%i", &seed_file_num);

  // next lines are x y z seed positions
  float x, y, z;
  seed_file_seeds = new VECTOR3[seed_file_num];
  for(int i=0; i<seed_file_num; i++)
  {
    fscanf(fileptr, "%f %f %f", &x, &y, &z);
    seed_file_seeds[i].Set(x, y, z);
  }

  fclose(fileptr);
}
//-----------------------------------------------------------------------
//
// returns true if any proc has any seed still in time group g
//
bool isSeedInTimeGroupTotal(int g)
{
  int hasSeeds = isSeedInTimeGroup(g);
  int hasSeedsAll;
  MPI_Allreduce(&hasSeeds, &hasSeedsAll, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

  return hasSeedsAll;
}
//-----------------------------------------------------------------------
//
// returns 1 if there is at least one seed in the time group, only looks
// at seeds local to this process.
//
int isSeedInTimeGroup(int g)
{
  if(tsize == 1)
  {
    // steady flow case
    vector< vector<Particle> >::iterator sl_iter;  // seed list iterator
    for(sl_iter=Seeds.begin(); sl_iter!=Seeds.end(); sl_iter++)
    {
      if(sl_iter->size() > 0)
      {
	return 1;
      }
    }
  }
  else
  {
    // time-varying case

    // edited TP
    int64_t tmin, tmax;
//     blocking->GetRealTimeBounds(g, &tmin, &tmax);
    int time_block;
    bb_t bb;
    for (int i = 0; i < nblocks; i++) {
      DIY_In_time_block(0, i, &time_block);
      if (time_block == g) {
	DIY_No_ghost_block_bounds(0, i, &bb);
	tmin = bb.min[3];
	tmax = bb.max[3];
	break;
      }
    }
    // end TP

    vector< vector<Particle> >::iterator sl_iter;  // seed list iterator
    for(sl_iter = Seeds.begin(); sl_iter != Seeds.end(); sl_iter++)
    {
      vector<Particle>::iterator seed_iter;
      for(seed_iter=sl_iter->begin(); seed_iter!=sl_iter->end(); seed_iter++)
      {
	float time = seed_iter->pt[3];
	if(time >= tmin && time <= tmax)
	{
	  return 1;
	}
      }
    }
  }

  return 0;
}
//-----------------------------------------------------------------------
//
// returns the number of seeds in the time group
//
// counts seeds that are currently within this process, and is in the time
// group g. this function is mainly used for debugging.
//
int getNumSeedsInTimeGroup(int g)
{
  int count = 0;
  if(tsize == 1)
  {
    // steady flow case
    vector< vector<Particle> >::iterator sl_iter;  // seed list iterator
    for(sl_iter=Seeds.begin(); sl_iter!=Seeds.end(); sl_iter++)
    {
      count += sl_iter->size();
    }
  }
  else
  {
    // time-varying case
    int64_t tmin, tmax;

    // edited TP
//     blocking->GetRealTimeBounds(g, &tmin, &tmax);
    int time_block;
    bb_t bb;
    for (int i = 0; i < nblocks; i++) {
      DIY_In_time_block(0, i, &time_block);
      if (time_block == g) {
	DIY_No_ghost_block_bounds(0, i, &bb);
	tmin = bb.min[3];
	tmax = bb.max[3];
	break;
      }
    }
    // end TP

    vector< vector<Particle> >::iterator sl_iter;  // seed list iterator
    for(sl_iter=Seeds.begin(); sl_iter!=Seeds.end(); sl_iter++)
    {
      vector<Particle>::iterator seed_iter;
      for(seed_iter=sl_iter->begin(); seed_iter!=sl_iter->end(); seed_iter++)
      {
	float time = seed_iter->pt[3];
	if(time >= tmin && time <= tmax)
	{
	  count++;
	}
      }
    }
  }

  return count;
}


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
// convert a VEC data file to a VTK data file (specifically vtkRectilinearGrid)
//
// added by Zhanping Liu on 06/02/2013
//
void	VEC2VTK( char * vec_name, char * vtk_name )
{
	// get the vector data from a BINARY file of type VEC
	int     dataSize[3];
	FILE  *	vec_file = fopen( vec_name, "rb" );
	fread(  dataSize,  sizeof( int ),  3,  vec_file  );
	float * pVecData = new float [ dataSize[0] * dataSize[1] * dataSize[2] * 3 ];
	fread(  pVecData,  sizeof( float ),  
				       dataSize[0] * dataSize[1] * dataSize[2] * 3,  
				       vec_file  );
	fclose( vec_file );
	

	// create a 3-component vtkFloatArray for the vector data
	vtkFloatArray * vecArray = vtkFloatArray::New();
	vecArray->SetName( "3DTornadoVectors" );
	vecArray->SetNumberOfComponents( 3 );
	vecArray->SetNumberOfTuples( dataSize[0] * dataSize[1] * dataSize[2] );
	float	      * theArray = vecArray->GetPointer( 0 );
	memcpy(  theArray,  pVecData,  
		   sizeof( float ) * dataSize[0] * dataSize[1] * dataSize[2] * 3  );
	delete [] pVecData;
	pVecData = NULL;
	theArray = NULL;

	
	// create a vtkRectilinearGrid from the vtkFloatArray
	int   i;
	vtkFloatArray       * x_coords = vtkFloatArray::New();
	vtkFloatArray       * y_coords = vtkFloatArray::New();
	vtkFloatArray       * z_coords = vtkFloatArray::New();
	x_coords->SetNumberOfComponents( 1 );
	y_coords->SetNumberOfComponents( 1 );
	z_coords->SetNumberOfComponents( 1 );
	x_coords->SetNumberOfTuples( dataSize[0] );
	y_coords->SetNumberOfTuples( dataSize[1] );
	z_coords->SetNumberOfTuples( dataSize[2] );
	for ( i = 0; i < dataSize[0]; i ++ ) x_coords->SetTuple1( i, i );
	for ( i = 0; i < dataSize[1]; i ++ ) y_coords->SetTuple1( i, i );
	for ( i = 0; i < dataSize[2]; i ++ ) z_coords->SetTuple1( i, i );
	
	vtkRectilinearGrid  * rectData = vtkRectilinearGrid::New();
	rectData->SetDimensions  ( dataSize );
	rectData->SetXCoordinates( x_coords );
	rectData->SetYCoordinates( y_coords );
	rectData->SetZCoordinates( z_coords );
	rectData->GetPointData()->SetVectors( vecArray );
	vecArray->Delete();	   vecArray = NULL;
	x_coords->Delete();	   x_coords = NULL;
	y_coords->Delete();	   y_coords = NULL;
	z_coords->Delete();	   z_coords = NULL;


	// save the vtkRectilinearGrid to a VTK file
	vtkRectilinearGridWriter * dataWrtr = vtkRectilinearGridWriter::New();
	dataWrtr->SetFileName(  vtk_name  );
	dataWrtr->SetFileType( VTK_BINARY );
	dataWrtr->SetInput( rectData );
	dataWrtr->Update();
	dataWrtr->Delete();	   dataWrtr = NULL;
	rectData->Delete();	   rectData = NULL;
}


// -----------------------------------------------------------------------
//
// convert a series of VEC data files to VTK data files (of type vtkRectilinearGrid)
//
// thePrefx: common part of the VEC file names (including the path and symbol '/')
//
// added by Zhanping Liu on 06/02/2013
//
void	TornadoDataFiles2VTK( char * thePrefx, int fileIdx0, int fileIdx1 )
{
	char	vec_name[200];
	char	vtk_name[200];

	for ( int i = fileIdx0; i <= fileIdx1; i ++ )
	{
		sprintf( vec_name, "%s%d.vec",   thePrefx, i );
		sprintf( vtk_name, "%s%02d.vtk", thePrefx, i );
		VEC2VTK( vec_name, vtk_name );
	}
}


// -----------------------------------------------------------------------
// 
// save a rectilinear grid dataset as a curvilinear grid RAW dataset
// (prerequisite: a rectilinear grid VECTOR dataset HAS BEEN loaded into memory)
//
// pathName: a path name but WITHOUT "/" at the end
// mainName: a file name but WITHOUT the path and the extension
//
// added by Zhanping Liu on 06/17/2013
//
void	Rectilinear2CurvilinearRaw
	( int * gridSizs, float * gridVecs, char * pathName, char * mainName )
{
	int	  araySize = gridSizs[0] * gridSizs[1] * gridSizs[2] * 3;


	// create the coordinates for the grid
	float  *  theCords = new float [ araySize ];
	float  *  thisCord = theCords;
	for ( int k = 0; k < gridSizs[2]; k ++ )
	for ( int j = 0; j < gridSizs[1]; j ++ )
	for ( int i = 0; i < gridSizs[0]; i ++, thisCord += 3 )
	{	
		thisCord[0] = i;
		thisCord[1] = j;
		thisCord[2] = k;
	}
	thisCord = NULL;


	// save the grid (coordinates)
	char	  gridFile[200];
	sprintf(  gridFile,  "%s/%s.grd",  pathName, mainName  );
	FILE   *  out_file = fopen( gridFile, "wb" );
	fwrite (  gridSizs,  sizeof(  int  ),  3,         out_file  );
	fwrite (  theCords,  sizeof( float ),  araySize,  out_file  );
	fclose (  out_file  );	out_file = NULL;
	delete [] theCords;	theCords = NULL;


	// save the vector data
	char	  vec_file[200];
	sprintf(  vec_file,  "%s/%s.vec",  pathName, mainName  );
	out_file = fopen( vec_file, "wb" );
	fwrite (  gridSizs,  sizeof(  int  ),  3,         out_file  );
	fwrite (  gridVecs,  sizeof( float ),  araySize,  out_file  );
	fclose (  out_file  );	out_file = NULL;
	     
 
	// generate a list file for the curvilinear RAW dataset
	char	  listFile[200];
	sprintf(  listFile,  "%s/%s.list",  pathName,  mainName  );
	out_file = fopen( listFile, "w" );
	fprintf(  out_file, "%d %d %d 1 # CURVILINEAR grid size (x, y, z, t; t: number of time steps)\n", 
		            gridSizs[0], gridSizs[1], gridSizs[2]  );
	fprintf(  out_file, "1.0 # vector scaling factor\n"  );
	fprintf(  out_file, "# %s.grd: the grid file name, do NOT turn on this line.\n",  mainName  );
	fprintf(  out_file, "%s.vec # list of vector file names, one per time step\n",    mainName  );
	fclose (  out_file  );	out_file = NULL;
}


// -----------------------------------------------------------------------
//
// save a 3D curvilinear grid VECTOR dataset to a vtkStructuredGrid file
//
// NOTE: curv_grd and curv_vec have the exact number of floats, which is
//       dimensns[0] * dimensns[1] * dimensns[2] * 3
//
// added by Zhanping Liu on 06/19/2013
//
void	CurvilinearGridData2VTK
	( int   * dimensns, float * curv_grd, 
	  float * curv_vec, char  * fileName, int fileType )
{
	int	    num_pnts = dimensns[0] * dimensns[1] * dimensns[2];


	// create a vtkPoints object and fill it with the grid (coordinates)
	float     * this_pnt = curv_grd;
	vtkPoints * the_pnts = vtkPoints::New();
	the_pnts->SetNumberOfPoints( num_pnts );
	for ( int i = 0; i < num_pnts; i ++, this_pnt += 3 )
	the_pnts->SetPoint( i, this_pnt );

	this_pnt = NULL;


	// create a 3-component point attribute for the array of 3D vectors 
	float	      * vecArray = NULL;
	vtkFloatArray * ptAttrib = vtkFloatArray::New();
	ptAttrib->SetName( "VelocityVectors" );
	ptAttrib->SetNumberOfComponents( 3 );	 // 3D vectors
	ptAttrib->SetNumberOfTuples( num_pnts );
	vecArray = (  float *  )  (  ptAttrib->GetPointer( 0 )  );
	memcpy(  vecArray,  curv_vec,  sizeof( float ) * 3 * num_pnts  );

	vecArray = NULL;


	// create a vtkStructuredGrid object and attach the vtkPoints 
	// as well as the vtkFloatArray (a point attributte) to it
	vtkStructuredGrid * strctGrd = vtkStructuredGrid::New();
	strctGrd->SetDimensions( dimensns );
	strctGrd->SetPoints( the_pnts );
	strctGrd->GetPointData()->SetVectors( ptAttrib );

	the_pnts->Delete();    the_pnts = NULL;
	ptAttrib->Delete();    ptAttrib = NULL;


	// write this vtkStructuredGrid object to a VTK file
	vtkStructuredGridWriter * sgWriter = vtkStructuredGridWriter::New();
	sgWriter->SetFileName( fileName );
	sgWriter->SetFileType( fileType );
	sgWriter->SetInput   ( strctGrd );
	sgWriter->Update();

	sgWriter->Delete();    sgWriter = NULL;
	strctGrd->Delete();    strctGrd = NULL;
}


// -----------------------------------------------------------------------
//
// Load a 3D Cartesian dataset in raw mode (by process #0 only)
//
// vec_buff: already allocated with the exact size
//
// added by Zhanping Liu on 05/23/2013
//
void	LoadCartesianRawData( char * raw_file, float * vec_buff )
{	
	int  thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );
	
	if ( thisProc == 0 )
	{			
		FILE  * diskFile = fopen( raw_file, "rb" );
		fread ( vec_buff, sizeof( float ), 
			int( size[0] ) * int( size[1] ) * int( size[2] ) * 3, 
			diskFile );
		fclose( diskFile );
	}
	
	MPI_Barrier( MPI_COMM_WORLD );
}


// -----------------------------------------------------------------------
//
// read one time step of a vtkRectilinearGrid dataset and use it as a 
// Cartesian dataset (invoked by process #0 only)
//
// vec_buff: already allocated with the exact size
//
// added by Zhanping Liu on 05/23/2013
//
int	LoadRectilinearVTKdata( char * rectFile, float * vec_buff )
{
	int  thisProc, numProcs, retValue = 1;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

	if ( thisProc == 0 )
	{	
		cout << "  reading vtkRectilinearGrid data file " 
		     << rectFile << " ...... " << endl;
		retValue = 0;
		vtkRectilinearGridReader * rgReader = vtkRectilinearGridReader::New();
		vtkRectilinearGrid       * rectGrid = vtkRectilinearGrid::New();
		rgReader->SetFileName( rectFile );
		rgReader->SetOutput  ( rectGrid );
		rgReader->Update();
		rgReader->Delete();


		int     i, j;
		int  *  dimensns = rectGrid->GetDimensions();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     data  dimensions: "; 
		for ( i = 0; i < 3; i ++ ) 
		{ cout << dimensns[i]; if ( i + 1 < 3 ) cout << " x "; } cout << endl;
		#endif
	

		vtkPointData * pnt_data = rectGrid->GetPointData();
		int	       numArays = pnt_data->GetNumberOfArrays();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     number of point data arrays: " << numArays << endl;
		#endif

		float         * tmp_data = NULL;
		vtkDataArray  * dataAray = NULL;
		vtkFloatArray * vecArray = NULL;
		if (  (  numArays     >   0               )  &&  
                      (  dimensns[0]  ==  int( size[0] )  )  &&
		      (  dimensns[1]  ==  int( size[1] )  )  &&
                      (  dimensns[2]  ==  int( size[2] )  )
                   )	
		{	
			j = -1;
			for ( i = 0; i < numArays; i ++ )
			{
				dataAray = pnt_data->GetArray( i );
				#ifdef _DUMP_VTK_FILE_INF_
				cout << "     point data array #" << i << ":  name = " 
			     	     << dataAray->GetName() 
			     	     << ";  number of components = " 
			     	     << dataAray->GetNumberOfComponents() << endl; 
				#endif
				if ( dataAray->GetNumberOfComponents() == 3 ) { j = i; break; }
			}

			if ( j >= 0 )
			{
				#ifdef _DUMP_VTK_FILE_INF_
				cout << "     point data array #" << j << " is a 3D vector" << endl;
				#endif
				vecArray = vtkFloatArray::SafeDownCast( dataAray );
				tmp_data = vecArray->GetPointer(0);
			
				// save the 3D Cartesian vector data to a disk file (for test purposes)
				#ifdef _GENERATE_RAW_DATA_ 
				FILE * diskFile = fopen( "vtkDataRaw.list", "w" );
				fprintf( diskFile, "%d %d %d 1 # grid size (x, y, z, t; t is number of timesteps)\n", 
				 	 dimensns[0], dimensns[1], dimensns[2] );
				fprintf( diskFile, "1.0 # vector scaling factor\n" );
				fprintf( diskFile, "vtkDataRaw.dat #list of file names, one per timestep\n" );
				fclose( diskFile );

				diskFile = fopen( "vtkDataRaw.dat", "wb" );
				fwrite( tmp_data, sizeof( float ), 
					dimensns[0] * dimensns[1] * dimensns[2] * 3, diskFile );
				fclose( diskFile );
				#endif


				memcpy(  vec_buff,  tmp_data, 
				 	 sizeof( float ) * dimensns[0] * dimensns[1] * dimensns[2] * 3  );
				retValue = 1;

				#ifdef _RECT_TO_CURV_DATA_
				if ( tsize == 1 ) Rectilinear2CurvilinearRaw( dimensns, vec_buff, 
						  "/home/zhanping/Integration/Data", "curvGridRaw" );
				#endif
			}
		}

		dimensns = NULL;
		pnt_data = NULL;
		dataAray = NULL;
		vecArray = NULL;
		tmp_data = NULL;
		rectGrid->Delete();

		cout << "  ...... done." << endl << endl;
	}

	MPI_Barrier( MPI_COMM_WORLD );

	return  retValue;
}


// -----------------------------------------------------------------------
//
// load the (point) coordinates of a curvilinear grid
//
// invoked by process #0 only, once for ALL time steps
//
// grd_buff: to be (re-)allocated here inside this function
//
// added by Zhanping Liu on 06/18/2013
//
void  	LoadCurvilinearRawGrid( char * gridFile, float * & grd_buff )
{
	cout << endl << "==== loading a curvilinear RAW grid ====" 
	     << endl << "     " << gridFile << endl;

	if ( grd_buff )  delete [] grd_buff;   grd_buff = NULL;

	int	dimensns[3];
	int	araySize = 0;
	FILE  * cordFile = fopen( gridFile, "rb" );
	fread(  dimensns,  sizeof( int ),  3,  cordFile  );
	assert( dimensns[0] == int( size[0] ) && 
	        dimensns[1] == int( size[1] ) && 
		dimensns[2] == int( size[2] ) );
	araySize = dimensns[0] * dimensns[1] * dimensns[2] * 3;
	grd_buff = new float [ araySize ];
	fread(  grd_buff,  sizeof( float ),  araySize,  cordFile  );
	fclose( cordFile );    cordFile = NULL;

	cout << endl << ".... curvilinear RAW grid loaded." << endl;
}


// -----------------------------------------------------------------------
//
// load the VECTOR (or 'SOLUTION') DATA, without the grid (coordinates), 
// defined on a curvilinear grid
//
// invoked by process #0 only, one for each time step
//
// vec_buff: has been allocated with the exact size by the caller
//
// added by Zhanping Liu on 06/19/2013
//
void  	LoadCurvilinearRawData( char * vec_file, float * vec_buff )
{
	int  thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

	if ( thisProc == 0 )
	{
		cout << endl << "==== loading a curvilinear RAW vector dataset ====" 
	     	     << endl << "     " << vec_file << endl;
	
		int     dimensns[3];
		FILE  *	dataFile = fopen( vec_file, "rb" );
		fread(  dimensns,  sizeof( int ),  3,  dataFile  );
		assert( dimensns[0] == int( size[0] ) &&
	  		dimensns[1] == int( size[1] ) &&
			dimensns[2] == int( size[2] ) );
		fread(  vec_buff,  sizeof( float ),
			dimensns[0] * dimensns[1] * dimensns[2] * 3,  dataFile  );
		fclose( dataFile );   dataFile = NULL;

		cout << endl << ".... curvilinear RAW vector dataset loaded." << endl;

		#ifdef _CURVILINEAR_2_VTK_

		// save the curvilinear dataset to a VTK file (of type vtkStructuredGrid)
		CurvilinearGridData2VTK( dimensns, grd_cord, vec_buff,
	                                 "/home/zhanping/Integration/Data/curvGridVTK.vtk" );

		// create a list file for this VTK file
		FILE   *  listFile = fopen( "/home/zhanping/Integration/Data/curvGridVTK.list", "w" );
		fprintf(  listFile,  "%d %d %d 1 # CURVILINEAR grid size (x, y, z, t; t: number of time steps)\n", 
		                     dimensns[0], dimensns[1], dimensns[2]  );
		fprintf(  listFile,  "1.0 # vector scaling factor\n"  );
		fprintf(  listFile,  "curvGridVTK.vtk # list of vector file names, one per time step\n" );
		fclose (  listFile  );	listFile = NULL;
		
		#endif
	}

	MPI_Barrier( MPI_COMM_WORLD );
}


// -----------------------------------------------------------------------
//
// load the (point) coordinates of a curvilinear grid
//
// invoked by process #0 only, once for ALL time steps
//
// grd_buff: to be (re-)allocated here inside this function
//
// added by Zhanping Liu on 06/19/2013
//
void	LoadCurvilinearVtkGrid( char * curvFile, float * & grd_buff )
{
	cout << endl << "==== loading a curvilinear VTK grid ====" 
	     << endl << "     " << curvFile << endl;

	if ( grd_buff )  delete [] grd_buff;   grd_buff = NULL;


	// obtain a vtkStructuredGrid object from the file
	vtkStructuredGrid       * curvGrid = vtkStructuredGrid::New();
	vtkStructuredGridReader * sgReader = vtkStructuredGridReader::New();
	sgReader->SetFileName( curvFile );
	sgReader->SetOutput  ( curvGrid );
	sgReader->Update();
	sgReader->Delete();    sgReader = NULL;


	// dump some grid size info
	int     i;
	int  *  dimensns = curvGrid->GetDimensions();
	#ifdef _DUMP_VTK_FILE_INF_
	cout << "     data  dimensions: "; 
	for ( i = 0; i < 3; i ++ ) 
	{ cout << dimensns[i]; if ( i + 1 < 3 ) cout << " x "; } cout << endl;
	#endif


	// make sure the grid size info is correct
	assert( dimensns[0] == int( size[0] ) && 
	        dimensns[1] == int( size[1] ) && 
		dimensns[2] == int( size[2] ) );


	// allocate the grid buffer and fill it with the point coordinates
	int  	    num_pnts = dimensns[0] * dimensns[1] * dimensns[2];
	float     * that_pnt = NULL;	// the destination
	double    * this_pnt = NULL;	// the source
	vtkPoints * the_pnts = curvGrid->GetPoints();
	grd_buff  = new float [ num_pnts * 3 ];
	for ( i = 0, that_pnt = grd_buff; i < num_pnts; i ++, that_pnt += 3 )
	{
		this_pnt = the_pnts->GetPoint( i );
		that_pnt[0] = float( this_pnt[0] );
		that_pnt[1] = float( this_pnt[1] );
		that_pnt[2] = float( this_pnt[2] );
	}

	this_pnt = NULL;             that_pnt = NULL;	      
	the_pnts = NULL;	     dimensns = NULL;
	curvGrid->Delete();	     curvGrid = NULL;
	
	cout << endl << ".... curvilinear VTK grid loaded." << endl;
}


// -----------------------------------------------------------------------
//
// load vector data from a vtk curvilinear grid (vtkStructuredGrid) file
//
// invoked by process #0 only, once for ALL time steps
//
// vec_buff: has been allocated with the exact size by the caller 
//
// added by Zhanping Liu on 06/20/2013
//
int	LoadCurvilinearVtkData( char * curvFile, float * vec_buff )
{
	int  thisProc, numProcs, retValue = 1;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );


	if ( thisProc == 0 )
	{
		retValue = 0;
		cout << endl << "==== loading a curvilinear VTK data ====" 
	     	     << endl << "     " << curvFile << endl;


		// obtain a vtkStructuredGrid object from the file
		vtkStructuredGrid       * curvGrid = vtkStructuredGrid::New();
		vtkStructuredGridReader * sgReader = vtkStructuredGridReader::New();
		sgReader->SetFileName( curvFile );
		sgReader->SetOutput  ( curvGrid );
		sgReader->Update();
		sgReader->Delete();    sgReader = NULL;


		int     i, j;
		int  *  dimensns = curvGrid->GetDimensions();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     data  dimensions: "; 
		for ( i = 0; i < 3; i ++ ) 
		{ cout << dimensns[i]; if ( i + 1 < 3 ) cout << " x "; } cout << endl;
		#endif
	

		vtkPointData * pnt_data = curvGrid->GetPointData();
		int	       numArays = pnt_data->GetNumberOfArrays();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     number of point data arrays: " << numArays << endl;
		#endif

		float         * tmp_data = NULL;
		vtkDataArray  * dataAray = NULL;
		vtkFloatArray * vecArray = NULL;
		if (  (  numArays     >   0               )  &&  
                      (  dimensns[0]  ==  int( size[0] )  )  &&
		      (  dimensns[1]  ==  int( size[1] )  )  &&
                      (  dimensns[2]  ==  int( size[2] )  )
                   )	
		{	
			j = -1;
			for ( i = 0; i < numArays; i ++ )
			{
				dataAray = pnt_data->GetArray( i );
				#ifdef _DUMP_VTK_FILE_INF_
				cout << "     point data array #" << i << ":  name = " 
			     	     << dataAray->GetName() 
			     	     << ";  number of components = " 
			     	     << dataAray->GetNumberOfComponents() << endl; 
				#endif
				if ( dataAray->GetNumberOfComponents() == 3 ) { j = i; break; }
			}

			if ( j >= 0 )
			{
				#ifdef _DUMP_VTK_FILE_INF_
				cout << "     point data array #" << j << " is a 3D vector" << endl;
				#endif
				vecArray = vtkFloatArray::SafeDownCast( dataAray );
				tmp_data = vecArray->GetPointer(0);


				memcpy(  vec_buff,     tmp_data, 
				 	 sizeof( float ) * 
				         dimensns[0] * dimensns[1] * dimensns[2] * 3  );
				retValue = 1;
			}
		}

		dimensns = NULL;
		pnt_data = NULL;
		tmp_data = NULL;
		dataAray = NULL;
		vecArray = NULL;
		curvGrid->Delete();

		cout << endl << ".... curvilinear VTK data loaded." << endl;
	}

	MPI_Barrier( MPI_COMM_WORLD );

	return retValue;
}


// -----------------------------------------------------------------------
//
// load a vtk curvilinear grid (vtkStructuredGrid) dataset that contains
// multiple INTERMEDIATE point attributes (density, 3-component momentum, 
// and energy) from which the 3D velocity vector is derived
//
// invoked by process #0 only, once for ALL time steps
//
// vec_buff: has been allocated with the exact size by the caller 
//
// added by Zhanping Liu on 06/20/2013
//
void	LoadCurvilinearVtkData_DME( char * curvFile, float * vec_buff )
{
	int  thisProc, numProcs, retValue = 1;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );


	if ( thisProc == 0 )
	{
		cout << endl << "==== loading a curvilinear VTK data ====" 
	     	     << endl << "     " << curvFile << endl;


		// obtain a vtkStructuredGrid object from the file
		vtkStructuredGrid       * curvGrid = vtkStructuredGrid::New();
		vtkStructuredGridReader * sgReader = vtkStructuredGridReader::New();
		sgReader->SetFileName( curvFile );
		sgReader->SetOutput  ( curvGrid );
		sgReader->Update();
		sgReader->Delete();    sgReader = NULL;


		// dump some data size info
		int            i;
		int       *    dimensns = curvGrid->GetDimensions();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     data  dimensions: "; 
		for ( i = 0; i < 3; i ++ ) 
		{ cout << dimensns[i]; if ( i + 1 < 3 ) cout << " x "; } cout << endl;
		#endif


		// make sure the data size size info is correct
		assert( dimensns[0] == int( size[0] ) && 
	        	dimensns[1] == int( size[1] ) && 
			dimensns[2] == int( size[2] ) );


		// access to the point data
		int	       denstyId = 0;	// index of the 'density'  attribute
		int	       momntmId = 0;	// index of the 'Momentum' attribute
		vtkDataArray * dataAray = NULL;
		vtkPointData * pnt_data = curvGrid->GetPointData();
		

		// dump the number of point attributes
		int	       numArays = pnt_data->GetNumberOfArrays();
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     number of point data arrays: " << numArays << endl;
		#endif
	

		// collect the meta info of each point attribute
		for ( i = 0; i < numArays; i ++ )
		{
			dataAray = pnt_data->GetArray( i );
			#ifdef _DUMP_VTK_FILE_INF_
			cout << "     point data array #" << i << ":  name = " 
			     << dataAray->GetName() 
			     << ";  number of components = " 
			     << dataAray->GetNumberOfComponents() << endl;
			#endif
			if (    strcmp(  dataAray->GetName(),  "Density"  )  
			        ==    0    )     		      denstyId = i;
			else
			if ( dataAray->GetNumberOfComponents() == 3 ) momntmId = i;
		}

			
		// identify the 'Density' and 'Momentum' attributes	
		#ifdef _DUMP_VTK_FILE_INF_
		cout << "     point data array #" << denstyId
		     << " is a 1D scalar attribute called 'Density' " << endl;
		if (  pnt_data->GetArray( denstyId )->IsA( "vtkFloatArray" )  )
		cout << "                     type: vtkFloatArray"    << endl;
		else
		cout << "                     type: vtkDoubleArray"   << endl;
		cout << "     point data array #" << momntmId 
		     << " is a 3D vector attribute called 'Momentum'" << endl;
		if (  pnt_data->GetArray( momntmId )->IsA( "vtkFloatArray" )  )
		cout << "                     type: vtkFloatArray"    << endl;
		else
		cout << "                     type: vtkDoubleArray"   << endl;
		#endif

		// derive the velocity vector
		int  	   num_pnts = dimensns[0] * dimensns[1] * dimensns[2];
		float	   ytisned1 = 1.0;     // reverse the letters: 1.0 / density
		float   *  this_vec = vec_buff;
		float	*  pDensity = ( float * ) 
		(    vtkFloatArray::SafeDownCast(  pnt_data->GetArray( denstyId )  )
		   ->GetPointer( 0 )    );
		float   *  pMomntum = ( float * )
		(    vtkFloatArray::SafeDownCast(  pnt_data->GetArray( momntmId )  )
		   ->GetPointer( 0 )    );

		for ( i = 0; i < num_pnts; 
		      i ++,  this_vec += 3, pDensity ++, pMomntum += 3 )
		{
			ytisned1    = 1.00000 / ( * pDensity );
			this_vec[0] = pMomntum[0] * ytisned1;
			this_vec[1] = pMomntum[1] * ytisned1;
			this_vec[2] = pMomntum[2] * ytisned1;
		}


		#ifdef _CURVILINEAR_2_VTK_

		// save the curvilinear dataset to a VTK file (of type vtkStructuredGrid)
		CurvilinearGridData2VTK( dimensns, grd_cord, vec_buff,
	                                 "/home/zhanping/Integration/Data/curvGrid_bluntFin.vtk" );

		// create a list file for this VTK file
		FILE   *  listFile = fopen( "/home/zhanping/Integration/Data/curvGrid_bluntFin.list", "w" );
		fprintf(  listFile,  "%d %d %d 1 # CURVILINEAR grid size (x, y, z, t; t: number of time steps)\n", 
		                     dimensns[0], dimensns[1], dimensns[2]  );
		fprintf(  listFile,  "1.0 # vector scaling factor\n"  );
		fprintf(  listFile,  "curvGrid_bluntFin.vtk # list of vector file names, one per time step\n" );
		fclose (  listFile  );	listFile = NULL;
		
		#endif
	
		
		// clean up
		dimensns = NULL;		    pnt_data = NULL;
		dataAray = NULL;		    this_vec = NULL;
		pDensity = NULL;		    pMomntum = NULL;
		curvGrid->Delete();	            curvGrid = NULL;

		cout << endl << ".... curvilinear VTK data loaded." << endl;
	}

	MPI_Barrier( MPI_COMM_WORLD );
}


// -----------------------------------------------------------------------
//
// given a 3-component volume data, extract a sub-volume to fill in a block
// 
// 3-component: either a 3D vector or a 3D coordinate (for curvilinear grids)
//
// invoked by process #0 only
//
// added by Zhanping Liu on 05/24/2013
//
void	FillBlock( int * volSizes, float * vol_data, int * blckMins, 
		   int * blkSizes, float * blckData )
{	
	int	sliceSiz = volSizes[0] * volSizes[1];	
	float *	des_data = blckData;
	float *	src_data = NULL;

	for ( int k = 0; k < blkSizes[2]; k ++ )
	{	
		int	z_offset = sliceSiz * ( blckMins[2] + k );
		for ( int j = 0; j < blkSizes[1]; j ++ )
		{
			int	y_offset = volSizes[0] * ( blckMins[1] + j );
			for ( int i = 0; i < blkSizes[0]; i ++, des_data += 3 )
			{
				src_data    = vol_data + 
					    ( z_offset + y_offset + blckMins[0] + i ) * 3;
				des_data[0] = src_data[0];
				des_data[1] = src_data[1];
				des_data[2] = src_data[2];
			}
		}
	}

	src_data = NULL;
	des_data = NULL;
}


//-----------------------------------------------------------------------
//
// init the spatial-temporal (GRID and DATA) blocks distributor by determining
// some global information that keeps unchanged throughout the parallel process
//
// also allocated are some arrays of which the sizes keep unchanged
//
// this process avoids frequent on-the-fly determination of some information
//
// added by Zhanping Liu on 06/02/2013 and last updated on 06/03/2013
//
void 	InitBlocksDistributor()
{
	int   	i, j, thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );


	// determine the max number of time steps of a single time group
  	for ( i = 0; i < nblocks; i ++ )
  	{
		int	 blockInf[8];
		DIY_Block_starts_sizes( 0, i, blockInf, blockInf + 4 );
		max_tblk = ( max_tblk >= blockInf[7] ) ? max_tblk : blockInf[7];
  	}


	// get the max number of spatial-temporal blocks on a single process
	MPI_Reduce( &nblocks,  &maxBlcks, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
	MPI_Bcast ( &maxBlcks, 1, MPI_INT, 0, MPI_COMM_WORLD );	

	
	// since LoadRectilinearVTKdata() references vec_data[i], we need to
	// allocate (the 1st-level array of) vec_data for ALL processes!!!
	//
  	vec_data = new  float * [ max_tblk ];
  	for ( i = 0; i < max_tblk; i ++ ) vec_data[i] = NULL;


	// only process #0 needs to actually load and maintain the volume data
  	if ( thisProc == 0 )
  	{
		for ( i = 0; i < max_tblk; i ++ )
  		vec_data[i] = new  float  
		[  int( size[0] )  *  int( size[1] )  *  int( size[2] )  *  3  ];
  	}


	// status of each spatial-temporal block on a single process
	// 0: inactive (not used by the 'current' time group)
	// 1:   active (    used by the 'current' time group)
	//
	beActive = new int [ maxBlcks ];  


  	// get the bounding-box info of each block on THIS process
	// (each entry: 4 starts or mins followed by 4 sizes)
	bb_infor = new int [ 8 * maxBlcks ];
	memset(  bb_infor,  0,  sizeof( int ) * 8 * maxBlcks  );
  	for ( i = 0; i < nblocks; i ++ )
  	{
  		DIY_Block_starts_sizes( 0, i, bb_infor + 8 * i, bb_infor + 8 * i + 4 );
  	}


  	// the number of blocks  on each  process --- collected and maintained by process #0
	// ALL bounding-box info of ALL processes --- collected and maintained by process #0
	// status of ALL spatial-temporal blocks on ALL processes:  maintained by process #0
  	if ( thisProc == 0 )
  	{
		numBlcks = new int [ numProcs ];
		BEACTIVE = new int [ numProcs * maxBlcks     ];
		BB_INFOR = new int [ numProcs * maxBlcks * 8 ];
  	}

  	MPI_Gather( &nblocks, 1,            MPI_INT, numBlcks, 1,            MPI_INT, 0, MPI_COMM_WORLD );
  	MPI_Gather( bb_infor, 8 * maxBlcks, MPI_INT, BB_INFOR, 8 * maxBlcks, MPI_INT, 0, MPI_COMM_WORLD );

	
	// get the time group index of each spatial-temporal block of THIS process
	tmGrpIds = new int [ maxBlcks ]; 
	memset(  tmGrpIds,  0,  sizeof( int ) * maxBlcks  );
	for ( i = 0; i < nblocks; i ++ )  DIY_In_time_block( 0, i, tmGrpIds + i );
	

	// allocate two arrays for maintaining the following information
	// 
	// tBlkMins[i]: index  of the first time step  of time group #i
	// tBlkSizs[i]: number of the       time steps of time group #i
	//
	// this information is the same across ALL processes  and hence
	// only a single process (#0) may be needed to obtain it before
	// broadcasing it to the other processes, while a parallel mode
	// is adopted here to avoid broadcast
	//
	tBlkMins = new int [ ntpart ];
	tBlkSizs = new int [ ntpart ];
	for ( i = 0; i < ntpart; i ++ )
	{
		// find any one block ('key block') of the time group
		//
		// the iterator below (j) MIGHT be incremented as different
		// time group indices are retrieved, while a conservative
		// scheme is used here, just in case the block decomposition
		// method is arbitrary, e.g., local block #n within time group
		// ( m + 1 ) v.s. local block #( n + 1 ) within time group
		// ( m - 1 )
		//
		int   keyBlkId = 0;
		for ( j = 0; j < nblocks; j ++ )
		{
			if ( tmGrpIds[j] == i )   {  keyBlkId = j;  break;  }	
		}

		// obtain the global index of the first time step of this time
 		// group and the number of time steps of this time group
		tBlkMins[i] = bb_infor[ keyBlkId * 8 + 3 ];
		tBlkSizs[i] = bb_infor[ keyBlkId * 8 + 7 ];
	}


	// load the coordinates of a curvilinear grid (for process #0 only)
	//
	// allocate grid (coordinates) blocks and establish the indexing scheme
	// for correspondance between SPATIAL blocks and SPATIAL-TEMPORAL blocks
	// (for ALL processes)
	//
	// note: the issue of 'nblocks', as annotated in DistributeDataBlocks( ... ),
	//       incurs an overlap / confusion between the spatial dimension and the
	//       temporal dimension
	//	
	//       a single copy of the grid volume (of coordinates) is maintained
	//	 (for curvilinear grid data) even though 'nblocks' is used herein 
	//	 (for consistency purposes)
	//       
	// executed once, even for a time-varying curvilinear grid dataset
	//
	// added by Zhanping Liu on 06/18/2013
	//
	#ifdef _CURVILINEAR_GRIDS_
	
	// load the coordinates and allocate memory for STBLKIDS: for process #0
	if ( thisProc == 0 )
	{
		#ifdef _LOAD_VTK_DATA_SET_
		LoadCurvilinearVtkGrid( dataset_files[0], grd_cord );
		#else
		// derive the grid file name based on the vector file name
		char	gridFile[200];
		char  * ext_name = gridFile + strlen( dataset_files[0] ) - 3;
		strcpy( gridFile,  dataset_files[0] );
		strcpy( ext_name, "grd" );    ext_name = NULL;  // "vec" to "grd"

		// load the whole grid (set of coordinates)
		LoadCurvilinearRawGrid( gridFile, grd_cord );
		//LoadCurvilinearVtkGrid( "/home/zhanping/Integration/Data/curvGridVTK.vtk", grd_cord );

		#endif

		// allocate STBLKIDS: maintained by process #0 but collected from ALL processes
		// index of the FIRST corresponding spatial-temporal block of each SPATIAL block
		STBLKIDS = new int [ nspart ]; 
	}

	// allocate grid (coordinates) blocks: for ALL processes
	int
	numRegds = 0;		     	     // number of the SPATIAL blocks that have been registered
	nsBlocks = nspart    /   numProcs;   // derived here to avoid any change to GetArgs( ... ) ^_^
	sBlckIds = new int     [ nblocks  ]; // SPATIAL index of each spatial-temporal block
	stBlkIds = new int     [ nsBlocks ]; // the FIRST spatial-temporal index of each SPATIAL block
	gridBlks = new float * [ nsBlocks ]; // only nsBlocks SPATIAL (coordinates) blocks are needed
	
	for ( i = 0; i < nblocks; i ++ )
	{
		int	targetId = -1;	     // ID   of the target SPATIAL block that has been registered
		int  *  thatInfo = NULL;     // info of the target SPATIAL block that has been registered
		int  *	thisInfo = bb_infor + 8 * i;// info of the spatial-temporal block being retrieved
		
		// check whether the corresponding SPATIAL block has been registered or not
		for ( j = 0; j < numRegds; j ++ )
		{
			thatInfo = bb_infor + 8 * stBlkIds[j];// turn to the first spatial-temporal block
			if ( thisInfo[0] == thatInfo[0] &&
			     thisInfo[1] == thatInfo[1] &&
			     thisInfo[2] == thatInfo[2] ) // min[0..2]: a signature of each SPATIAL block
			{
				targetId = j;		  // found  the SPATIAL block  registered therein
				break;
			}
		}
		thatInfo = NULL;

		// to add a link to an existing SPATIAL block or register a new one		
		if ( targetId >= 0 )
		{
			sBlckIds[ i        ] = targetId;	// SPATIAL-TEMPORAL ---> SPATIAL
		}
		else
		{
			stBlkIds[ numRegds ] = i;	  	// SPATIAL ---> SPATIAL-TEMPORAL
			sBlckIds[ i        ] = numRegds;	// SPATIAL-TEMPORAL ---> SPATIAL
			gridBlks[ numRegds ] = new float 	// allocate memory for the SPATIAL block
			        [ thisInfo[4] * thisInfo[5] * thisInfo[6] * 3 ];
			numRegds  ++;
		}
		thisInfo = NULL;
	}

	
	// gather the SPATIAL ---> SPATIAL-TEMPORAL map info across ALL processes
	MPI_Gather( stBlkIds, nsBlocks, MPI_INT, STBLKIDS, nsBlocks, MPI_INT, 0, MPI_COMM_WORLD );
	#endif

	
	// an array of 2D spatial-temporal data-block pointers for THIS process
	dataBlks = new float ** [ nblocks ];
	for ( i = 0; i < nblocks; i ++ ) dataBlks[i] = NULL;
}


// -----------------------------------------------------------------------
//
// distribute DIY GRID blocks from the master process to the slave processes
// (for curvilinear grids of which the coordinates are involved & distributed)
//
// this function is invoked only once, AFTER InitBlocksDistributor( ... ),
// even for a time-varying dataset
//
// added by Zhanping Liu on 06/18/2013
//
void 	DistributeGridBlocks()
{
	int     thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

	
	if ( thisProc == 0 )
	cout << "==== distributing DIY (curvilinear) grid blocks ====" << endl;


	int     i, j;
	int  *  blockInf    = NULL;
	int	volSizes[3] = {  int( size[0] ),  int( size[1] ),  int( size[2] )  };


	// directly fill in the blocks for process #0 only, without inter-process send-receive
	if ( thisProc == 0 )
	{
		for ( i = 0; i < nsBlocks; i ++ )
		{
			blockInf =  bb_infor + stBlkIds[i] * 8;
			FillBlock(  volSizes,  grd_cord,  
				    blockInf,  blockInf + 4,  gridBlks[i]  );
		}
	}


	// fill in the blocks for the slave processes
	//
	// msg_tags, sendReqs, and sendBufs are all excessively allocated  just to
	// make things simple, whereas sendBufs  does NOT waste memory much  since
	// each of its entries (initialized with NULL) is allocated ONLY ON DEMAND
	//
	// msg_tags and sendBufs both involve exact send-receive correspondance
	//
	// in fact, ( numProcs * nsBlocks )  ==  nspart, while the former is intuitive ^_^
	//
	int           commFlag = 0; // whether a communication is completed (1) or not (0)
	int           flag_sum = 0; // number of the requests that have been COMPLETED
	int	      num_reqs = 0; // number of the requests that have been ISSUED
	int	    * msg_tags = new  int  [ numProcs * nsBlocks ];  // comm IDs
	for ( i = 0;  i < numProcs * nsBlocks; i ++ )   msg_tags[i] = i;
	float      ** sendBufs = NULL;
	MPI_Request * sendReqs = NULL;
	MPI_Request * recvReqs = NULL;
	MPI_Status    comStats;	    // status of a request (PENDING or COMPLETED)
	
	if ( thisProc == 0 )
	{
		// sending blocks ......
		
		sendReqs = new MPI_Request [ numProcs * nsBlocks ];
		sendBufs = new float     * [ numProcs * nsBlocks ];	
		for ( i = 0;  i < numProcs * nsBlocks;  i ++ ) sendBufs[i] = NULL;

		for ( i = 1;  i < numProcs;  i ++ ) // beginning with process #1
        	{	
			// s2stBase: the SPATIAL block ---> SPATIAL-TEMPORAL block map
			int   blckBase = nsBlocks * i;   
			int * s2stBase = STBLKIDS + i * nsBlocks; // nsBlocks IDs for each
			int * infoBase = BB_INFOR + i * maxBlcks * 8;

			for ( j = 0; j < nsBlocks; j ++ )
			{
				blockInf = infoBase  +  s2stBase[j] * 8;
				sendBufs[ blckBase + j ] =
				new float [ blockInf[4] * blockInf[5] * blockInf[6] * 3 ];

				FillBlock(  volSizes, grd_cord,  blockInf,  blockInf + 4,  
					    sendBufs[ blckBase + j ]  );
				MPI_Isend(  sendBufs[ blckBase + j ],
					    blockInf[4] * blockInf[5] * blockInf[6] * 3,
					    MPI_FLOAT,   i,  msg_tags[  blckBase + j  ],
					    MPI_COMM_WORLD,  sendReqs + num_reqs  );
				num_reqs ++;
				blockInf = NULL;
			}
			s2stBase = NULL;
			infoBase = NULL;
		}

		// make sure all of the 'num_reqs' GRID blocks are sent
		// out  from the master process  to the slave processes
		while ( flag_sum < num_reqs )
		{
			flag_sum = 0;
			for ( i = 0; i < num_reqs; i ++ )
			{
				MPI_Test( sendReqs + i, &commFlag, &comStats );
				flag_sum += commFlag ? 1 : 0;
			}
		}

		// the while-loop might be substituted by the following line
		// while the former is able to mix computation and communication
		//
		// for ( i = 0; i < num_reqs; i ++ ) MPI_Wait( sendReqs + i, &comStats );
	}
	else
	{
		// receiving blocks ......
		
		int // the first nsBlocks MSGs are virtual and un-used (for process #0)
		msg_base =  thisProc * nsBlocks;
		recvReqs = new MPI_Request [ nsBlocks ];

		for ( j = 0; j < nsBlocks; j ++ )
		{
			blockInf =  bb_infor + stBlkIds[j] * 8;	
			MPI_Irecv(  gridBlks[j],  
				    blockInf[4] * blockInf[5] * blockInf[6] * 3,
				    MPI_FLOAT,   0,  msg_tags[  msg_base + j  ],
				    MPI_COMM_WORLD,  recvReqs + num_reqs  );
			num_reqs ++;
			blockInf = NULL;
		}

		// make sure 'num_reqs' GRID blocks are all received by THIS process
		while ( flag_sum < num_reqs )
		{
			flag_sum = 0;
			for ( i  = 0; i < num_reqs; i ++ )
			{
				MPI_Test( recvReqs + i, &commFlag, &comStats );
				flag_sum += commFlag ? 1 : 0;
			}
		}

		// the while-loop might be substituted by the following line
		// while the former is able to mix computation and communication
		//
		// for ( i = 0; i < num_reqs; i ++ ) MPI_Wait( recvReqs + i, &comStats );
	}
	

	// now grid blocks have been ready on every process (specifically in gridBlks)
	MPI_Barrier( MPI_COMM_WORLD ); 

	
	// clean up
	if ( thisProc == 0 )
	{
		for ( i = 0; i < numProcs * nsBlocks; i ++ )
		{  if ( sendBufs[i] ) delete [] sendBufs[i];  sendBufs[i] = NULL;  }

		delete [] sendBufs;  sendBufs = NULL;
		delete [] sendReqs;  sendReqs = NULL;
  	}
	else	// for each slave process
	{
		delete [] recvReqs;  recvReqs = NULL;
	}

	delete [] msg_tags;  	msg_tags = NULL;

	if ( thisProc == 0 ) 
        cout << ".... DIY (curvilinear) grid blocks distributed." << endl;
}


// -----------------------------------------------------------------------
//
// attach (ghost and real) physical bounds to each DIY GRID block
// 
// this function MUST be called for NON-Cartesian grids  where the physical 
// coordinates (x, y, z) differ from the logical coordinates (computational
// coordinates: i, j, k)
// 
// prior to the call to this function, the (ghost and real) physical bounds
// of each block are 'initialized' with the logical bounds so as to support
// Cartesian grids by default
//
// this function is invoked only once,  AFTER  DistributeGridBlocks( ... ),
// even for a time-varying dataset
//
// added by Zhanping Liu on 07/02/2013
//
void	AttachPhysicalBounds()
{
	// prepare a 2D array of block size information
	int     i;
	int     stBlkIdx = 0;	
	int   * blockInf = NULL;
	int  **	blkSizes = new int * [ nsBlocks ];
	for ( i = 0; i < nsBlocks; i ++ )
	{
		blkSizes[i]    = new int [3];
		stBlkIdx       = stBlkIds[i];
		blockInf       = bb_infor + stBlkIdx * 8;

		blkSizes[i][0] = blockInf[4];
		blkSizes[i][1] = blockInf[5];
		blkSizes[i][2] = blockInf[6];
	}
	blockInf = NULL;

	
	// update the physical bounds with the REAL ones
	// (here REAL is NOT the one relative to ghost)
	DIY_Attach_physical_block_bounds
        ( nsBlocks, blkSizes, gridBlks, stBlkIds, sBlckIds );


	// release memory
	for ( i = 0; i < nsBlocks; i ++ ) 
	{ 
		delete  []  blkSizes[i]; 
		blkSizes[i] = NULL; 
	}
	delete  blkSizes;   blkSizes = NULL; 
}


// -----------------------------------------------------------------------
//
// distribute DIY DATA blocks from the master process to the slave processes
//
// this function is invoked, AFTER InitBlocksDistributor( ... ), for EACH 
// time group (of several time steps)
//
// added by Zhanping Liu on 05/27/2013 and last updated on 06/04/2013
//
void	DistributeDataBlocks( int tmGrpIdx )
{	
	int     thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

	
	if ( thisProc == 0 )
	cout << "==== distributing DIY data blocks for time group #"
	     << tmGrpIdx << " (time steps #" << tBlkMins[tmGrpIdx] << " ~ #" 
	     << tBlkMins[tmGrpIdx] + tBlkSizs[tmGrpIdx] - 1 << ") ====" << endl;
	

	int     i, j, k;
	int	volSizes[3] = {  int( size[0] ),  int( size[1] ),  int( size[2] )  };
	

	// ------------------------------------------------------------------------------ //
	//  										  //
	// 'nblocks' in fact refers to the number of spatial-TEMPORAL blocks of THIS      //
	// process and hence there are multiple such blocks (except for the static case   //
	// where the whole dataset contains only a single time step) with exactly the     //
	// same bounding box but only different time steps (of a single time group or     //
	// in other words called 'time block')                                            //
	//                                                                                //
	// according to the working mechanism of Blocks::BilLoadTimeGroupBlocks( ... )    //
	// and Blocks::LoadBlocks4D( ... ), there are 'nblocks' (1st-level) entries,      //
	// of which each maintains nTmSteps (2nd-level) entries, of which further each    //
	// maintains one (3rd-level) array of data values for a specific SPATIAL-TEMPORAL //
	// block, even though 'nblocks' itself already involves the time dimension        //
	//                                                                                //
	// in other words, there is a redundancy / confusion between 'nblocks' (number of //
	// 1st-level entries) and nTmSteps (number of 2nd-level entries, i.e., the number //
	// of time steps of a time group) that are both used for 'data' --- one argument  //
	// of Blocks::BilLoadTimeGroupBlocks( ... ) and Blocks::LoadBlocks4D( ... )       //
	//                                                                                //
	// to avoid modifications to class Blocks while aseamlessly attaching DIY data    //
	// blocks to the the OSUFlow engine(s), 'nblocks' remains to be used below even   //
	// though it is confusing and inefficient                                         //
	//										  //
	// ------------------------------------------------------------------------------ //


	// detach each 2D (spatial-temporal) data-block pointer ( dataBlks[blockLID] )
	// from Blocks::block_osuflow[blockLID] while the latter releases the memory
	//
	// NOTE: as mentioned above, must NOT release the actual data blocks here !!!
	//       --- "once something is given away, you will not be reponsible for it"
	//
	for ( i = 0; i < nblocks; i ++ ) dataBlks[i] = NULL;


	// ---------------------------------------------------------------------------- //
	// 'maxBlcks' is used below, in place of 'nlocks', to make things simpler while //
	// only a few extra bytes are needed						//
	// ---------------------------------------------------------------------------- //

	
	// clear (re-set) the status of each spatial-temporal block on THIS process
	// (no need to clear BEACTIVE as it will be overwritten by the 'beActive's)
	//
	// fill in this array with the up-to-date status based on the THIS time group
	//
	// allocate 2nd-level (2D) and 3rd-level (1D) arrays for the ACTIVE blocks 
	//
	memset(  beActive,  0,  sizeof( int ) * maxBlcks  );
	int	 nTmSteps = tBlkSizs[ tmGrpIdx ]; // number of time steps of THIS group
	int  *	 blockSiz = NULL;		  // accessing the spatial size of blocks
	for ( i = 0; i < nblocks; i ++ )
	{
		if ( tmGrpIds[i] == tmGrpIdx )	  // this is an ACTIVE block
		{
			beActive[i] = 1;
			blockSiz    = bb_infor + 8 * i + 4;
			dataBlks[i] = new float * [ nTmSteps ];
			for ( j = 0; j < nTmSteps; j ++ )
			dataBlks[i][j] = new float 
			[ blockSiz[0] * blockSiz[1] * blockSiz[2] * 3 ];
		}
	}
	blockSiz = NULL;


	// let process #0 gather the status of each spatial-temporal block from EACH process
	//
	// this information is used to guide the send-receive correspondance between process
	// #0 and each slave process
	//
	MPI_Gather( beActive, maxBlcks, MPI_INT, 
		    BEACTIVE, maxBlcks, MPI_INT, 0, MPI_COMM_WORLD );

	
	// directly fill in the blocks for process #0 only, without inter-process send-receive
	//
	// note: this part is separated from the allocation of "2nd-level (2D) and 3rd-level (1D)
	//	 arrays for ACTIVE blocks" part (see above) to save the use of if-statements
	//
	if ( thisProc == 0 )
	{
		for ( i = 0; i < nblocks; i ++ )
		{
			if ( beActive[i] == 0 ) continue;	// skip INactive blocks
		
			for ( j = 0; j < nTmSteps; j ++ )
			{
				FillBlock(  volSizes,  vec_data[j],  bb_infor + 8 * i,
				            bb_infor + 8 * i + 4,    dataBlks[i][j]  );
			}
		}
	}

	MPI_Barrier( MPI_COMM_WORLD );

	
	// the following code is currently USELESS, but remains here to demnstrate the use of 
	// MPI_Exscan in case that in the future 'maxBlcks' might be replaced with 'nblocks' 
	//
	// if ALL spatial-temporal blocks (of ALL processes) are numbered (and placed in a single
	// list) in a process-wise manner, blockId0 refers to the "global" Id of the first block 
	// of THIS process in such a list --- unnecessarily equal to the standard (or DIY) GID
	//
	// int	blockId0 = 0;
	// MPI_Exscan( &nblocks, &blockId0, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  	

	// fill in the data blocks for the slave processes
	//
	// msg_tags, sendReqs, and sendBufs are all excessively allocated  just to
	// make things simple, whereas sendBufs  does NOT waste memory much  since
	// each of its entries (initialized with NULL) is allocated ONLY ON DEMAND
	//
	// msg_tags and sendBufs both involve exact send-receive correspondance
	// sendReqs does NOT, but it seems unsafe to estimate the number of send-cooms
	//
	int           commFlag = 0; // whether a communication is completed (1) or not (0)
	int           flag_sum = 0; // number of the requests that have been COMPLETED
	int	      num_reqs = 0; // number of the requests that have been ISSUED
	int	    * msg_tags = new int  [ numProcs * maxBlcks * max_tblk ]; // comm IDs
	for ( i = 0;  i < numProcs * maxBlcks * max_tblk; i ++ )  msg_tags[i] = i;
	float      ** sendBufs = NULL;
	MPI_Request * sendReqs = NULL;
	MPI_Request * recvReqs = NULL;
	MPI_Status    comStats;	    // status of a request (PENDING or COMPLETED)
	
	if ( thisProc == 0 )
	{
		// sending blocks ......

		sendReqs = new MPI_Request [ numProcs * maxBlcks * max_tblk ];
		sendBufs = new float     * [ numProcs * maxBlcks * max_tblk ];	
		for ( i = 0;  i < numProcs * maxBlcks * max_tblk;  i ++ ) sendBufs[i] = NULL;

		for ( i = 1;  i < numProcs; i ++ ) // beginning with process #1
        	{
			int   blckBase = maxBlcks * max_tblk * i;   
			int * actvFlag = BEACTIVE + maxBlcks * i;
			int * infoBase = BB_INFOR + maxBlcks * 8 * i;
			for ( j = 0; j < numBlcks[i]; j ++, 
			      actvFlag ++, infoBase += 8, blckBase += max_tblk )
			{
				if (  ( *actvFlag )  ==  0  )  continue;
			
				for ( k = 0; k < nTmSteps; k ++ ) 
				{
					sendBufs[ blckBase + k ] 
					= new float [ infoBase[4] * infoBase[5] * infoBase[6] * 3 ];
					FillBlock(  volSizes,       vec_data[k],  infoBase, 
						    infoBase + 4,   sendBufs[ blckBase + k ]  );
					MPI_Isend(  sendBufs[ blckBase + k ],
					   	    infoBase[4] *   infoBase[5] * infoBase[6] * 3,
					   	    MPI_FLOAT,  i,  msg_tags[ blckBase + k ],
					   	    MPI_COMM_WORLD, sendReqs + num_reqs  );
					num_reqs ++;
				}
			}
			actvFlag = NULL;
			infoBase = NULL;
		}

		// make sure all of the 'num_reqs' data blocks are sent
		// out from the master process to the slave processes
		while ( flag_sum < num_reqs )
		{
			flag_sum = 0;
			for ( i = 0; i < num_reqs; i ++ )
			{
				MPI_Test( sendReqs + i, &commFlag, &comStats );
				flag_sum += commFlag ? 1 : 0;
			}
		}

		// the while-loop might be substituted by the following line
		// for ( i = 0; i < num_reqs; i ++ ) MPI_Wait( sendReqs + i, &comStats );
	}
	else
	{
		// receiving blocks ......
	
		recvReqs = new MPI_Request [ nblocks * max_tblk ]; // to make things simple

		for ( j = 0; j < nblocks; j ++ )
		{
			if ( beActive[j] == 0 ) continue; // INactive and will NOT receive

			int   
			msg_base = thisProc * maxBlcks * max_tblk + max_tblk * j;
			blockSiz = bb_infor + 8 * j + 4; // the same for the time steps below
			for ( k = 0; k < nTmSteps; k ++ ) 
			{	
				MPI_Irecv(  dataBlks[j][k],  
					    blockSiz[0] * blockSiz[1] * blockSiz[2] * 3,
				   	    MPI_FLOAT,  0,   msg_tags[ msg_base + k ],
				   	    MPI_COMM_WORLD,  recvReqs + num_reqs  );
				num_reqs ++;
			}
			blockSiz = NULL;
		}

		// make sure 'num_reqs'data blocks are all received by THIS process
		while ( flag_sum < num_reqs )
		{
			flag_sum = 0;
			for ( i = 0; i < num_reqs; i ++ )
			{
				MPI_Test( recvReqs + i, &commFlag, &comStats );
				flag_sum += commFlag ? 1 : 0;
			}
		}

		// the while-loop might be substituted by the following line
		// for ( i = 0; i < num_reqs; i ++ ) MPI_Wait( recvReqs + i, &comStats );
	}
	

	// now data blocks have been ready on every process (specifically in dataBlks)
	MPI_Barrier( MPI_COMM_WORLD ); 

	
	// clean up
	if ( thisProc == 0 )
	{
		for ( i = 0; i < numProcs * maxBlcks * max_tblk; i ++ )
		{  if ( sendBufs[i] ) delete [] sendBufs[i];  sendBufs[i] = NULL;  }

		delete [] sendBufs;  sendBufs = NULL;
		delete [] sendReqs;  sendReqs = NULL;
  	}
	else	// for each slave process
	{
		delete [] recvReqs;  recvReqs = NULL;
	}

	delete [] msg_tags;  	msg_tags = NULL;

	if ( thisProc == 0 ) 
	cout << ".... DIY data blocks distributed for time group #" 
             << tmGrpIdx << "." << endl;
}


// -----------------------------------------------------------------------
// 
// import seeds from a file (for detecting & fixing a bug with CurvilinearGrid.C)
// (bug: nearly every particle is integrated only one step)
// 
//
// added by Zhanping Liu on 06/14/2013
//
int	InitTraces4CurvilinearBlocks( char * seedFile )
{
	int	 thisProc;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );

	// load a set of seeds from the file
	int	 numSeeds = 0;
	float *  theSeeds = NULL;
	FILE  *  thisFile = fopen( seedFile, "rb" );
	fread(  &numSeeds,  sizeof(  int  ),  1,             thisFile  );
	theSeeds = new float [ numSeeds * 3 ];
	fread (  theSeeds,  sizeof( float ),  numSeeds * 3,  thisFile  );	
	fclose(  thisFile  );  thisFile = NULL;

	// assign (part of the) seeds to the blocks
	int	 totSeeds = 0;
	Particle particle;
	for ( int i = 0; i < nblocks; i ++ )
	{
		bb_t 	 boundBox;
		DIY_No_ghost_block_bounds( 0, i, &boundBox );

		Seeds[i].clear();
		float *  thisSeed = theSeeds;
		for ( int j = 0; j < numSeeds; j ++, thisSeed += 3 )
		{
			if ( thisSeed[0] >= boundBox.min[0] &&
			     thisSeed[0] <  boundBox.max[0] &&
			     thisSeed[1] >= boundBox.min[1] &&
			     thisSeed[1] <  boundBox.max[1] &&
			     thisSeed[2] >= boundBox.min[2] &&
			     thisSeed[2] <  boundBox.max[2] )
			{
				particle.pt.Set( thisSeed[0], thisSeed[1], 
					   	 thisSeed[2], 0 );
				particle.steps = 0;
				Seeds[i].push_back( particle );	
			}

			if (  int( Seeds[i].size() )  >=  tf  )  break;
		}
		thisSeed  = NULL;
		totSeeds += int( Seeds[i].size() );
	}

	delete [] theSeeds;	theSeeds = NULL;

	return totSeeds;
}


// -----------------------------------------------------------------------
// 
// generate a set of seeds based on the grid points of each block of a
// curvilinear grid
//
// added by Zhanping Liu on 06/14/2013
//
void	InitTraces4CurvilinearBlocks()
{
	int	  totSeeds = 0;
	Particle  particle;

	for ( int i = 0; i < nblocks; i ++ )
	{	
		int     spaceIdx = sBlckIds[i];
		int	pntIndex = 0;
		int   * blockInf = bb_infor + 8 * i;
		int	sliceSiz = blockInf[4] * blockInf[5];
		float * this_pnt = NULL;
		
		Seeds[i].clear();
		for ( int t = 2; t < blockInf[6] - 2; t += 5 )
		for ( int s = 2; s < blockInf[5] - 2; s += 5 )
		for ( int r = 2; r < blockInf[4] - 2; r += 5 )
		{
			pntIndex = t * sliceSiz + s * blockInf[4] + r;
			this_pnt = gridBlks[spaceIdx] + pntIndex * 3;
			particle.pt.Set( this_pnt[0], this_pnt[1], this_pnt[2], 0 );
			particle.steps = 0;
			Seeds[i].push_back( particle );	
			//if (  int( Seeds[i].size() )  >=  tf  )  break;
		}
		blockInf = NULL;
		this_pnt = NULL;

		totSeeds += int( Seeds[i].size() );
	}

	cout << endl << endl << "number of seeds: " << totSeeds << endl;
}


// -----------------------------------------------------------------------
//
// add a flow line to the constituents of a vtkPolyData (process #0)
//
// possibly add some flow line / cell stributes: process Id, trace Id, trace length
//
// to turn OFF any option: NULL for the pointer (procIdxs / strmIdxs / strmLens)
//                     and -1   for the value   (procIndx / strmIndx)
//
// added by Zhanping Liu on 05/30/2013
//
void	AddFlowLine2Polydata
	( int         numSamps, VECTOR4      * strmSmps,
          vtkPoints * vtk_pnts, vtkCellArray * vtkLines,
	  int         procIndx, vtkIntArray  * procIdxs,
          int         strmIndx, vtkIntArray  * strmIdxs, vtkIntArray * strmLens )
{
	VECTOR4   * thisSamp = strmSmps;
	vtkIdType * pointIds = new vtkIdType [ numSamps ];
	vtkIdType * thisPtId = pointIds;


	// insert point coordinates and their Ids
	for ( int i = 0; i < numSamps; i ++, thisSamp ++, thisPtId ++ )
	{
		*thisPtId = vtk_pnts->InsertNextPoint(  ( float * ) thisSamp  );
	}
	thisSamp = NULL;
	thisPtId = NULL;


	// insert the flow line as a cell --- the inter-point connectivity
	vtkIdType vtkLinId = vtkLines->InsertNextCell( numSamps, pointIds );
	delete [] pointIds;  pointIds = NULL;


	// options: insert the flow line / cell attributes
	if ( procIdxs ) procIdxs->InsertValue( vtkLinId, procIndx );
	if ( strmIdxs ) strmIdxs->InsertValue( vtkLinId, strmIndx );
	if ( strmLens ) strmLens->InsertValue( vtkLinId, numSamps );
}


// -----------------------------------------------------------------------
//
// place flow lines in a vtkPolyData (process #0)
//
// pIdxInfo: to save process   Ids    (1) or not (0)
// sIdxInfo: to save flow line Ids    (1) or not (0)
// sLenInfo: to save flow line length (1) or not (0)
// 
// actvScal: 0 ("ProcessIndex") / 1 ("FlowLineIndex") / 2 ("NumberOfSamples")
//
// added by Zhanping Liu on 05/30/2013
//
void	FlowLines2PolyData
	( vtkPolyData * polyData, int   pIdxInfo, 
	  int           sIdxInfo, int   sLenInfo,
	  int           numProcs, int * numStrms, int actvScal )
{
	// geometry and topology --- points and lines
	vtkPoints     * vtk_pnts = vtkPoints::New();
	vtkCellArray  * vtkLines = vtkCellArray::New();
	vtk_pnts->SetDataTypeToFloat(   );
	vtk_pnts->Allocate( 10000, 5000 );
	vtkLines->Allocate( 1000,  500  );


	// flow line attributes --- 3 cell data arrays
	// process Id, flow line Id, flow line length
	//
	vtkIntArray * pIdArray = NULL;
	vtkIntArray * sIdArray = NULL;
	vtkIntArray * lenArray = NULL;
	

	// process Id
	if ( pIdxInfo && numProcs && numStrms ) // numStrms: number of flow lines from each process 
	{
		pIdArray = vtkIntArray::New();
		pIdArray->Allocate( 1000, 500 );
		pIdArray->SetName( "ProcessIndex" );
	}

	// flow line Id
	if ( sIdxInfo )
	{
		sIdArray = vtkIntArray::New();
		sIdArray->Allocate( 1000,  500 );
		sIdArray->SetName( "FlowLineIndex" );
	}

	// flow line length (in number of samples)
	if ( sLenInfo )
	{
		lenArray = vtkIntArray::New();
		lenArray->Allocate( 1000, 500 );
		lenArray->SetName( "NumberOfSamples" );
	}


	// add each flow line to the vtkPolyData object
	int       i, j;	
	int	  strmIndx = 0;
	int     * numSamps = npt;
	VECTOR4 * thisStrm = pt;
	if ( pIdArray )
	{	
		#ifdef   _TEST_VTKIO_RESULT_ // verify the data distribution version 
		int      allStrms = 0;
		for ( i = 0; i < numProcs; i ++ ) allStrms += numStrms[i];
		FILE  *  out_file = fopen( "flowlines_WITH_proc.dat", "wb" );
		fwrite( &allStrms, sizeof( int ), 1, out_file );
		#endif
		

		for ( i = 0; i < numProcs;    i ++ )
		for ( j = 0; j < numStrms[i]; j ++, strmIndx ++, numSamps ++ )
		{	
			#ifdef _TEST_VTKIO_RESULT_
			fwrite( numSamps, sizeof( int     ), 1,         out_file );
			fwrite( thisStrm, sizeof( VECTOR4 ), *numSamps, out_file );
			#endif	
			
			AddFlowLine2Polydata
			( *numSamps, thisStrm, 
			   vtk_pnts, vtkLines,
	  		   i,        pIdArray,
          		   strmIndx, sIdArray, lenArray );

			thisStrm += ( *numSamps );	// next flow line
		}

		#ifdef _TEST_VTKIO_RESULT_
		fclose( out_file );
		#endif
	}
	else
	{	
		#ifdef   _TEST_VTKIO_RESULT_ // verify the data distribution version
		FILE  *  out_file = fopen( "flowlines_NONE_proc.dat", "wb" );
		fwrite( &tot_ntrace, sizeof( int ), 1, out_file );
		#endif

		for ( i = 0; i < tot_ntrace; i ++, numSamps ++ )
		{	
			#ifdef _TEST_VTKIO_RESULT_
			fwrite( numSamps, sizeof( int     ), 1,         out_file );
			fwrite( thisStrm, sizeof( VECTOR4 ), *numSamps, out_file );
			#endif

			AddFlowLine2Polydata
			( *numSamps, thisStrm, 
			   vtk_pnts, vtkLines,
	  		   -1,       NULL,
          		   i,        sIdArray, lenArray );

			thisStrm += ( *numSamps );	// next flow line
		}

		#ifdef _TEST_VTKIO_RESULT_
		fclose( out_file );
		#endif
	}
	numSamps = NULL;
	thisStrm = NULL;


	// attach the components to the vtkPolyData
	if ( actvScal < 0 || actvScal > 2 ) actvScal = 0;
	polyData->SetPoints( vtk_pnts );
	polyData->SetLines ( vtkLines );
	if ( pIdArray )
	{
		polyData->GetCellData()->AddArray( pIdArray ); 
		if ( actvScal == 0 ) 
		polyData->GetCellData()->SetActiveScalars( "ProcessIndex"    );
	}
	if ( sIdArray )
	{
		polyData->GetCellData()->AddArray( sIdArray );
		if ( actvScal == 1 )
		polyData->GetCellData()->SetActiveScalars( "FlowLineIndex" );
	}
	if ( lenArray ) 
	{
		polyData->GetCellData()->AddArray( lenArray );
		if ( actvScal == 2 )
		polyData->GetCellData()->SetActiveScalars( "NumberOfSamples" );
	}
	polyData->Squeeze();


	// release memory
	vtk_pnts->Delete();			vtk_pnts = NULL;
	vtkLines->Delete();			vtkLines = NULL;
	if ( pIdArray ) pIdArray->Delete();	pIdArray = NULL;
	if ( sIdArray ) sIdArray->Delete();	sIdArray = NULL;
	if ( lenArray ) lenArray->Delete();	lenArray = NULL;
}


// -----------------------------------------------------------------------
//
// place flow lines in a vtkPolyData and possibly save it to a disk file
//
// pIdxInfo: to save process   Ids    (1) or not (0)
// sIdxInfo: to save flow line Ids    (1) or not (0)
// sLenInfo: to save flow line length (1) or not (0)
//
// actvScal: 0 ("ProcessIndex") / 1 ("FlowLineIndex") / 2 ("NumberOfSamples")
//
// fileType: VTK_ASCII or VTK_BINARY
//
// added by Zhanping Liu on 05/30/2013
//
void	FlowLines2PolyData( int pIdxInfo, int    sIdxInfo, int sLenInfo,
	                    int actvScal, char * vtk_file, int fileType )
{
	int   thisProc, numProcs;
  	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
  	MPI_Comm_size( MPI_COMM_WORLD, &numProcs );


	if ( thisProc == 0 )
	cout << endl << "==== accessing and converting flow lines to VTK data ====" << endl;


	// let process #0 know the number of traces from each process
	// 
	// changes made to ParFlow.h/C
	//
	int * numStrms = NULL;
	int   n_traces = parflow->GetNumberOfFlowLines(); // for a single process
	if ( thisProc == 0 ) numStrms = new int [ numProcs ];
	
	//MPI_Barrier( MPI_COMM_WORLD ); // unnecessary
	MPI_Gather( &n_traces, 1, MPI_INT, numStrms, 1, MPI_INT, 0, MPI_COMM_WORLD );
	cout << endl << "number of traces from process #" 
	             << thisProc << ": " << n_traces << endl;

	// let process #0 convert flow lines to a vtkPolyData object
	if ( thisProc == 0 )
	{	
		polyData = vtkPolyData::New();	// a global object

		// fill in the vtkPolyData
		FlowLines2PolyData( polyData, pIdxInfo, sIdxInfo, sLenInfo,
	          		    numProcs, numStrms, actvScal );
		delete [] numStrms;   numStrms = NULL;
		
		cout << endl << ".... VTK data generated." << endl;

		// write the vtkPolyData to a VTK file
		if ( vtk_file )
		{
			vtkPolyDataWriter  *   pdWriter = vtkPolyDataWriter::New();
			pdWriter->SetFileName( vtk_file );
			pdWriter->SetFileType( fileType );
			pdWriter->SetInput   ( polyData );
			pdWriter->Update();
			pdWriter->Delete();    pdWriter = NULL;

			cout << endl << ".... flow lines also saved to " << vtk_file << "." << endl;
		}

		cout << endl << "  there are " << polyData->GetNumberOfCells() <<  " flow lines and " 
			     << polyData->GetNumberOfPoints() << " points." << endl;
	}

	//MPI_Barrier( MPI_COMM_WORLD ); // unnecessary
}


// -----------------------------------------------------------------------
//
// render flow lines (in the form of a vtkPolyData object) that are generated
// by processes in parallel but eventually collected and maintained by the
// master process
//
// rendered by process #0 only
//
// added by Zhanping Liu on 05/29/2013
//
void	RenderFlowLines( int be_tubed, float tube_rad, int tube_res )
{
	int	thisProc;
	MPI_Comm_rank( MPI_COMM_WORLD, &thisProc );
	if ( thisProc ) return;


	cout << endl << "==== rendering flow lines using VTK ====" << endl << endl;
   

	// map (the constituents of) the vtkPolyData object to 
	// renderable / graphical primitives
	vtkTubeFilter     * strmTube = NULL;
	vtkPolyDataMapper * strmMapr = vtkPolyDataMapper::New();
	if ( be_tubed )
	{	// tubes
		strmTube = vtkTubeFilter::New();
		strmTube->SetRadius( tube_rad );
		strmTube->SetNumberOfSides( tube_res );
		strmTube->SetInput( polyData );
		strmTube->Update();
  		strmMapr->SetInput( strmTube->GetOutput() );
	}
	else	strmMapr->SetInput( polyData );	// lines
	strmMapr->ImmediateModeRenderingOn();
	strmMapr->SetScalarVisibility( 1 );
	strmMapr->SetScalarModeToUseCellData();
	strmMapr->SetScalarRange( polyData->GetScalarRange() );
  	

    	// the vtkPolyData object is an actor to show up in the 3D scene
    	vtkActor * strmActr = vtkActor::New();
    	strmActr->SetMapper( strmMapr );


    	// a renderer takes up (part or all of) a window on the screen
    	// to draw the constituent actors
    	vtkRenderer * renderer = vtkRenderer::New();
    	renderer->AddActor( strmActr );
    	renderer->SetBackground( 0.1, 0.2, 0.4 );


    	// create a render window to for the renderer
    	vtkRenderWindow * rendrWin = vtkRenderWindow::New();
	rendrWin->SetWindowName( "Argonne National Lab:    VTK ---> [DIY + OSUFlow] ---> VTK  ----  towards Parallel Visualization of Large Flow Data" );
    	rendrWin->AddRenderer( renderer );
    	rendrWin->SetSize( 1000, 1000 );


	// create an interactor for the user to play with the window (actually the renderer)
	#ifdef _USE_SELF_ROTATION_
    	while ( 1 ) 	    // let the flow lines rotate by themselves
   	{
      		rendrWin->Render();
      		renderer->GetActiveCamera()->Azimuth( 0.5 );
		renderer->GetActiveCamera()->Roll( 0.5 );
    	}
	#else		    // allow the user to manipulate the flow lines
	vtkRenderWindowInteractor * intactor = vtkRenderWindowInteractor::New();
	intactor->SetRenderWindow( rendrWin );
	intactor->Start();
	rendrWin->Render(); // this line MUST be AFTER intactor->SetRenderWindow( rendrWin )
	#endif
  

	// memory deallocation
	// do NOT destroy polyData here!!! let's destroy it in a more obvious place 
  	#ifndef _USE_SELF_ROTATION_
	intactor->Delete();
	#endif
    	rendrWin->Delete();
	renderer->Delete();
	strmActr->Delete();
	strmMapr->Delete();
	if ( strmTube )
	strmTube->Delete();
}


// -----------------------------------------------------------------------
//
// release memory used for data loading and distribution
//
// added by Zhanping Liu on 05/28/2013 and updated on 06/04/2013
//
void 	ReleaseMemory()
{
	cout << "==== releasing memory (for DIY blocks distribution and VTK rendering) ====" << endl;


	#ifdef _CURVILINEAR_GRIDS_	
	if ( sBlckIds ) delete [] sBlckIds;    	sBlckIds = NULL;	    // for ALL processes
	if ( stBlkIds ) delete [] stBlkIds;	stBlkIds = NULL;	    // fpr ALL processes
	if ( STBLKIDS ) delete [] STBLKIDS;	STBLKIDS = NULL;	    // fpr process #0 only
	if ( grd_cord ) delete [] grd_cord;    	grd_cord = NULL;	    // for process #0 only

	// destroy the curvilinear grid BLOCK (sub-volume) data (coordinates): for ALL processes
	if ( gridBlks ) 
	{
		for ( int i = 0; i < nsBlocks; i ++ )
		{
			if ( gridBlks[i] ) delete [] gridBlks[i];  gridBlks[i] = NULL;
		}

		delete [] gridBlks;	       	gridBlks = NULL;
	}
	#endif


	// destroy the global input data
	if ( vec_data ) 
	{
		for ( int i = 0; i < max_tblk; i ++ )		  // for process #0 only
		{
			if ( vec_data[i] ) delete [] vec_data[i]; 
			vec_data[i] = NULL;
		}

		delete [] vec_data;	       	vec_data = NULL;  // for ALL processes
	}

	
	if ( numBlcks ) delete [] numBlcks;     numBlcks = NULL;  // for process #0 only
	if ( BEACTIVE ) delete [] BEACTIVE;	BEACTIVE = NULL;  // for process #0 only
	if ( BB_INFOR ) delete [] BB_INFOR;     BB_INFOR = NULL;  // for process #0 only
	if ( tmGrpIds ) delete [] tmGrpIds;	tmGrpIds = NULL;  // for ALL processes
	if ( beActive ) delete [] beActive;	beActive = NULL;  // for ALL processes
	if ( bb_infor ) delete [] bb_infor;     bb_infor = NULL;  // for ALL processes
  	if ( tBlkMins ) delete [] tBlkMins;	tBlkMins = NULL;  // for ALL processes
	if ( tBlkSizs ) delete [] tBlkSizs;	tBlkSizs = NULL;  // for ALL processes


	// for the global output (rendering) data
	if ( polyData ) polyData->Delete();	polyData = NULL;  // for process #0 only


	// NOTE: Here we must NOT deallocate the 2nd-level and 3rd-level arrays of dataBlks,
	// of which the 2nd-level arrays have been delivered and attached to block_osuflow[i]
	// in Blocks.C through LoadBlocks4D( ...... ) in the form of the last argument and
	// hence these 2nd-level arrays (and the 3rd-level arrays) will be released by
	// block_osuflow[i] as block_osuflow[i] is destroyed.
        //
	if ( dataBlks ) delete [] dataBlks;	dataBlks = NULL;  // for ALL processes	

	cout << endl << ".... memory released." << endl << endl;
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
