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
#include "bil.h"
#include "assignment.hpp"
#include "blocking.hpp"
#include "neighborhoods.hpp"

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

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void Cleanup();
void Run();
void Header(char *filename, float *size, int *tsize, float *vec_scale);
void LoadSeedsFromFile();

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
int end_steps; // final number of points each particle should travel
DataMode data_mode; // data format
Blocks *blocks; // block class object
ParFlow *parflow; // parallel flow class object
Blocking *blocking; // blocking class object
Assignment *assign; // assignment class object
Neighborhoods *nbhds; // neighborhoods class object
int compute_begin, compute_end; // jumpshot states
int rank; // mpi rank
float vec_scale; // vector scaling factor
char seed_file[256]; // seed file name
int seed_file_num = 0; // number of seeds read from the seed file
VECTOR3* seed_file_seeds = NULL; // seeds read from the seed file
float wf = 0.1; // wait factor for nonblocking communication
                // wait for this portion of messages to arrive each round
const int ghost = 1;  // a ghost layer of at least 1 is required for correct
		      // advection between blocks

// integration parameters
const float maxError = 0.001;
const float initialStepSize = 1.0;
const float minStepSize = 0.01;
const float maxStepSize = 5.0;
const float lowerAngleAccuracy = 3.0;
const float upperAngleAccuracy = 15.0;

const INTEG_ORD integrationOrder = RK45;
const bool useAdaptiveStepSize = true;

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
#ifdef USE_BIL
  BIL_Init(MPI_COMM_WORLD);
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  Init();

  // run
  MPI_Barrier(MPI_COMM_WORLD);
  TotTime = MPI_Wtime();
  Run();
  MPI_Barrier(MPI_COMM_WORLD);
  TotTime = MPI_Wtime() - TotTime;

  // print the performance stats
  parflow->PrintPerf(TotTime, TotInTime, TotOutTime,
		     TotCompCommTime, TotParticles, size);

#ifdef GRAPHICS

  if (rank == 0) {
    VECTOR3 min, max;
    min = VECTOR3(0.0f, 0.0f, 0.0f);
    max = VECTOR3((float)(size[0] - 1), (float)(size[1] - 1),
		  (float)(size[2] - 1));
    DrawInit(pt, npt, tot_ntrace, argc, argv, min, max, 0);
  }

#endif

  Cleanup();
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef USE_BIL
  BIL_Finalize();
#endif
  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// Run
//
void Run() {

  int i, j;
  int g; // current group
  double time; // time to load a block
  double t0;

#ifdef REPARTITION
  assert(nblocks <= MAX_BLK);
  for (i = 0; i < nblocks; i++)
    wgts[i] = 0;
#endif

  parflow->InitTraces(Seeds, tf, nblocks, tsize, ntpart, seed_file_seeds, 
		      seed_file_num);

#ifdef USE_BIL
  float ***bil_data = NULL;
#endif

  // for all groups
  for (g = 0; g < ntpart; g++) {

    // synchronize before starting I/O
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();

    // delete blocks from previous time group
    if (g > 0) { 
      blocks->DeleteBlocks(g, tsize, ntpart, nblocks);
    }

    // todo: change seeds to vector in repartition
    // #ifdef REPARTITION
    //     parflow->Repartition(g, &nblocks, &Seeds, 0, &osuflow,
    // 		       block, OSUFLOW, MPI_COMM_WORLD, wgts);
    // #endif

    // load blocks for this time group
#ifdef USE_BIL
    bil_data = blocks->BilLoadTimeGroupBlocks(g, nblocks, size, tsize, ntpart);
    blocks->LoadBlocks4D(g, &time, nblocks, size, tsize, ntpart, bil_data);
    delete[] bil_data;
#else
    blocks->LoadBlocks4D(g, &time, nblocks, size, tsize, ntpart);
#endif

    // synchronize after I/O
    MPI_Barrier(MPI_COMM_WORLD);
    TotInTime += (MPI_Wtime() - t0);
    t0 = MPI_Wtime();

    // scale blocks to improve visibility
    for (i = 0; i < nblocks; i++) {
      if (blocks->GetLoad(i))
	osuflow[i]->ScaleField(vec_scale);
    }

#ifdef REPARTITION
    assert(nblocks <= MAX_BLK);
    for (i = 0; i < nblocks; i++)
      wgts[i] = 0;
#endif

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

      // for all blocks
      for (i = 0; i < nblocks; i++) {

#ifdef MPE
	MPE_Log_event(compute_begin, 0, NULL);
#endif

	// compute fieldlines
	if (tsize > 1) {
#ifdef REPARTITION
	  assert(i < MAX_BLK);
	  parflow->ComputePathlines(Seeds[i], i, pf, end_steps, &wgts[i]);
#else
	  parflow->ComputePathlines(Seeds[i], i, pf, end_steps);
#endif
	}
	else {
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
      parflow->FlushNeighbors(Seeds);

    } // for all rounds

    // flush any remaining messages
    //parflow->FlushNeighbors(Seeds);

#ifdef REPARTITION
    AdvanceWeights(g);
#endif

    // end time group synchronized to get accurate timing
    MPI_Barrier(MPI_COMM_WORLD);
    TotCompCommTime += (MPI_Wtime() - t0);

  } // for all groups

  // synchronize prior to gathering
  MPI_Barrier(MPI_COMM_WORLD);
  TotOutTime = MPI_Wtime();

  // gather fieldlines for rendering
  parflow->GatherFieldlines(nblocks, size, tsize);
  MPI_Barrier(MPI_COMM_WORLD);
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
  strncpy(filename, argv[1], sizeof(filename));
  Header(filename, size, &tsize, &vec_scale);
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

  assert(nspart * ntpart >= nproc);
  assert(ntpart <= tsize);
  if(seed_file[0] == '!')
    assert(tf > 0);

  // partition domain
  // todo: don't create partition if there is a part file?
  int64_t data_size[4] = {size[0], size[1], size[2], tsize};
  int maxb; // unused
  int64_t given[4] = {0, 0, 0, ntpart};
  assign = new Assignment(nspart * ntpart, nblocks, maxb, MPI_COMM_WORLD);

  // which dimensions to apply ghost cells in
  int ghost_dim[4] = {1, 1, 1, 1};

  // which directions to apply ghost cells, for each dimension, 
  // minimum sides only (-1), maximum sides only (1), 
  // or all sides equally (0)
  int ghost_dir[4] = {0, 0, 0, 0};

  // if computing streamlines, then the time dimension needs to be restricted
  // so that a ghost cell layer is not added on. if computing pathlines, only
  // add a ghost layer to the maximum side in the time dimension
  if(tsize == 1) {
    ghost_dim[3] = 0;
  }
  else {
    ghost_dir[3] = 1;
  }

  blocking = new Blocking(4, nspart * ntpart, data_size, 1, ghost, ghost_dir,
			  ghost_dim, given, assign, MPI_COMM_WORLD); 

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
  blocks = new Blocks(blocking, assign, (void *)osuflow, OSUFLOW, 
		      dataset_files, num_dataset_files, data_mode, ghost);
  parflow = new ParFlow(blocking, assign, blocks, osuflow, sl_list, 
			&pt, &npt, &tot_ntrace, nblocks, 0);
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
  delete [] sl_list;

  for (i = 0; i < nblocks; i++) {
    if (osuflow[i] != NULL)
      delete osuflow[i];
  }
  for (i = 0; i < Seeds.size(); i++)
    Seeds[i].clear();
  Seeds.clear();

  delete blocks;
  delete parflow;
  delete blocking;
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

  assert((fp = fopen(filename, "r")) != NULL);

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
