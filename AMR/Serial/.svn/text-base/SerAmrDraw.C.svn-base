//------------------------------------------------------------------------------
//
// serial AMR test draw
//
// Han-Wei Shen
// The Ohio State University
// Columbus, OH
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

#include <stdio.h>
#include <stdlib.h> 
#include <list>
#include <iterator>
#include "OSUFlow.h"
#include "LatticeAMR.h"
#include "Draw.h"
#include "Blocks.h"
#include "ParFlow.h"

// drawing data
VECTOR4 *pt = NULL; // points in everyone's traces
int *npt = NULL; // everyone's number of points in their traces
int tot_ntrace; // total number of everyone's traces

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void Cleanup();
void Run();
void Header(char *filename, int *tsize, float *data_scale, float *vec_scale);

// globals
static char filename[256]; // dataset file name
static char **dataset_files = NULL; // files in the dataset
int num_dataset_files = 0;
int dims[3]; // number of cells in an AMR block (eg 16x16x16)
VECTOR3 size; // spatial domain size
static int tsize; // temporal domain size
vector< vector<Particle> > Seeds; // list of seeds lists
VECTOR3 *seeds; // one temporary list of (3d) seeds
OSUFlow **osuflow; // one flow object for each block
list<vtListTimeSeedTrace*> *sl_list = NULL; // pathlines list
int nspart; // global total number of spatial blocks
int ntpart; // global total number of temporal blocks
int nblocks; // my number of blocks
LatticeAMR* lat; // lattice
int tf; // max number of traces per block
const int pf = 1000; // max number of points per trace
const int max_rounds = 100; // max number of rounds
int end_steps; // final number of points each particle should travel
DataMode data_mode; // data format
float data_scale; // data scaling factor
float vec_scale; // vector scaling factor
char vx[256], vy[256], vz[256]; // names of velocity variables
Blocks *blocks; // block class object
ParFlow *parflow; // parallel class object

//----------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  char buf[256];

  GetArgs(argc, argv);
  Init();
  blocks = new Blocks(lat, (void *)osuflow, OSUFLOW, dataset_files, 
		    num_dataset_files, data_mode);
  parflow = new ParFlow(lat, osuflow, sl_list, &pt, &npt,
		    &tot_ntrace, nblocks);
  Run();

#ifdef GRAPHICS

  float fmin[4], fmax[4];
  lat->GetExtents(fmin, fmax);
  VECTOR3 min, max;
  min = VECTOR3(fmin[0], fmin[1], fmin[2]);
  max = VECTOR3(fmax[0], fmax[1], fmax[2]);
  DrawInit(pt, npt, tot_ntrace, argc, argv, min, max, 1);

#endif

  Cleanup();

}
//-----------------------------------------------------------------------
//
// Run
//
void Run() {

  int i, j, k;
  int g; // current group
  double time; // unused

  // init
  parflow->InitTraces(Seeds, tf, nblocks, tsize, ntpart);

  // for all groups
  for (g = 0; g < ntpart; g++) {

    if (g > 0) // delete blocks from earlier groups
      blocks->DeleteBlocks(g, tsize, ntpart, nblocks);
    blocks->LoadBlocksAMR(g, &time, data_mode); // load blocks for this time group

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

      // for all blocks
      for (i = 0; i < nblocks; i++) {

	if (lat->GetLoad(i))
	  osuflow[i]->ScaleField(vec_scale);

	// compute fieldlines
	if (tsize > 1)
	  parflow->ComputePathlines(Seeds[i], i, pf, end_steps);
	else
	  parflow->ComputeStreamlines(Seeds[i], i, pf, end_steps);

      } // for all blocks

      parflow->SerExchangeNeighbors(Seeds);

    } // for all rounds

  } // for all groups

  // gather fieldlines for rendering
  parflow->SerialGatherFieldlines(nblocks, size, tsize);

}
//-----------------------------------------------------------------------
//
// GetArgs
//
// gets command line args
//
void GetArgs(int argc, char *argv[]) {

  assert(argc >= 8);
  strncpy(filename,argv[1],sizeof(filename));
  Header(filename, &tsize, &data_scale, &vec_scale);
  size[0] = size[1] = size[2] = 1.0f;
  ntpart = (tsize == 1 ? 1 : atoi(argv[2])); // total time partitions
  tf = atoi(argv[3]); // total number of particles
  switch(atoi(argv[4])) {
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
  strncpy(vx, argv[5], sizeof(vx));
  strncpy(vy, argv[6], sizeof(vy));
  strncpy(vz, argv[7], sizeof(vz));

}
//-----------------------------------------------------------------------
//
// Init
//
// inits the app
//
void Init() {

  int nt; // max total number of traces
  int np; // max total number of points
  int i;

  assert(ntpart <= tsize);

  // create the lattice
  lat = new LatticeAMR(dataset_files, tsize, ntpart, vx, vy, vz,
		       data_mode, data_scale); 
  nblocks = lat->GetMyNumPartitions();
  nspart = lat->GetTotalNumPartitions(); // all spatial for now
  tf /= nspart; // number of traces per block
  assert(tf > 0);
  lat->GetBlockDims(dims);

  // init osuflow
  assert((osuflow = new OSUFlow*[nblocks]) != NULL);
  for (i = 0; i < nblocks; i++)
    assert((osuflow[i] = new OSUFlow) != NULL);

  // seeds and fieldline list
  Seeds.resize(nblocks);
  sl_list = new list<vtListTimeSeedTrace*>[nspart * ntpart];

}
//-----------------------------------------------------------------------
//
// Cleanup
//
// frees memory and such
//
void Cleanup() {

  int i;

  if (pt != NULL)
    delete [] pt;
  if (npt != NULL)
    delete [] npt;
  if (sl_list != NULL)
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
  delete [] osuflow;
  delete lat;

}
//-----------------------------------------------------------------------
//
// reads header file
//
// filename: filename (input)
// tsize: number of timesteps (output)
// data_scale: scaling factor of entire dataset (output)
// vec_scale: extra vector scaling factor (output)
//
void Header(char *filename, int *tsize, float *data_scale, float *vec_scale) {

  FILE *fp;
  char line[1024];
  char *token;
  const char delims[] = " \t\n";
  int n = 0; // number of tokens parsed so far

  assert((fp = fopen(filename, "r")) != NULL);

  while (fgets(line, sizeof(line), fp) != NULL) {

    token = strtok(line, delims);

    while (*line != '\0' && *line != '\n' && *token != '#') {

      if (n == 0) {
	*tsize = atoi(token);
	num_dataset_files = *tsize;
	assert((dataset_files = (char **)malloc(sizeof(char *) * 
						num_dataset_files)) != NULL);
	n++;
      }

      else if (n == 1) {
	*data_scale = atof(token);
	n++;
      }

      else if (n == 2) {
	*vec_scale = atof(token);
	n++;
      }

      else if (n >= 3) {
	if (n - 3 >= *tsize)
	  break;
	dataset_files[n - 3] = strdup(token);
	assert(dataset_files[n - 3] != NULL);
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
