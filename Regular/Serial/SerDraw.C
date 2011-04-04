//------------------------------------------------------------------------------
//
// serial test draw
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
#include "Lattice4D.h"
#include "Draw.h"
#include "ParFlow.h"

// drawing data
VECTOR4 *pt; // points in everyone's traces
int *npt; // everyone's number of points in their traces
int tot_ntrace; // total number of everyone's traces

// function prototypes
void GetArgs(int argc, char *argv[]);
void Init();
void Cleanup();
void Run();
void Header(char *filename, float *size, int *tsize, float *vec_scale);

// globals
static char filename[256]; // dataset file name
static char **dataset_files = NULL; // files in the dataset
int num_dataset_files = 0;
float size[3]; // spatial domain size
static int tsize; // temporal domain size
vector< vector<Particle> > Seeds; // list of seeds lists
VECTOR3 *seeds; // one temporary list of (3d) seeds
OSUFlow **osuflow; // one flow object for each block
list<vtListTimeSeedTrace*> *sl_list; // pathlines list
int nspart; // global total number of spatial blocks
int ntpart; // global total number of temporal blocks
int nblocks; // my number of blocks
Lattice4D* lat; // lattice
int tf; // max number of traces per block
const int pf = 1000; // max number of points per trace in each round
const int max_rounds = 100; // max number of rounds
int end_steps; // final number of points each particle should travel
DataMode data_mode; // data format
Blocks *blocks; // blocks class object
ParFlow *parflow; // parallel flow class object
float vec_scale; // vector scaling factor

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

  VECTOR3 min, max;
  min = VECTOR3(0.0f, 0.0f, 0.0f);
  max = VECTOR3((float)(size[0] - 1), (float)(size[1] - 1),
		(float)(size[2] - 1));
  DrawInit(pt, npt, tot_ntrace, argc, argv, min, max, 0);

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

  parflow->InitTraces(Seeds, tf, nblocks, tsize, ntpart);

  // for all groups
  for (g = 0; g < ntpart; g++) {

    if (g > 0) // delete blocks from earlier groups
      blocks->DeleteBlocks(g, tsize, ntpart, nblocks);

    blocks->LoadBlocks4D(g, &time, nblocks, size, tsize, ntpart);

    // scale blocks to improve visibility
    for (i = 0; i < nblocks; i++) {
      if (lat->GetLoad(i))
        osuflow[i]->ScaleField(vec_scale);
    }

    // for all rounds
    for (j = 0; j < max_rounds; j++) {

      // for all blocks
      for (i = 0; i < nblocks; i++) {

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
  parflow->SerialGatherFieldlines(nblocks);

}
//-----------------------------------------------------------------------
//
// GetArgs
//
// gets command line args
//
void GetArgs(int argc, char *argv[]) {

  assert(argc >= 6);

  strncpy(filename,argv[1],sizeof(filename));
  Header(filename, size, &tsize, &vec_scale);
  nspart = atoi(argv[2]); // total space partitions
  ntpart = (tsize == 1 ? 1 : atoi(argv[3])); // total time partitions
  tf = atoi(argv[4]) / nspart; // traces per block
  switch(atoi(argv[5])) {
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
  assert(tf > 0);

  // init lattice and osuflow
  lat = new Lattice4D((int)size[0], (int)size[1], (int)size[2], tsize, 
		      nspart, &ntpart);

  nblocks = nspart * ntpart;
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
  delete [] osuflow;
  delete lat;

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
	dataset_files[n - 5] = strdup(token);
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
