//------------------------------------------------------------------------------
//
// mpi test main
//
// Copyright (c) 2008 Han-Wei Shen and Tom Peterka
//
// Contact:
//
// Han-Wei Shen
// The Ohio State University
// Columbus, OH
//
// Tom Peterka
// MCS Radix Lab
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// All rights reserved. May not be used, modified, or copied
// without permission
//
//------------------------------------------------------------------------------

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

#include "calc_subvolume.h"
#include "Lattice.h"

#define MAX_RECV_ATTEMPTS 10
#define MAX_ITERATIONS 10

// function prototypes
void Config(int argc, char *argv[]);
void PostPoint(VECTOR3 p, int neighbor);
void PrintPost();
void PrintRecv();
void SendNeighbors(Lattice* lat);
int ReceiveNeighbors(Lattice* lat);
void Reinit();

// globals
static char filename[256]; // datset file name
static VECTOR3 minLen, maxLen; // data bounds

// send and receive points
int NumSendPoints[6]; // number of points ready to send
int SizeSendPoints[6]; // size of sending points list (bytes)
float *SendPoints[6]; // sending points list
int NumRecvPoints[6]; // number of points received
int SizeRecvPoints[6]; // size of receiving points list (bytes)
float *RecvPoints[6]; // receiving points list
int NumSeeds; // number of seeds
int SizeSeeds; // size of seeds list (bytes)
VECTOR3 *Seeds; // seeds list

VECTOR3 size; // domain size

//---------------------------------------------------------------------------------

main(int argc, char *argv[]) {

  VECTOR3 minB, maxB; // subdomain bounds
  volume_bounds_type *vb_list;  // volume bounds list
  int nproc;  // mpi groupsize
  int rank; // mpi rank
  float from[3], to[3]; // seed points limits
  list<vtListSeedTrace*> list; // seed trace
  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  std::list<VECTOR3*>::iterator pnIter; 
  VECTOR3 p; // for debug purposes
  Lattice* lat; // lattice
  int lattice_xdim, lattice_ydim, lattice_zdim; // lattice bounds
  int neighbor; // neighbor's number (0-5)
  int ei, ej, ek; // neighbor's lattice position
  int i, j;
  int ghost; // number of ghost cells per side
  int init_seed_num = 10; // initial number of seeds

  if (rank == 0)
    printf("Initializing...\n"); 

  // initialize mpi and the app
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Barrier(MPI_COMM_WORLD);
  Config(argc, argv);
  OSUFlow *osuflow = new OSUFlow(); 
  ghost = 1;

  // subdivide domain and compute boundary
  lat = new Lattice(size[0], size[1], size[2], ghost, nproc);
  vb_list = lat->GetBoundsList(); 
  lat->InitSeedLists(); 
  minB.Set(vb_list[rank].xmin, vb_list[rank].ymin, vb_list[rank].zmin);
  maxB.Set(vb_list[rank].xmax, vb_list[rank].ymax, vb_list[rank].zmax);

  // read data
  osuflow->ReadData(filename, true, minB, maxB, size); 
  if (rank == 0)
    fprintf(stderr,"read file %s\n", filename); 

  // init seeds
  from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
  to[0]   = maxB[0];   to[1]   = maxB[1];   to[2]   = maxB[2]; 
  osuflow->SetRandomSeedPoints(from, to, init_seed_num); 
  Seeds = osuflow->GetSeeds(NumSeeds); 
  SizeSeeds = NumSeeds * sizeof(VECTOR3);

  // trace particles until they remain constant for some period of time

  for (i = 0; i < MAX_ITERATIONS; i++) {

    if (NumSeeds) {

      lat->ResetSeedLists();    // clear up the lattice seed lists

      // perform the integration
      // does this integrate in both directions?
      // if not, should it?
      osuflow->SetIntegrationParams(1, 5); 
      osuflow->GenStreamLines(Seeds, FORWARD_DIR, NumSeeds, 50, list); 

      // debug
      printf("rank %d completed integration i = %d\n", rank, i); 

      // redistribute boundary points to neighbors

      pIter = list.begin(); 

      // for all boundary points
      for (; pIter!=list.end(); pIter++) {

	trace = *pIter; 
	if (trace->size() == 0) continue; 
	pnIter = trace->end(); 
	pnIter--;
	p = **pnIter;

	// find which neighbor the point is in
	neighbor = lat->GetNeighbor(rank, p[0], p[1], p[2], ei, ej, ek); 

	// post the point to the send list
	if (neighbor != -1)
	  PostPoint(p, neighbor);

      } // for all boundary points

	// send boundary list to neighbors
      SendNeighbors(lat);

    } // if (num_seeds)

    // receive boundary list from neighbors
    for (j = 0; j < MAX_RECV_ATTEMPTS; j++) {
      if (ReceiveNeighbors(lat))
	break;
      usleep(1000000);
    }
    if (j == MAX_RECV_ATTEMPTS)
      Error("Rank %d giving up waiting for more work - program terminating\n",rank);

      // prepare for next iteration
      Reinit();

  } // for all iterations

  // debug: print current seeds
  fprintf(stderr,"\n");
  for (i = 0; i < NumSeeds; i++)
    fprintf(stderr,"Rank %d Seed %d: %.3f\t%.3f\t%.3f\n", rank, i, Seeds[i][0], Seeds[i][1], Seeds[i][2]);
  fprintf(stderr,"\n");

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}
//-----------------------------------------------------------------------
//
// Config
//
// gets command line args
//
void Config(int argc, char *argv[]) {

  int rank;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 5)
    Error("Error: Insufficient number of command line args\n");

  strncpy(filename,argv[1],sizeof(filename));

  // hard code the minimum corner to be 0,0,0
  // need to allow for variable data origin in the future
  minLen[0] = minLen[1] = minLen[2] = 0.0;

  maxLen[0] = atof(argv[2]) - 1.0f;
  maxLen[1] = atof(argv[3]) - 1.0f;
  maxLen[2] = atof(argv[4]) - 1.0f;

  size[0] = maxLen[0] - minLen[0] + 1.0f;
  size[1] = maxLen[1] - minLen[1] + 1.0f;
  size[2] = maxLen[2] - minLen[2] + 1.0f;

  if (rank == 0)
    printf("Volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
	   minLen[0], maxLen[0], minLen[1], maxLen[1], minLen[2], maxLen[2]); 

  // init sending and receiving points
  for (i = 0; i < 6; i++){

    NumSendPoints[i] = NumRecvPoints[i] = 0;
    if ((SendPoints[i] = (float *)malloc(3 * sizeof(float))) == NULL)
      Error("Error: Config() cannot allocate memory\n");
    if ((RecvPoints[i] = (float *)malloc(3 * sizeof(float))) == NULL)
      Error("Error: Config() cannot allocate memory\n");
    SizeSendPoints[i] = SizeRecvPoints[i] = 3 * sizeof(float);

  }

}
//------------------------------------------------------------------------
//
// PostPoint
//
// posts a point for sending to a neighbor
//
void PostPoint(VECTOR3 p, int neighbor) {


  while (SizeSendPoints[neighbor] < 
      (NumSendPoints[neighbor] + 1) * 3 * sizeof(float)) {

    SendPoints[neighbor] = (float *)realloc(SendPoints[neighbor],
      SizeSendPoints[neighbor] * 2);

    if (SendPoints[neighbor] == NULL)
      Error("Error: PostPoint() cannot reallocate memory\n");

    SizeSendPoints[neighbor] *= 2;

  }

  SendPoints[neighbor][NumSendPoints[neighbor] + 0] = p[0];
  SendPoints[neighbor][NumSendPoints[neighbor] + 1] = p[1];
  SendPoints[neighbor][NumSendPoints[neighbor] + 2] = p[2];

  NumSendPoints[neighbor]++;

}
//------------------------------------------------------------------------
//
// SendNeigbors()
//
// sends points to all neighbors
//
void SendNeighbors(Lattice* lat) {

  int myrank;
  int ranks[6];
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  lat->GetNeighborRanks(myrank, ranks);

  for (i = 0; i < 6; i++) {

    if (ranks[i] >= 0) {

      MPI_Send(&(NumSendPoints[i]), 1, MPI_INT, ranks[i], 0, MPI_COMM_WORLD);
      if (NumSendPoints[i])
	MPI_Send(SendPoints[i], NumSendPoints[i] * 3, MPI_FLOAT, ranks[i], 0 , MPI_COMM_WORLD);

    }

  }

}
//--------------------------------------------------------------------------
//
// ReceiveNeighbors
//
// receives points from all neighbors
//
// returns total number of points received
//
int ReceiveNeighbors(Lattice* lat) {

  MPI_Status status;
  int myrank;
  int ranks[6];
  int num = 0;
  int i;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  lat->GetNeighborRanks(myrank, ranks);

  for (i = 0; i < 6; i++) {

    if (ranks[i] >= 0) {

      MPI_Recv(&(NumRecvPoints[i]), 1, MPI_INT, ranks[i], 0, MPI_COMM_WORLD, &status);

      if (NumRecvPoints[i]) {

	while (SizeRecvPoints[i] < NumRecvPoints[i] * 3 * sizeof(float)) {

	  RecvPoints[i] = (float *)realloc(RecvPoints[i],
					      SizeRecvPoints[i] * 2);
	  if (RecvPoints[i] == NULL)
	    Error("Error: ReceivePoints() cannot reallocate memory\n");

	  SizeRecvPoints[i] *= 2;

	}

	MPI_Recv(RecvPoints[i], NumRecvPoints[i] * 3, MPI_FLOAT, ranks[i], 0, MPI_COMM_WORLD, &status);
	num += NumRecvPoints[i];

      }

    }

  }

  return num;

}
//------------------------------------------------------------------------
//
// PrintPost
//
// prints the posted points
// (debug)
//
void PrintPost() {

  int rank;
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 0; i < 6; i++) {

    if (NumSendPoints[i]) {
      fprintf(stderr, "Rank %d posted %d points for neighbor %d:\n", rank, NumSendPoints[i], i);
      for (j = 0; j < NumSendPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\n",SendPoints[i][j]);
    }

  }

  fprintf(stderr,"\n");

}
//--------------------------------------------------------------------------
//
// PrintRecv
//
// prints the received points
// (debug)
//
void PrintRecv() {

  int rank;
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i = 0; i < 6; i++) {

    if (NumRecvPoints[i]) {
      fprintf(stderr, "Rank %d received %d points from neighbor %d:\n", rank, NumRecvPoints[i], i);
      for (j = 0; j < NumRecvPoints[i]; j++)
	fprintf(stderr, "%.3f\t%.3f\t%.3f\n", RecvPoints[i][j]);
    }

  }

  fprintf(stderr,"\n");

}
//--------------------------------------------------------------------------
//
// Reinit
//
// reinitializes for next iteration
// copies seeds from received points
// clears communication lists
//
void Reinit() {

  int num = 0;
  int i, j;

  // count the total number of seeds and allocate space
  for (i = 0; i < 6; i++)
    num += NumRecvPoints[i];

  while (SizeSeeds < num * sizeof(VECTOR3)) {

    Seeds = (VECTOR3 *)realloc(Seeds, SizeSeeds * 2);
    if (Seeds == NULL)
      Error("Error: Reinit() cannot reallocate memory\n");
    SizeSeeds *= 2;

  }

  // copy seeds
  // it would be good to be able to avoid this step in the future
  num = 0;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < NumRecvPoints[i]; j++)
      Seeds[num++].Set(RecvPoints[i][j + 0], RecvPoints[i][j + 1], 
		     RecvPoints[i][j + 2]);
  }

  // clear communication lists
  for (i = 0; i < 6; i++) {
    NumSendPoints[i] = 0;
    NumRecvPoints[i] = 0;
  }

  NumSeeds = num;

}
//--------------------------------------------------------------------------
