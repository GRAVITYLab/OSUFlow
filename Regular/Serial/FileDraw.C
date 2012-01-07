//------------------------------------------------------------------------------
//
// field line drawing program
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

// ADD-BY-LEETEN 01/06/2012-BEGIN
#include <vector>
using namespace std;
// ADD-BY-LEETEN 01/06/2012-END

#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>
#include <math.h>

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#ifdef LINUX
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

#include "Draw.h"
#include "leeten_Draw.h"	// ADD-BY-LEETEN 01/06/2012

// whether or not we have opposite endianness in the file
int byte_swap;

// drawing data
vector<float> pt; // points in traces
vector<int> npt; // number of points in each trace
int tot_ntrace; // total number of traces
float min_ex[4], max_ex[4]; // extent of domain

// function prototypes
void ReadFieldlines(char *filename);

//----------------------------------------------------------------------------
//
int main(int argc, char *argv[]) {

  float rgb[3];
  int tubes;

  // parse args
  if (argc < 4) {
    fprintf(stderr, "Usage: draw file_name tubes/lines<t/l> swap<0/1> [r<0.6> g<0.75> b<0.9>]\n");
    exit(0);
  }
  if (*argv[2] =='t')
    tubes = 1;
  else
    tubes = 0;
  byte_swap = atoi(argv[3]);
  if (argc >= 7) {
    rgb[0] = atof(argv[4]);
    rgb[1] = atof(argv[5]);
    rgb[2] = atof(argv[6]);
  }

  ReadFieldlines(argv[1]);

  if (argc >= 7)
    DrawInit(pt, npt, tot_ntrace, argc, argv, min_ex, max_ex, tubes, rgb);
  else
    DrawInit(pt, npt, tot_ntrace, argc, argv, min_ex, max_ex, tubes, NULL);

}
//------------------------------------------------------------------------------
//
// ReadFieldlines
//
void ReadFieldlines(char *filename) {

  FILE *fd;
  int i, j, k;
  int n;
  float p[4];

  assert((fd = fopen(filename, "r")) != NULL);

  // extents
  assert(fread(min_ex, sizeof(float), 4, fd) == 4);
  assert(fread(max_ex, sizeof(float), 4, fd) == 4);
  if (byte_swap) {
    for (i = 0; i < 4; i++) {
      swap4((char *)&min_ex[i]);
      swap4((char *)&max_ex[i]);
    }
  }

  // number of points in each trace
  for (i = 0; ; i++) {
    assert(fread(&n, sizeof(int), 1, fd) == 1);
    if (byte_swap)
      swap4((char *)&n);
    if (n < 0)
      break;
    npt.push_back(n);
  }
  tot_ntrace = i;

  // points
  for (i = 0; i < tot_ntrace; i++) {
    for (j = 0; j < npt[i]; j++) {
      assert(fread(&p, sizeof(float), 4, fd) == 4);
      if (byte_swap) {
	for (k = 0; k < 4; k++)
	  swap4((char *)&p[k]);
      }
      for (k = 0; k < 4; k++)
	pt.push_back(p[k]);
    }
  }

  fclose(fd);

}
//-----------------------------------------------------------------------
//
// this is needed to link in the correct order under linux
// apparently I need to access an object and a function in the anlcom lib
// to link renderer and Core in correct order.
// this does not get called or do anything at runtime
//
#include "Blocks.h"
Blocks *blocks;
void noop(){
  delete blocks;
}
//------------------------------------------------------------------------------
