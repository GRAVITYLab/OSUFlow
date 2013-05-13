//------------------------------------------------------------------------------
//
// drawing header
// performs drawing functions
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

#ifndef _DRAW_H_
#define _DRAW_H_

#ifdef GRAPHICS

#include <list>
#include <iterator>

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#ifdef LINUX
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

#include "VectorMatrix.h"
#include "libgcb/gcb.h"
#include "OSUFlow.h"
#include "LineRendererInOpenGL.h"
#include "TubeRendererInOpenGL.h"

// function prototypes
void DrawInit(VECTOR4 *pt, int *npt, int ntrace, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb = NULL);
void DrawInit(float *pt, int *npt, int ntrace, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb = NULL);
#if	0	// MOD-BY-LEETEN 12/08/201-FROM:
	void DrawInit(vector<float> pt, vector<int> np, int ntrace, int argc,
			  char **argv, float *min, float *max, int tubes, float *rgb = NULL);
#else	// MOD-BY-LEETEN 12/08/201-TO:
void DrawInit(const vector<float>& pt, const vector<int>& np, int ntrace, int argc,
	      char **argv, float *min, float *max, int tubes, float *rgb = NULL);
#endif	// MOD-BY-LEETEN 12/08/201-END

void display();
void key(unsigned char key, int x, int y);
void render_init();
void gl_init();
void quit();
void trace_lists();

#endif

#endif
