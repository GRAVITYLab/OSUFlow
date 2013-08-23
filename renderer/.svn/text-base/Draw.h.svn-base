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

class CDynamicColorLineRenderer: public CLineRendererInOpenGL {

  vector<bool> vbTraceFlags;
  vector<VECTOR4> vv4TraceColors;
  virtual void _CheckTrace(int iTrace, bool& bIsDrawingTrace) {
    bIsDrawingTrace = vbTraceFlags[iTrace];
  }

  virtual void _GetTraceColor(int iTrace, float& fR, float& fG, float& fB, 
			      float& fA) {
    VECTOR4 v4Color = vv4TraceColors[iTrace];
    fR = v4Color[0];
    fG = v4Color[1];
    fB = v4Color[2];
    fA = v4Color[3];
  }

 public:
  void _TraverseLinesBegin(int iNrOfTraces) {
    vbTraceFlags.resize(iNrOfTraces);
    vv4TraceColors.resize(iNrOfTraces);
    CLineRendererInOpenGL::_TraverseLinesBegin(iNrOfTraces);
  }
	
  void _GenerateRandomFlags() {
    for(int t = 0; t < vbTraceFlags.size(); t++)
	vbTraceFlags[t] = (rand()%2)?true:false;
  }

  void _SetAllFlags() {
    for(int t = 0; t < vbTraceFlags.size(); t++)
	vbTraceFlags[t] = true;
  }

  void _GenerateRandomColor() {
    for(int t = 0; t < vbTraceFlags.size(); t++) {
      VECTOR4 v4Color;
      int iRand = rand() % 7;
      switch(iRand) {
      case 0: v4Color = VECTOR4(0.8f, 0.0f, 0.0f, 1.0f); break;
      case 1: v4Color = VECTOR4(0.0f, 0.8f, 0.0f, 1.0f); break;
      case 2: v4Color = VECTOR4(0.0f, 0.0f, 0.8f, 1.0f); break;
      case 3: v4Color = VECTOR4(0.8f, 0.8f, 0.0f, 1.0f); break;
      case 4: v4Color = VECTOR4(0.8f, 0.0f, 0.8f, 1.0f); break;
      case 5: v4Color = VECTOR4(0.0f, 0.8f, 0.8f, 1.0f); break;
      case 6: v4Color = VECTOR4(0.8f, 0.8f, 1.0f, 1.0f); break;
      }
      vv4TraceColors[t] = v4Color;
    }
  }

};

class CDynamicColorTubeRenderer: public CTubeRendererInOpenGL {

  vector<bool> vbTraceFlags;
  vector<VECTOR4> vv4TraceColors;
  virtual void _CheckTrace(int iTrace, bool& bIsDrawingTrace) {
    bIsDrawingTrace = vbTraceFlags[iTrace];
  }

  virtual void _GetTraceColor(int iTrace, float& fR, float& fG, float& fB, 
			      float& fA) {
    VECTOR4 v4Color = vv4TraceColors[iTrace];
    fR = v4Color[0];
    fG = v4Color[1];
    fB = v4Color[2];
    fA = v4Color[3];
  }

 public:
  void _TraverseLinesBegin(int iNrOfTraces) {
    vbTraceFlags.resize(iNrOfTraces);
    vv4TraceColors.resize(iNrOfTraces);
    CTubeRendererInOpenGL::_TraverseLinesBegin(iNrOfTraces);
  }

  void _GenerateRandomFlags() {
    for(int t = 0; t < vbTraceFlags.size(); t++)
      vbTraceFlags[t] = (rand()%2)?true:false;
  }

  void _SetAllFlags() {
    for(int t = 0; t < vbTraceFlags.size(); t++)
	vbTraceFlags[t] = true;
  }

  void _GenerateRandomColor() {
    for(int t = 0; t < vbTraceFlags.size(); t++) {
      VECTOR4 v4Color;
      int iRand = rand() % 7;
      switch(iRand)
	{
	case 0: v4Color = VECTOR4(0.8f, 0.0f, 0.0f, 1.0f);	break;
	case 1: v4Color = VECTOR4(0.0f, 0.8f, 0.0f, 1.0f);	break;
	case 2: v4Color = VECTOR4(0.0f, 0.0f, 0.8f, 1.0f);	break;
	case 3: v4Color = VECTOR4(0.8f, 0.8f, 0.0f, 1.0f);	break;
	case 4: v4Color = VECTOR4(0.8f, 0.0f, 0.8f, 1.0f);	break;
	case 5: v4Color = VECTOR4(0.0f, 0.8f, 0.8f, 1.0f);	break;
	case 6: v4Color = VECTOR4(0.8f, 0.8f, 0.8f, 1.0f);	break;
	}
      vv4TraceColors[t] = v4Color;
    }
  }

};

// function prototypes
void DrawInit(VECTOR4 *pt, int *npt, int ntrace, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb = NULL);
void DrawInit(float *pt, int *npt, int ntrace, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb = NULL);
void DrawInit(const vector<float>& pt, const vector<int>& np, int ntrace, int argc,
	      char **argv, float *min, float *max, int tubes, float *rgb = NULL);

void display();
void key(unsigned char key, int x, int y);
void render_init();
void gl_init();
void quit();
void trace_lists();

#endif

#endif
