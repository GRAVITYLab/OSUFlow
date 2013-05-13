//------------------------------------------------------------------------------
//
// drawing module
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

#ifdef GRAPHICS

#include "Draw.h"

// program args
static int argn;
static char **args;

// drawing data
static vector<VECTOR4> vpt; // points in everyone's traces, vec4 format
static vector<float> fpt; // points in everyone's traces, float format
static vector<int> npt; // everyone's number of points in their traces
static int tot_ntrace; // total number of everyone's traces
static list<vtListSeedTrace*> sl_lists; //list of traces needed by renderer lib
static list<VECTOR4> cl_list; // list of colors, one for each trace
static VECTOR3 minLen, maxLen;
static CLineRendererInOpenGL LineRenderer;
static CTubeRendererInOpenGL TubeRenderer;
static bool tube = 0; // tubes or lines
static float col[3]; // rendering color
//----------------------------------------------------------------------------
//
// initializes drawing, vector4 array format
//
// p: points, vector4 format
// np: number of points in each trace
// nt: number of traces
// argc, argv: program args
// min, max: domain boundary
// tubes: draw tubes or lines
//
void DrawInit(VECTOR4 *p, int *np, int nt, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb) {

  int i, j, n;

  n = 0;
  for (i = 0; i < nt; i++) {
    npt.push_back(np[i]);
    for (j = 0; j < npt[i]; j++)
      vpt.push_back(p[n++]);
  }
  tot_ntrace = nt;
  argn = argc;
  args = argv;
  minLen = min;
  maxLen = max;
  tube = tubes;
  if (rgb != NULL) {
    col[0] = rgb[0]; col[1] = rgb[1]; col[2] = rgb[2];
  }
  else { // default color
    col[0] = 0.6; col[1] = 0.75; col[2] = 0.9;
  }
  render_init();
  trace_lists();
  glutPostRedisplay();

}
//----------------------------------------------------------------------------
//
// initializes drawing, vector4 array format
//
// p: points, vector4 format
// np: number of points in each trace
// nt: number of traces
// argc, argv: program args
// min, max: domain boundary
// tubes: draw tubes or lines
//
void DrawInit(float *p, int *np, int nt, int argc,
	      char **argv, VECTOR3 min, VECTOR3 max, int tubes, float *rgb) {

  VECTOR4 vp;
  int i, j, k, n;

  n = 0;
  for (i = 0; i < nt; i++) {
    npt.push_back(np[i]);
    for (j = 0; j < npt[i]; j++) {
      vp.Set(p[n], p[n + 1], p[n + 2], p[n + 3]);
      vpt.push_back(vp);
      n += 4;
    }
  }
  tot_ntrace = nt;
  argn = argc;
  args = argv;
  minLen = min;
  maxLen = max;
  tube = tubes;
  if (rgb != NULL) {
    col[0] = rgb[0]; col[1] = rgb[1]; col[2] = rgb[2];
  }
  else { // default color
    col[0] = 0.6; col[1] = 0.75; col[2] = 0.9;
  }
  render_init();
  trace_lists();
  glutPostRedisplay();

}
//----------------------------------------------------------------------------
//
// initializes drawing, float array format
//
// p: points, float format
// np: number of points in each trace
// nt: number of traces
// argc, argv: program args
// min, max: domain boundary
// tubes: draw tubes or lines
//
#if	0	// MOD-BY-LEETEN 12/08/201-FROM:
void DrawInit(vector<float>p, vector<int>np, int nt, int argc,
	      char **argv, float *min, float *max, int tubes, float *rgb) {
#else	// MOD-BY-LEETEN 12/08/201-TO:
void DrawInit(const vector<float>& p, const vector<int>& np, int nt, int argc,
	      char **argv, float *min, float *max, int tubes, float *rgb) {
#endif	// MOD-BY-LEETEN 12/08/201-END

  fpt = p;
  npt = np;
  tot_ntrace = nt;
  argn = argc;
  args = argv;
  minLen = VECTOR3(min[0], min[1], min[2]);
  maxLen = VECTOR3(max[0], max[1], max[2]);
  tube = tubes;
  if (rgb != NULL) {
    col[0] = rgb[0]; col[1] = rgb[1]; col[2] = rgb[2];
  }
  else { // default color
    col[0] = 0.6; col[1] = 0.75; col[2] = 0.9;
  }
  render_init();
  trace_lists();
  glutPostRedisplay();

}
//----------------------------------------------------------------------------
//
void display() {

  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
  glPushAttrib(GL_LIGHTING_BIT | 0);
  int lighting;
  if (tube)
    TubeRenderer._GetInteger(CTubeRenderer::ENABLE_LIGHTING, &lighting);
  else
    LineRenderer._GetInteger(CLineRenderer::ENABLE_LIGHTING, &lighting);
  if(lighting && tube) {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glPushMatrix();
    glLoadIdentity();
    static GLfloat pfLightPos[4] = {0.0f, 0.0f, 1.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, pfLightPos);
    glPopMatrix();
  }
  if (tube)
    TubeRenderer._Draw();
  else
    LineRenderer._Draw();
  glPopAttrib();
  glutSwapBuffers(); 

}
//--------------------------------------------------------------------------
//
void key(unsigned char key, int x, int y) {

  switch(key) {

  case 'q':
    exit(1);
    break; 
  case 'h':
    int halo;
    if (tube)
      TubeRenderer._GetInteger(CTubeRenderer::ENABLE_HALO, &halo);
    else
      LineRenderer._GetInteger(CLineRenderer::ENABLE_HALO, &halo);
    if (tube)
      TubeRenderer._SetInteger(CTubeRenderer::ENABLE_HALO, !halo);
    else
      LineRenderer._SetInteger(CLineRenderer::ENABLE_HALO, !halo);
    glutPostRedisplay();
    break;
  case 'l':
    int lighting;
    if (tube)
      TubeRenderer._GetInteger(CTubeRenderer::ENABLE_LIGHTING, &lighting);
    else
      LineRenderer._GetInteger(CLineRenderer::ENABLE_LIGHTING, &lighting);
    if(tube)
      TubeRenderer._SetInteger(CTubeRenderer::ENABLE_LIGHTING, !lighting);
    else
      LineRenderer._SetInteger(CLineRenderer::ENABLE_LIGHTING, !lighting);
    glutPostRedisplay();
    break;
  case 's':
    trace_lists();
    glutPostRedisplay();
    break;
  default:
    glutPostRedisplay();
    break;
  }

}
//--------------------------------------------------------------------------
//
void render_init() {

  glutInit(&argn, args); 
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | 
		      GLUT_ALPHA | GLUT_STENCIL);

  if (tube) {
    TubeRenderer._SetBoundingBox(minLen[0], minLen[1], minLen[2], 
				  maxLen[0], maxLen[1], maxLen[2]);
    TubeRenderer._SetDataSource(&sl_lists);
    TubeRenderer._SetColorSource(&cl_list);
    TubeRenderer._SetInteger(CTubeRenderer::COLOR_SCHEME, 
			      CTubeRenderer::CColorScheme::COLOR_PER_TRACE);
	// ADD-BY-LEETEN 12/08/201-BEGIN
	int iMaxWindowSize = 256;
	float fLineWidth = (maxLen[0] - minLen[0])/(float)iMaxWindowSize;
	TubeRenderer._SetFloat(CLineRenderer::LINE_WIDTH, fLineWidth);
	// ADD-BY-LEETEN 12/08/201-END
  }

  else {
    LineRenderer._SetBoundingBox(minLen[0], minLen[1], minLen[2], 
				  maxLen[0], maxLen[1], maxLen[2]);
    LineRenderer._SetDataSource(&sl_lists);
    LineRenderer._SetColorSource(&cl_list);
    LineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, 
			      CLineRenderer::CColorScheme::COLOR_PER_TRACE);
  }

  glutCreateWindow("Fieldlines"); 
  gcbInit(gl_init, quit);
  gcbDisplayFunc(display); 
  gcbKeyboardFunc(key);

  glutMainLoop(); // run

}
//--------------------------------------------------------------------------
//
void gl_init() {

  trace_lists();

  glEnable(GL_DEPTH_TEST);
  static GLfloat pfLightAmbient[4] =	{0.1f, 0.1f, 0.1f, 1.0f};
  static GLfloat pfLightDiffuse[4] =	{0.6f, 0.6f, 0.6f, 1.0f};
  static GLfloat pfLightSpecular[4] =	{0.6f, 0.6f, 0.6f, 1.0f};
  glLightfv(GL_LIGHT0, GL_AMBIENT,		pfLightAmbient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,		pfLightDiffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR,		pfLightSpecular);
  glLightf(GL_LIGHT0, GL_SPOT_EXPONENT,	4.0f);
  if (!tube)
    LineRenderer._UpdateLighting();

}
//--------------------------------------------------------------------------
void quit() {}
//--------------------------------------------------------------------------
//
// packages traces into lists amenable for renderer lib
//
void trace_lists() {

  int i, j;
  int n = 0;
  VECTOR3 *pt3;
  vtListSeedTrace *trace;
  VECTOR4 color;

  color = VECTOR4(col[0], col[1], col[2], 1.0f); // constant color for now
  sl_lists.clear();

  for (i = 0; i < tot_ntrace; i++) {

    trace = new vtListSeedTrace;

    for (j = 0; j < npt[i]; j++) {
      if (vpt.size()) {
	if (vpt[n][0] >= minLen[0] && vpt[n][0] <= maxLen[0] &&
	    vpt[n][1] >= minLen[1] && vpt[n][1] <= maxLen[1] &&
	    vpt[n][2] >= minLen[2] && vpt[n][2] <= maxLen[2]) {
	  pt3 = new VECTOR3(vpt[n][0], vpt[n][1], vpt[n][2]);
	  trace->push_back(pt3);
	}
	n++;
      }
      else {
	if (fpt[n] >= minLen[0] && fpt[n] <= maxLen[0] &&
	    fpt[n + 1] >= minLen[1] && fpt[n + 1] <= maxLen[1] &&
	    fpt[n + 2] >= minLen[2] && fpt[n + 2] <= maxLen[2]) {
	  pt3 = new VECTOR3(fpt[n], fpt[n + 1], fpt[n + 2]);
	  trace->push_back(pt3);
	}
	n += 4;
      }

    }

    sl_lists.push_back(trace);
    cl_list.push_back(color);

  }

  if (tube)
    TubeRenderer._Update();
  else
    LineRenderer._Update();

}
//--------------------------------------------------------------------------

#endif // GRAPHICS

