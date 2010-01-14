//------------------------------------------------------------------------------
//
// field line drawing program
//
// Copyright (c) 2009 Han-Wei Shen and Tom Peterka
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
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>
#include <math.h>

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#ifdef LINUX
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

#ifdef MAC_OSX_10_4
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#endif

// whether or not we have opposite endianness in the file
int swap;

// defines related to drawing
#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 
#define PI 3.1415928						

// drawing state
int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0;
float scale_size = 1.0;
int xform_mode = 0; 
bool toggle_bounds = true;

// drawing data
float *pt; // points in traces
int *npt; // number of points in each trace
int tot_ntrace; // total number of traces
float min[4], max[4]; // extent of domain
float size[4]; // size of domain
float cen[3]; // center of spatial extents
int npt_alloc; // allocated space for this number of traces
int pt_alloc; // allocated space for this number of points
GLUquadricObj *cyl; // cylinder for streamlines
float *rgb; // color for each trace

// function prototypes
void Cleanup();
void draw_bounds(float *from, float *to);
void draw_cube(float r, float g, float b);
void display();
void timer(int val);
void mymouse(int button, int state, int x, int y);
void mymotion(int x, int y);
void mykey(unsigned char key, int x, int y);
void idle();
void init_scene();
void draw_cyl(float *p0, float *p1, float *rgb);
void ReadFieldlines(char *filename);
void swap4(char *n);

//----------------------------------------------------------------------------
//
int main(int argc, char *argv[]) {

  // init
  assert((npt = (int *)malloc(256 * sizeof(int))) != NULL);
  npt_alloc = 256; // initially allocate 256 traces
  assert((pt = (float *)malloc(256 * 4 *sizeof(float))) != NULL);
  pt_alloc = 256; // initially allocate 256 points (each x, y, z, t)

  if (argc < 2) {
    fprintf(stderr, "Usage: draw file_name [swap]\n");
    exit(0);
  }
  if (argc > 2)
    swap = 1;
  else
    swap = 0;

  ReadFieldlines(argv[1]);

  // event loop for drawing
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  glutCreateWindow("Field lines"); 
  glutDisplayFunc(display); 
  glutIdleFunc(idle); 
  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 

  Cleanup();

}
//-----------------------------------------------------------------------
//
// ReadFieldlines
//
void ReadFieldlines(char *filename) {

  FILE *fd;
  int i, j, k, n;

  assert((fd = fopen(filename, "r")) != NULL);

  // extents
  fread(min, sizeof(float), 4, fd);
  fread(max, sizeof(float), 4, fd);
  if (swap) {
    for (i = 0; i < 4; i++) {
      swap4((char *)&min[i]);
      swap4((char *)&max[i]);
    }
  }
  size[0] = max[0] - min[0];
  size[1] = max[1] - min[1];
  size[2] = max[2] - min[2];
  size[3] = max[3] - min[3]  + 1;
  cen[0] = (min[0] + max[0]) * 0.5;
  cen[1] = (min[1] + max[1]) * 0.5;
  cen[2] = (min[2] + max[2]) * 0.5;

  // number of points in each trace
  for (i = 0; ; i++) {
    // check allocated space
    if (i + 1 > npt_alloc) {
      assert((npt = (int *)realloc(npt, npt_alloc * sizeof(int) * 2)) != NULL);
      npt_alloc *= 2;
    }
    fread(&npt[i], sizeof(int), 1, fd);
    if (swap)
      swap4((char *)&npt[i]);
    if (npt[i] < 0)
      break;
  }
  tot_ntrace = i;

  // points
  n = 0;
  for (i = 0; i < tot_ntrace; i++) {
    for (j = 0; j < npt[i]; j++) {
      // check allocated space
      if (n + 1 > pt_alloc) {
	assert((pt = (float *)realloc(pt, pt_alloc * sizeof(float) * 4 * 2))
	       != NULL);
	pt_alloc *= 2;
      }
      fread(&pt[n * 4], sizeof(float), 4, fd);
      if (swap) {
	for (k = 0; k < 4; k++)
	  swap4((char *)&pt[n * 4 + k]);
      }
      n++;
    }
  }

  fclose(fd);

}
//-----------------------------------------------------------------------
//
// Cleanup
//
// frees memory and such
//
void Cleanup() {

  free(npt);
  free(pt);

}
//-----------------------------------------------------------------------
//
void draw_bounds(float *from, float *to) {

  float xmin = from[0];
  float ymin = from[1];
  float zmin = from[2];
  float xmax = to[0];
  float ymax = to[1];
  float zmax = to[2];

  glPushMatrix(); 

  glScalef(1.0f / (float)size[0], 1.0f / (float)size[0], 1.0f / (float)size[0]);
  glTranslatef(-cen[0], -cen[1], -cen[2]);

  glColor3f(1,0,0); 
  glBegin(GL_LINES); 

  glVertex3f(xmin, ymin, zmin); glVertex3f(xmax, ymin, zmin); 
  glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymax, zmin); 
  glVertex3f(xmax, ymax, zmin); glVertex3f(xmin, ymax, zmin); 
  glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, ymin, zmin); 

  glVertex3f(xmin, ymin, zmax); glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax); 
  glVertex3f(xmax, ymax, zmax); glVertex3f(xmin, ymax, zmax); 
  glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, ymin, zmax); 

  glVertex3f(xmin, ymin, zmin); glVertex3f(xmin, ymin, zmax); 
  glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, ymax, zmin); 

  glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmax, ymax, zmax); glVertex3f(xmax, ymax, zmin); 

  glEnd(); 

  glPopMatrix();

}
//--------------------------------------------------------------------------
//
void DrawFieldlines() {
  
  static int step = 0; // timestep number
  static int frame_num = 0; // number of frames at this timestep
  static int frames_per_step; // hold each time step for this many frames
  int tsize = 1; // default for now
  int start; // start new trace
  int i, j, k;
  float *p0, *p1; // pointers to points

  // compute hold time
  frames_per_step = 2000 - tsize;

  glPushMatrix(); 

  glScalef(1.0f / size[0], 1.0f / size[0], 1.0f / size[0]);
  glTranslatef(-cen[0], -cen[1], -cen[2]);
  glColor3f(1.0, 1.0, 0.0); 

  // lines

//   k = 0;

//   for (i = 0; i < tot_ntrace; i++) {

//     glBegin(GL_LINE_STRIP); 
//     for (j = 0; j < npt[i]; j++) {
// //       if (pt[k * 4 + 3] <= step) {

// 	// artificial boundary to see only center region
// 	if (pt[k * 4 + 0] >= -3.072e8 && pt[k * 4 + 0] <= 3.072e8 && 
// 	pt[k * 4 + 1] >= -5.12e7 && pt[k * 4 + 1] <= 5.12e7 && 
// 	pt[k * 4 + 2] >= -3.072e8 && pt[k * 4 + 2] <= 3.072e8)

// 	  glVertex3f(pt[k * 4 + 0], pt[k * 4 + 1], pt[k * 4 + 2]);

// //       }
//       k++;
//     }
//     glEnd();

//   }

  // cylinders

  k = 0;

  for (i = 0; i < tot_ntrace; i++) {

    start = 1;
    for (j = 0; j < npt[i]; j++) {
//       if (pt[k * 4 + 3] <= step) {

      // artificial boundary to see only center region
      if (pt[k * 4 + 0] >= -3.072e8 && pt[k * 4 + 0] <= 3.072e8 && 
	  pt[k * 4 + 1] >= -5.12e7 && pt[k * 4 + 1] <= 5.12e7 && 
	  pt[k * 4 + 2] >= -3.072e8 && pt[k * 4 + 2] <= 3.072e8) {

	if (start) {
	  p0 = &(pt[k * 4]);
	  start = 0;
	}
	else {
	  p1 = &(pt[k * 4]);
	  draw_cyl(p0, p1, &rgb[i]);
	  p0 = p1;
	}

      } // artificial boundary


//       }
      k++;
    }

  }

  glPopMatrix(); 

  if (frame_num == frames_per_step) {
    frame_num = 0;
    if (step == tsize - 1)
      step = 0;
    else
      step++;
  }
  else
    frame_num++;

}
//-------------------------------------------------------------------------
//
void draw_cyl(float *p0, float *p1, float *rgb) {

  float len; // length of cylinder
  float lenxy; // length projected to xy plane
  float zrot, xrot; // rotation angles about z, x (degrees)
  float rad; // cylinder radius
  
  // some math prerequisites
  len = sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + 
	     (p1[1] - p0[1]) * (p1[1] - p0[1]) + 
	     (p1[2] - p0[2]) * (p1[2] - p0[2]));  
  lenxy = sqrt((p1[0] - p0[0]) * (p1[0] - p0[0]) + 
	       (p1[1] - p0[1]) * (p1[1] - p0[1]));
  xrot = atan2((p1[2] - p0[2]), lenxy); // degrees
  zrot = atan2((p1[1] - p0[1]), (p1[0] - p0[0]));
  xrot = xrot * 180.0 / PI - 90.0; // rads
  zrot = zrot * 180.0 / PI - 90.0;
  rad = (size[0] + size[1] + size[2]) / 3.0 / 7000.0;
//   rad = (size[0] + size[1] + size[2]) / 3.0 / 700.0;
  
  glColor3f(rgb[0], rgb[1], rgb[2]);

  glPushMatrix();
  glTranslatef(p0[0], p0[1], p0[2]);
  glRotatef(zrot, 0.0, 0.0, 1.0);
  glRotatef(xrot, 1.0, 0.0, 0.0);
  gluCylinder(cyl, rad, rad, len, 10, 1);
  glPopMatrix();

}
//-------------------------------------------------------------------------
//
void display() {

  static int first = 1;

  // init
  if (first) {
    init_scene();
    first = 0;
  }

  // transform
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 
  gluPerspective(60, 1, .1, 100); 

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(0,0,5,0,0,0,0,1,0); 

  glRotatef(x_angle, 0, 1,0); 
  glRotatef(y_angle, 1,0,0);
  glScalef(scale_size, scale_size, scale_size); 

  // clear
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 

  // field lines
  DrawFieldlines(); 

  // axes and bounds
  if (toggle_bounds) {  
    glBegin(GL_LINES); 
    glColor3f(1,0,0); // red x-axis
    glVertex3f(0,0,0); 
    glVertex3f(1,0,0);
    glColor3f(0,1,0);  // green y-axis
    glVertex3f(0,0,0);
    glVertex3f(0,1,0); 
    glColor3f(0,0,1);  // blue z-axis
    glVertex3f(0,0,0); 
    glVertex3f(0,0,1); 
    glEnd(); 
    draw_bounds(min, max);
  }

  glutSwapBuffers(); 

}
//--------------------------------------------------------------------------
//
void init_scene() {

  int i;

  assert((rgb = (float *)malloc(tot_ntrace * 3 * sizeof(float))) != NULL);

  // background
  glClearColor(0.0, 0.0, 0.0, 1.0); 

  // modes
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_DEPTH_TEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_NORMALIZE);

  // lights
  GLfloat light_ambient[4] = {0.05, 0.05, 0.05, 1.0};  
  GLfloat light_diffuse[4] = {0.6, 0.6, 0.6, 1.0};  
  GLfloat light_specular[4] = {0.8, 0.8, 0.8, 1.0};
  GLfloat light_position[4] = {1.0, 1.0, 1.0, 0.0};
  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  light_position[0] = -1.0;
  light_position[1] = -0.2;
  light_position[2] = -0.1;
  glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_POSITION, light_position);

  // init quadric
  cyl = gluNewQuadric();
  gluQuadricDrawStyle(cyl, GLU_FILL);
  gluQuadricNormals(cyl, GLU_SMOOTH);

  // generate color randomly for now for each trace
  for (i = 0; i < tot_ntrace; i++) {
    rgb[3 * i]     = 0.33 + rand() / (float)RAND_MAX * 0.66;
    rgb[3 * i + 1] = 0.33 + rand() / (float)RAND_MAX * 0.66;
    rgb[3 * i + 2] = 0.33 + rand() / (float)RAND_MAX * 0.66;
  }

}
//--------------------------------------------------------------------------
//
void timer(int val) {

  glutTimerFunc(10, timer, 0); 

}
//--------------------------------------------------------------------------
//
void idle() {

  glutPostRedisplay();

}
//--------------------------------------------------------------------------
//
void mymouse(int button, int state, int x, int y) {

  if (state == GLUT_DOWN) {

    press_x = x; press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE;
    else if (button == GLUT_RIGHT_BUTTON)
      xform_mode = XFORM_SCALE; 

  }

  else if (state == GLUT_UP)
    xform_mode = XFORM_NONE;

}
//--------------------------------------------------------------------------
//
void mymotion(int x, int y) {

  if (xform_mode==XFORM_ROTATE) {

    x_angle += (x - press_x)/5.0; 
    if (x_angle > 180)
      x_angle -= 360; 
    else if (x_angle <-180)
      x_angle += 360; 
    press_x = x; 
	   
    y_angle += (y - press_y)/5.0; 
    if (y_angle > 180)
      y_angle -= 360; 
    else if (y_angle <-180)
      y_angle += 360; 
    press_y = y; 

  }

  else if (xform_mode == XFORM_SCALE){

    float old_size = scale_size;
    scale_size *= (1+ (y - press_y)/60.0); 
    if (scale_size <0)
      scale_size = old_size; 
    press_y = y; 

  }

  glutPostRedisplay(); 

}
//--------------------------------------------------------------------------
//
void mykey(unsigned char key, int x, int y) {

  switch(key) {

  case 'q':
    exit(1);
    break; 
  case 'b':
    toggle_bounds = !toggle_bounds;
    break;
  default:
    break;

  }

}
//--------------------------------------------------------------------------
//
// swap4(n)
//
// Swaps 4 bytes from 1-2-3-4 to 4-3-2-1 order.
// cast the input as a char and use on any 4 byte variable
//
void swap4(char *n) {

  char *n1;
  char c;

  n1 = n + 3;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//----------------------------------------------------------------------------
