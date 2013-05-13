////////////////////////////////////////////////////////
//
// 3D sample program
//
// Han-Wei Shen
//
////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#else
#include <GL/glut.h> 
#include <GL/gl.h>
#endif

#include <list>
#include <iterator>

#include "OSUFlow.h"
#include "calc_subvolume.h"
#include "Lattice.h"

#include "ComputeFieldLines.h"

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0;
float scale_size = 1;
int xform_mode = 0; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 

float center[3], len[3]; 
int first_frame = 1; 

OSUFlow *osuflow; 
OSUFlow **osuflow_list; 
list<vtListSeedTrace*> *sl_list; 
Lattice* lat; 
volume_bounds_type *vb_list; 

int nproc = 16;   // number of subdomains we will create 

////////////////////////////////////////////////////////

void compute_streamlines() {

  sl_list = ComputeStreamlines(osuflow_list, lat, vb_list, nproc, 40); 

}


////////////////////////////////////////////////

void draw_bounds(float xmin, float xmax, float ymin, float ymax, 
		 float zmin, float zmax)
{
  glColor3f(1,0,0); 
  glBegin(GL_LINES); 

    glVertex3f(xmin, ymin, zmin); glVertex3f(xmax, ymin, zmin); 
    glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymax, zmin); 
    glVertex3f(xmax, ymax, zmin); glVertex3f(xmin, ymax, zmin); 
    glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, xmin, zmin); 

    glVertex3f(xmin, ymin, zmax); glVertex3f(xmax, ymin, zmax); 
    glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax); 
    glVertex3f(xmax, ymax, zmax); glVertex3f(xmin, ymax, zmax); 
    glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, xmin, zmax); 

    glVertex3f(xmin, ymin, zmin); glVertex3f(xmin, ymin, zmax); 
    glVertex3f(xmin, ymin, zmax); glVertex3f(xmin, ymax, zmax); 
    glVertex3f(xmin, ymax, zmax); glVertex3f(xmin, ymax, zmin); 
    glVertex3f(xmin, ymax, zmin); glVertex3f(xmin, ymin, zmin); 

    glVertex3f(xmax, ymin, zmin); glVertex3f(xmax, ymin, zmax); 
    glVertex3f(xmax, ymin, zmax); glVertex3f(xmax, ymax, zmax); 
    glVertex3f(xmax, ymax, zmax); glVertex3f(xmax, ymax, zmin); 
    glVertex3f(xmax, ymax, zmin); glVertex3f(xmax, ymin, zmin); 

  glEnd(); 
}

void draw_streamlines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

  std::list<vtListSeedTrace*>::iterator pIter; 

  for (int i=0; i<nproc; i++) {   // looping through all subdomains 
    pIter = sl_list[i].begin(); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      vtListSeedTrace *trace = *pIter; 
      std::list<VECTOR3*>::iterator pnIter; 
      pnIter = trace->begin(); 
      glColor3f(0.3,0.3,0.3); 
      glBegin(GL_LINE_STRIP); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
	//printf(" %f %f %f ", p[0], p[1], p[2]); 
	glVertex3f(p[0], p[1], p[2]); 
      }
      glEnd(); 
    }
    volume_bounds_type vb = vb_list[i];     
    draw_bounds(vb.xmin, vb.xmax, vb.ymin, vb.ymax, vb.zmin, vb.zmax); 
  }
  glPopMatrix(); 
}

void animate_streamlines() {
}

/*
void animate_streamlines() {

  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  static std::list<VECTOR3*>::iterator *pnIter; 
  static int frame = 0; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

  glColor3f(1,1,0); 

  pIter = sl_list.begin(); 
  int num_lines = sl_list.size(); 
  printf(" animate %d streamlines\n", num_lines); 
  if (first_frame==1) {
    pnIter = new std::list<VECTOR3*>::iterator[num_lines]; 
  }
  int count = 0; 
  int max_len = 0; 
  for (; pIter!=sl_list.end(); pIter++) {
    trace = *pIter; 
    int sz = trace->size(); 
    if (sz> max_len) {
      max_len = sz;
    }
    pnIter[count] = trace->begin(); 
    count++; 
  }
  if (first_frame ==1) {
    frame = 0; 
  }
  else frame = (frame+1)%max_len; 
  printf(" *** max len = %d frame time = %d \n", max_len, frame); 

  pIter = sl_list.begin(); 

  count = 0; 
  for (; pIter!=sl_list.end(); pIter++) {
    trace = *pIter; 
    int sz = trace->size(); 
    //    if (frame >sz) {count++; continue; }
    int frame_count = 0; 
    glBegin(GL_LINE_STRIP); 
    for (; pnIter[count]!= trace->end(); pnIter[count]++) {
      VECTOR3 p = **pnIter[count]; 
      //printf(" %f %f %f ", p[0], p[1], p[2]); 
      glVertex3f(p[0], p[1], p[2]); 
      frame_count++; 
      if (frame_count > frame) break; 
    }
    glEnd(); 
    count++; 
  }
  
  glPopMatrix(); 
  frame++; 
  if (first_frame == 1) first_frame = 0; 
  sleep(.5); 
}
*/

////////////////////////////////////////////// 

void draw_cube(float r, float g, float b)
{
  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube 
}

//////////////////////////////////////////////////////

void display()
{
  glEnable(GL_DEPTH_TEST); 
  glClearColor(1,1,1,1); 
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); 
  gluPerspective(60, 1, .1, 100); 

  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 
  gluLookAt(0,0,5,0,0,0,0,1,0); 

  glRotatef(x_angle, 0, 1,0); 
  glRotatef(y_angle, 1,0,0); 
  glScalef(scale_size, scale_size, scale_size); 

  if (toggle_draw_streamlines == true)
    draw_streamlines(); 
  else if (toggle_animate_streamlines == true)
    animate_streamlines(); 

  glPushMatrix(); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 
  

  glBegin(GL_LINES); 
  glColor3f(1,0,0); 
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0);
  glColor3f(0,1,0);  
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1);  
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 


  glutSwapBuffers(); 
}

void timer(int val) {
  if (toggle_animate_streamlines == true) {
    //    animate_streamlines(); 
    glutPostRedisplay(); 
  }
  glutTimerFunc(10, timer, 0); 
}

///////////////////////////////////////////////////////////

void mymouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    press_x = x; press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE; 
	 else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode = XFORM_SCALE; 
  }
  else if (state == GLUT_UP) {
	  xform_mode = XFORM_NONE; 
  }
}

/////////////////////////////////////////////////////////

void mymotion(int x, int y)
{
    if (xform_mode==XFORM_ROTATE) {
      x_angle += (x - press_x)/5.0; 
      if (x_angle > 180) x_angle -= 360; 
      else if (x_angle <-180) x_angle += 360; 
      press_x = x; 
	   
      y_angle += (y - press_y)/5.0; 
      if (y_angle > 180) y_angle -= 360; 
      else if (y_angle <-180) y_angle += 360; 
      press_y = y; 
    }
	else if (xform_mode == XFORM_SCALE){
      float old_size = scale_size;
      scale_size *= (1+ (y - press_y)/60.0); 
      if (scale_size <0) scale_size = old_size; 
      press_y = y; 
    }
	glutPostRedisplay(); 
}

///////////////////////////////////////////////////////////////

void mykey(unsigned char key, int x, int y)
{
        switch(key) {
	case 'q': exit(1);
	  break; 
	case 's': compute_streamlines(); 
	  glutPostRedisplay(); 
	  break; 
	case 'd': 
	  toggle_draw_streamlines = !toggle_draw_streamlines; 
	  toggle_animate_streamlines = false; 
	  break; 
	case'a': 
	  toggle_animate_streamlines = !toggle_animate_streamlines; 
	  toggle_draw_streamlines = false; 
	  first_frame = 1; 
	}
}

///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////

int main(int argc, char** argv) 
{
  VECTOR3 minLen, maxLen; 

  osuflow_list = new OSUFlow*[nproc];  

  OSUFlow *osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 

  //loading the whole dataset just to get dims. 
  //obviously not very smart. need to change. 
  osuflow->LoadData((const char*)argv[1], true); 
                                                       
  osuflow->Boundary(minLen, maxLen);    // query the dims
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
	 minLen[0], maxLen[0], minLen[1], maxLen[1], 
	 minLen[2], maxLen[2]); 

  // now subdivide the entire domain into nproc subdomains

  // partition the domain and create a lattice
  lat = new Lattice(maxLen[0]-minLen[0+1], maxLen[1]-minLen[1]+1, 
			     maxLen[2]-minLen[2]+1, 0, nproc);
  vb_list = lat->GetBoundsList(); 
  lat->InitSeedLists(); 

  // now create a list of flow field for the subdomains 
  for (int i=0; i<nproc; i++) {
    osuflow_list[i] = new OSUFlow(); 
    printf("Domain(%d):  %d %d %d : %d %d %d\n", i, vb_list[i].xmin,  
	   vb_list[i].ymin,  vb_list[i].zmin, vb_list[i].xmax,  
	   vb_list[i].ymax,  vb_list[i].zmax); 

    // load subdomain data into OSUFlow
    VECTOR3 minB, maxB; 
    minB[0] = vb_list[i].xmin;  
    minB[1] = vb_list[i].ymin;     
    minB[2] = vb_list[i].zmin; 
    maxB[0] = vb_list[i].xmax;  
    maxB[1] = vb_list[i].ymax;     
    maxB[2] = vb_list[i].zmax; 
    osuflow_list[i]->LoadData((const char*)argv[1], true, minB, maxB); 
  }
  vb_list = lat->GetBoundsList(); 

  center[0] = (minLen[0]+maxLen[0])/2.0; 
  center[1] = (minLen[1]+maxLen[1])/2.0; 
  center[2] = (minLen[2]+maxLen[2])/2.0; 
  len[0] = maxLen[0]-minLen[0]; 
  len[1] = maxLen[1]-minLen[1]; 
  len[2] = maxLen[2]-minLen[2]; 

  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  
  glutCreateWindow("Display streamlines"); 
  glutDisplayFunc(display); 
  //  glutIdleFunc(idle); 
  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}


