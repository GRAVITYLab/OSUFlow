////////////////////////////////////////////////////////
//
// 3D sample program
//
// Han-Wei Shen
//
////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#endif

#include "OSUFlow.h"

#include <list>
#include <iterator>

int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0; 
float scale_size = 1; 

int xform_mode = 0; 

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListTimeSeedTrace*> sl_list; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 
float center[3], len[3]; 
int first_frame = 1; 

int num_timesteps; 
int num_frames = 50; 
int current_frame = 0; 
float time_incr; 


////////////////////////////////////////////////////////

void compute_pathlines() {

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 128); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  float* tarray = new float[nSeeds]; 
  for (int i=0;i<nSeeds; i++) 
    tarray[i] = (float)(i % num_timesteps); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(1, 5);
  //  osuflow->GenPathLines(seeds, sl_list , FORWARD, nSeeds, 5000); 
  osuflow->GenPathLines(seeds, sl_list , FORWARD, nSeeds, 5000, tarray); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 
}

void draw_pathlines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

  printf("draw pathlines.\n"); 
  glColor3f(1,1,0); 
  std::list<vtListTimeSeedTrace*>::iterator pIter; 

  pIter = sl_list.begin(); 
  for (; pIter!=sl_list.end(); pIter++) {
    vtListTimeSeedTrace *trace = *pIter; 
    std::list<VECTOR4*>::iterator pnIter; 
    pnIter = trace->begin(); 
    glBegin(GL_LINE_STRIP); 
    for (; pnIter!= trace->end(); pnIter++) {
      VECTOR4 p = **pnIter; 
      //printf(" %f %f %f ", p[0], p[1], p[2]); 
      glVertex3f(p[0], p[1], p[2]); 
    }
    glEnd(); 
  }
  glPopMatrix(); 
}

void animate_pathlines() {

  std::list<vtListTimeSeedTrace*>::iterator pIter; 
  vtListTimeSeedTrace *trace; 
  static std::list<VECTOR4*>::iterator pnIter; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

  float min_time = current_frame * time_incr; 
  float max_time = (current_frame+1) * time_incr; 

  glColor3f(1,1,0); 
  pIter = sl_list.begin(); 


  for (; pIter!=sl_list.end(); pIter++) {
    trace = *pIter; 
    pnIter = trace->begin(); 
    glBegin(GL_LINE_STRIP); 
    for (; pnIter!= trace->end(); pnIter++) {
      VECTOR4 p = **pnIter; 
      if (p[3] >= min_time && p[3] < max_time) {
	glColor3f(0,0,1); 
      }
      else glColor3f(0,0,0); 
      glVertex3f(p[0], p[1], p[2]); 
    }
    glEnd(); 
  }


  glPopMatrix(); 
  current_frame = (current_frame+1) % num_frames; 

}

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
  glClearColor(0,0,0,1); 
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
    draw_pathlines(); 
  else if (toggle_animate_streamlines == true)
    animate_pathlines(); 

  printf(" len %f %f %f\n", len[0], len[1], len[2]); 
  glPushMatrix(); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 
  
  glColor3f(1,1,1); 
  glBegin(GL_LINES); 
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0); 
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
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
  glutTimerFunc(50, timer, 0); 
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
	case 's': compute_pathlines(); 
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

int main(int argc, char** argv) 
{
  // read in the vector field 

  VECTOR3 minB, maxB; 

  osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 

  osuflow->LoadData((const char*)argv[1], false); //false : a time-varying flow field 

  osuflow->Boundary(minLen, maxLen); // get the boundary 
  minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
  maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
  osuflow->SetBoundary(minB, maxB);  // set the boundary. just to test
                                     // the subsetting feature of OSUFlow
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                minLen[0], maxLen[0], minLen[1], maxLen[1], 
                                minLen[2], maxLen[2]); 


  num_timesteps = osuflow->NumTimeSteps(); 
  printf(" reading in %d time steps.\n", num_timesteps); 

  time_incr = num_timesteps/(float) num_frames; 

  //  osuflow->NormalizeField(true); 
  osuflow->SetIntegrationParams(1, 5);
  osuflow->ScaleField(30.0);
  osuflow->SetMaxError(0.0001);
  osuflow->SetIntegrationParams(1, 0.01, 5);

  center[0] = (minLen[0]+maxLen[0])/2.0; 
  center[1] = (minLen[1]+maxLen[1])/2.0; 
  center[2] = (minLen[2]+maxLen[2])/2.0; 
  printf("center is at %f %f %f \n", center[0], center[1], center[2]); 
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


