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

bool toggle_draw_streaklines = false; 
bool toggle_animate_streaklines = false; 
float center[3], len[3]; 

OSUFlow *osuflow; 
vtStreakTraces sl_list; 

int first_frame = 1; 
int num_timesteps; 
int num_frames; 
int start_time, end_time; 
int current_frame = 0; 
float time_incr; 

VECTOR3 lMin, lMax; 
VECTOR3 gMin, gMax; 

////////////////////////////////////////////////////////

void compute_streaklines() {

  float from[3], to[3]; 

   from[0] = lMin[0];   from[1] = lMin[1];   from[2] = lMin[2]; 
   to[0] = lMax[0];   to[1] = lMax[1];   to[2] = lMax[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 500); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  float ctime = start_time; 
  sl_list.clear(); 
  printf("compute streaklines..\n"); 
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreakLines(sl_list , FORWARD, ctime); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", sl_list.size()); 

}

/////////////////////////////////////////////////////////////////////////
void draw_streaklines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  //  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  glColor3f(1,1,1); 
  glBegin(GL_POINTS); 
  for (int i=0; i<nSeeds; i++) 
    glVertex3f(seeds[i][0], seeds[i][1], seeds[i][2]); 
  glEnd(); 


  glColor3f(1,1,0); 
  vtStreakTracesIter pIter; 
  pIter = sl_list.begin(); 
  for (; pIter!=sl_list.end(); pIter++) {
    vtListStreakParticle *trace = *pIter; 
    vtStreakParticleIter pnIter; 
    pnIter = trace->begin(); 
    glBegin(GL_LINE_STRIP); 
    for (; pnIter!= trace->end(); pnIter++) {
      vtStreakParticle p = **pnIter; 
      float x = p.itsPoint.phyCoord[0]; 
      float y = p.itsPoint.phyCoord[1]; 
      float z = p.itsPoint.phyCoord[2]; 
      glVertex3f(x,y,z); 
    }
    glEnd(); 
  }
  glPopMatrix(); 
}

void animate_streaklines() {
 
  vtStreakTracesIter pIter; 
  vtListStreakParticle *trace; 
  static vtStreakParticleIter *pnIter; 

  //  glPointSize(2); 
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  //  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  float min_time = current_frame * time_incr; 
  float max_time = (current_frame+1) * time_incr; 

  glColor3f(1,1,0); 
  pIter = sl_list.begin(); 

  glBegin(GL_POINTS); 
  for (; pIter!=sl_list.end(); pIter++) {
    vtListStreakParticle *trace = *pIter; 
    vtStreakParticleIter pnIter; 
    pnIter = trace->begin(); 
    for (; pnIter!= trace->end(); pnIter++) {
      vtStreakParticle p = **pnIter; 
      if (p.itsTime >= min_time && p.itsTime < max_time) {
	float x = p.itsPoint.phyCoord[0]; 
	float y = p.itsPoint.phyCoord[1]; 
	float z = p.itsPoint.phyCoord[2]; 
	glVertex3f(x,y,z); 
      }
    }
  }
  glEnd(); 

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

  if (toggle_draw_streaklines == true)
    draw_streaklines(); 
  else if (toggle_animate_streaklines == true)
    animate_streaklines(); 

  // just draw a simple bounding box, scaled to match the aspect 
  // ration of the field 
  glPushMatrix(); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 
  
  // draw axes 
  glBegin(GL_LINES); 
  glColor3f(1,0,0);   // red: X
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0);
  glColor3f(0,1,0);   // green: Y
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1);   // blue: Z
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 

  glutSwapBuffers(); 
}

void timer(int val) {
  if (toggle_animate_streaklines == true) {
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
	case 's': compute_streaklines(); 
	  glutPostRedisplay(); 
	  break; 
	case 'd': 
	  toggle_draw_streaklines = !toggle_draw_streaklines; 
	  toggle_animate_streaklines = false; 
	  break; 
	case'a': 
	  toggle_animate_streaklines = !toggle_animate_streaklines; 
	  toggle_draw_streaklines = false; 
	  first_frame = 1; 
	}
}
///////////////////////////////////////////////////////////////

int main(int argc, char** argv) 
{
  // create the osuflow object for reading data and computing 
  // streaklines 
  osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 

  // hack the subdomain where the data will be read and 
  // streaklines will be computed 
  // 
  // the spatial subdomain 
  lMin[0] = 0; lMin[1] = 0; lMin[2] = 0; 
  lMax[0] = 127; lMax[1] = 127; lMax[2] = 127; 

  // the temporal subdomain 
  start_time = 0; 
  end_time = 49; 

  // 
  // Now read the data 
  // 
  // osuflow->LoadData((const char*)argv[1], false); //read the whole data 
  osuflow->LoadData((const char*)argv[1], false, lMin, lMax, start_time, 
		    end_time); //read the spatio-temporal subdomain

  //  Set up the animation frame time
  //  num_timesteps = osuflow->NumTimeSteps(); 
  num_timesteps = 50; 
  num_frames = 50; 
  printf(" Animate %d time steps in %d frames.\n", num_timesteps, 
	 num_frames); 
  time_incr = num_timesteps/(float) num_frames; 

  // a hack here. scale the flow field so that particles will run farther 
  osuflow->ScaleField(10.0); 

  // get the global domain bound 
  // osuflow->GetGlobalBounds(gMin, gMax); 
  // hack the global bounds for now 
  gMin[0] = gMin[1] = gMin[2] = 0; 
  gMax[0] = gMax[1] = gMax[2] = 127; 
  printf(" global boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                gMin[0], gMax[0], gMin[1], gMax[1], 
                                gMin[2], gMax[2]); 

  // Parameters for setting up OpenGL transformations
  center[0] = (gMin[0]+gMax[0])/2.0; 
  center[1] = (gMin[1]+gMax[1])/2.0; 
  center[2] = (gMin[2]+gMax[2])/2.0; 
  len[0] = gMax[0]-gMin[0]+1; 
  len[1] = gMax[1]-gMin[1]+1; 
  len[2] = gMax[2]-gMin[2]+1; 

  // glut initialization 
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  glutCreateWindow("Display streaklines"); 
  glutDisplayFunc(display); 
  glutTimerFunc(50, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}


