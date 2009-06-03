
#include <stdio.h>
#include <stdlib.h> 

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#else
#include <GL/glut.h> 
#include <GL/gl.h>
#endif 

#include "FlashAMR.h" 

int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0;
float scale_size = 1; 

int xform_mode = 0; 

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

FlashAMR *amr; 

int nb; 
int num_levels;

float len[3];     // domain dimensions 
float center[3];  // domain center 

int draw_level=0; 

////////////////////////////////////////////// 

void draw_cube(float r, float g, float b)
{
  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube 
}

///////////////////////////////////////////////

void draw_bbx(int i) {

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 


  float xsize, ysize, zsize; 
  int its_level = amr->GetLevel(i); 

  if (its_level%3 == 0) glColor3f(1,0,0); 
  if (its_level%3 == 1) glColor3f(0,1,0); 
  if (its_level%3 == 2) glColor3f(0,0,1); 

  float *center = amr->GetBlockCenter(i); 
  float *length = amr->GetBlockLengths(i); 

  glPushMatrix(); 
  glTranslatef(center[0], center[1], center[2]); 
  glScalef(length[0], length[1], length[2]); 

  glutWireCube(1.0); 
  glPopMatrix(); 

  glPopMatrix(); 

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

  glPushMatrix(); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 

  for (int i=0; i<nb; i++) {
    if (amr->GetLevel(i) >= draw_level) 
      draw_bbx(i); 
  }
    
  glBegin(GL_LINES); 
  glColor3f(1,0,0); 
  glVertex3f(0,0,0); 
  glVertex3f(1,0,0); 
  glColor3f(0,1, 0); 
  glVertex3f(0,0,0);
  glVertex3f(0,1,0); 
  glColor3f(0,0,1); 
  glVertex3f(0,0,0); 
  glVertex3f(0,0,1); 
  glEnd(); 

  glutSwapBuffers(); 
}


///////////////////////////////////////////////////////////////

void mykey(unsigned char key, int x, int y)
{
        switch(key) {
	case 'q': exit(1);
	  break; 
	case 'l': 
	  draw_level++; 
	  draw_level = draw_level % (num_levels); 
	  break; 
	}
	glutPostRedisplay(); 
}

////////////////////////////////////////////////////////

main(int argc, char* argv[]) {

  float min[3], max[3]; 
  
  amr = new FlashAMR; 
  amr->LoadData(argv[1], min, max); 

  nb = amr->GetNumBlocks(); 
  num_levels = amr->GetNumLevels(); 

  len[0] = max[0]-min[0]; 
  len[1] = max[1]-min[1]; 
  len[2] = max[2]-min[2]; 
  center[0] = (min[0]+max[0])/2.0; 
  center[1] = (min[1]+max[1])/2.0; 
  center[2] = (min[2]+max[2])/2.0; 

  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  
  glutCreateWindow("Display streamlines"); 
  glutDisplayFunc(display); 

  //    glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}
