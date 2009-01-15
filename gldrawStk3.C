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
#include "calc_subvolume.h"
#include "Lattice4D.h"

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

vtStreakTraces *sl_list; 

bool toggle_draw_streaklines = false; 
bool toggle_animate_streaklines = false; 
float center[3], len[3]; 
int first_frame = 1; 

int num_timesteps; 
int num_frames; 
int start_time, end_time; 
int current_frame = 0; 
float time_incr; 

VECTOR3 lMin, lMax; 
VECTOR3 gMin, gMax; 

int nsp = 2, ntp = 5; 
int npart; 
int nproc = 4; 
int total_seeds = 500; 
volume_bounds_type *vb_list; 

OSUFlow **osuflow_list; 
Lattice4D *lat; 

int **plist; 
int *num_partitions; 

////////////////////////////////////////////////////////

void compute_streaklines() {

  float from[3], to[3]; 

  for (int i=0; i<npart; i++) {
    from[0] = vb_list[i].xmin;  
    from[1] = vb_list[i].ymin;     
    from[2] = vb_list[i].zmin; 
    to[0] = vb_list[i].xmax;  
    to[1] = vb_list[i].ymax;     
    to[2] = vb_list[i].zmax; 
    osuflow_list[i]->SetRandomSeedPoints(from, to, total_seeds/npart); 
    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 
       for (int j=0; j<nSeeds; j++) 
	 printf(" domain[%d] seed no. %d : [%f %f %f]\n", i, j, seeds[j][0], 
    	     seeds[j][1], seeds[j][2]); 
  }
  float ctime;
  for (int proc = 0; proc<nproc; proc++) 
     for (int np=0; np<num_partitions[proc]; np++) {
       int i = plist[proc][np]; 
       sl_list[i].clear(); 
       osuflow_list[i]->SetIntegrationParams(1, 5); 
       ctime = vb_list[i].tmin; 
       osuflow_list[i]->GenStreakLines(sl_list[i], FORWARD, ctime, false); 
       printf("domain %d done integrations", i); 
       printf(" %d streaklines. \n", sl_list[i].size()); 
     }
}

/////////////////////////////////////////////////////////////////////////
void draw_streaklines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  //  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  for (int i=0; i<npart; i++) {
    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 
    glColor3f(1,1,1); 
    glBegin(GL_POINTS); 
    for (int j=0; j<nSeeds; j++) 
      glVertex3f(seeds[j][0], seeds[j][1], seeds[j][2]); 
    glEnd(); 

    glColor3f(1,1,0); 
    vtStreakTracesIter pIter; 
    pIter = sl_list[i].begin(); 
    for (; pIter!=sl_list[i].end(); pIter++) {
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
  for (int i=0; i<npart; i++) {
    pIter = sl_list[i].begin(); 
    glBegin(GL_POINTS); 
    for (; pIter!=sl_list[i].end(); pIter++) {
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
  VECTOR3 minLen, maxLen; 

  npart = nsp * ntp; 

  osuflow_list = new OSUFlow*[npart]; 

  // hardcoded dimeisons. the tornado dat is 128^3 

  minLen[0] = minLen[1] = minLen[2] = 0; 
  maxLen[0] = maxLen[1] = maxLen[2] = 127; 

  int ntime = 50; 
  // partition the domain and create a lattice
  lat = new Lattice4D(maxLen[0]-minLen[0]+1, maxLen[1]-minLen[1]+1, 
		      maxLen[2]-minLen[2]+1, 50, 1, 
		      nsp, ntp);  //1 is ghost layer

  vb_list = lat->GetBoundsList(); 
  lat->InitSeedLists();  
  // assign partitions to processors (npart partitions -> nproc processors) 
  lat->RoundRobin_proc(nproc); 

  plist = new int*[nproc]; 
  num_partitions = new int[nproc]; 
  // get each processor's partitions 
  for (int i=0; i<nproc; i++)
    lat->GetPartitions(i, &(plist[i]), num_partitions[i]);

  sl_list = new vtStreakTraces[npart]; //streakline lists, one per subdomain 

  // 
  // now create a list of flow field for the subdomains 
  for (int i=0; i<npart; i++) {
    osuflow_list[i] = new OSUFlow(); 
    printf("Domain(%d):  %d %d %d %d : %d %d %d %d\n", i, vb_list[i].xmin,  
	   vb_list[i].ymin,  vb_list[i].zmin, vb_list[i].tmin, 
	   vb_list[i].xmax,  vb_list[i].ymax,  vb_list[i].zmax, 
	   vb_list[i].tmax); 

    // load subdomain data into OSUFlow
    VECTOR3 minB, maxB; 
    minB[0] = vb_list[i].xmin;  maxB[0] = vb_list[i].xmax;  
    minB[1] = vb_list[i].ymin;  maxB[1] = vb_list[i].ymax;     
    minB[2] = vb_list[i].zmin;  maxB[2] = vb_list[i].zmax; 

    // load the time-varying subdomain data between start_time and end_time 
    osuflow_list[i]->LoadData((const char*)argv[1], false, minB, maxB, 
			      vb_list[i].tmin, vb_list[i].tmax); 

    //    osuflow_list[i]->NormalizeField(true); 
    osuflow_list[i]->ScaleField(10.0); 
  }

  // set up the animation frame time 
  num_timesteps = 50;  // number of time steps in the data 
  num_frames = 50;     // number of frames to loop through the data time 
  time_incr = num_timesteps/(float) num_frames; 

  //  osuflow->GetGlobalBounds(gMin, gMax); 
  gMin[0] = gMin[1] = gMin[2] = 0; 
  gMax[0] = gMax[1] = gMax[2] = 127; 

  // Parameters for setting up some OpenGL transformations 
  center[0] = (gMin[0]+gMax[0])/2.0; 
  center[1] = (gMin[1]+gMax[1])/2.0; 
  center[2] = (gMin[2]+gMax[2])/2.0; 
  len[0] = gMax[0]-gMin[0]+1; 
  len[1] = gMax[1]-gMin[1]+1; 
  len[2] = gMax[2]-gMin[2]+1; 

  // glut functions 
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


