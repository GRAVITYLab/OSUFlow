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

// ADD-BY-LEETEN 12/20/2011-BEGIN
#ifdef	WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif	// #ifdef WIN32
// ADD-BY-LEETEN 12/20/2011-END

// VTK
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkSmartPointer.h"

#include "OSUFlow.h"
#include "VectorFieldVTK.h"


// VTK

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
list<vtListSeedTrace*> sl_list; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 
float center[3], len[3]; 
int first_frame = 1; 

////////////////////////////////////////////////////////

void compute_streamlines() {

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 300);
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(.001, .001);
  osuflow->GenStreamLines(sl_list , BACKWARD_AND_FORWARD, 500, 0);
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 
}

void draw_streamlines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]);

  printf("draw streamlines.\n"); 
  glColor3f(1,1,0); 
  std::list<vtListSeedTrace*>::iterator pIter; 

  pIter = sl_list.begin(); 
  for (; pIter!=sl_list.end(); pIter++) {
    vtListSeedTrace *trace = *pIter; 
    std::list<VECTOR3*>::iterator pnIter; 
    pnIter = trace->begin(); 
    glBegin(GL_LINE_STRIP); 
    for (; pnIter!= trace->end(); pnIter++) {
      VECTOR3 p = **pnIter; 
      //printf(" %f %f %f ", p[0], p[1], p[2]); 
      glVertex3f(p[0], p[1], p[2]); 
    }
    glEnd(); 
  }
  glPopMatrix(); 
}

void animate_streamlines() {

  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  static std::list<VECTOR3*>::iterator *pnIter; 
  static int frame = 0; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]);

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
	// ADD-BY-LEETEN 12/20/2011-BEGIN
	#ifdef	WIN32
	Sleep(500); 
	#else
	// ADD-BY-LEETEN 12/20/2011-END
	sleep(.5); 
	// ADD-BY-LEETEN 12/20/2011-BEGIN
	#endif	// #ifdef	WIN32
	// ADD-BY-LEETEN 12/20/2011-END
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
    draw_streamlines(); 
  else if (toggle_animate_streamlines == true)
    animate_streamlines(); 

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
  //printf("call idle....\n");
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

int main(int argc, char** argv) 
{

	// VTK
	// Start by loading some data.
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	// set data
	{
		char path[256];
		sprintf(path, "%s/curvilinear/combxyz.bin", SAMPLE_DATA_DIR);
		printf("%s\n", path);
		pl3dReader->SetXYZFileName(path);
		sprintf(path, "%s/curvilinear/combq.bin", SAMPLE_DATA_DIR);
		pl3dReader->SetQFileName(path);
	}
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();


	// random points
	//vtkStructuredGrid *grid = pl3dReader->GetOutput();
	//int *dim = grid->GetDimensions();
	vtkSmartPointer<vtkDataSet> sData = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));
	double *bounds = sData->GetBounds();
	printf("bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);


	osuflow = new OSUFlow;
	CVectorField *field = new VectorFieldVTK( sData );
	osuflow->SetFlowField( field );

	// openmp
#ifdef _OPENMP
	osuflow->initOpenMP(1);
#endif

	// gen seeds
	minLen[0] = bounds[0]; maxLen[0] = bounds[1];
	minLen[1] = bounds[2]; maxLen[1] = bounds[3];
	minLen[2] = bounds[4]; maxLen[2] = bounds[5];

	center[0] = (minLen[0]+maxLen[0])/2.0;
	center[1] = (minLen[1]+maxLen[1])/2.0;
	center[2] = (minLen[2]+maxLen[2])/2.0;
	printf("center is at %f %f %f \n", center[0], center[1], center[2]);
	len[0] = maxLen[0]-minLen[0];
	len[1] = maxLen[1]-minLen[1];
	len[2] = maxLen[2]-minLen[2];

#if 0

	// read in the vector field

  VECTOR3 minB, maxB; 

  osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 

  osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 

  osuflow->Boundary(minLen, maxLen); // get the boundary 
  minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
  maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
  //  osuflow->SetBoundary(minB, maxB);  // set the boundary. just to test
                                         // the subsetting feature of OSUFlow
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                minLen[0], maxLen[0], minLen[1], maxLen[1], 
                                minLen[2], maxLen[2]); 

  center[0] = (minLen[0]+maxLen[0])/2.0; 
  center[1] = (minLen[1]+maxLen[1])/2.0; 
  center[2] = (minLen[2]+maxLen[2])/2.0; 
  printf("center is at %f %f %f \n", center[0], center[1], center[2]); 
  len[0] = maxLen[0]-minLen[0]; 
  len[1] = maxLen[1]-minLen[1]; 
  len[2] = maxLen[2]-minLen[2]; 
#endif
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  glutCreateWindow("Display streamlines"); 
  glutDisplayFunc(display); 
  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}


