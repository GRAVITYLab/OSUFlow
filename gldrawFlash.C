
#include <stdio.h>
#include <stdlib.h> 

#ifdef MAC_OSX
#include <GLUT/glut.h> 
#include <OpenGL/gl.h>
#else
#include <GL/glut.h> 
#include <GL/gl.h>
#endif 


int press_x, press_y; 
int release_x, release_y; 
float x_angle = 0.0; 
float y_angle = 0.0; 
float scale_size = 1; 

int xform_mode = 0; 

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 


int nb; 
int dims[3]; 
float *mins, *maxs;  
float *centers;   // block center 
float *lengths;   // block dimensions 
int   *indices;   // block index 
float *vectors;   //data 
float len[3];     // domain dimensions 
float center[3];  // domain center 
float *element_xsize;  //voxel x size in each block 
float *element_ysize;  //voxel y size
float *element_zsize;  //voxel z size 
int *its_level;        //the level that each block is in 
int num_levels;        //total number of levels 

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
  //  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 


  float xsize, ysize, zsize; 
  if (its_level[i]%3 == 0) glColor3f(1,0,0); 
  if (its_level[i]%3 == 1) glColor3f(0,1,0); 
  if (its_level[i]%3 == 2) glColor3f(0,0,1); 

  glPushMatrix(); 
  glTranslatef(centers[i*3], centers[i*3+1], centers[i*3+2]); 
  glScalef(lengths[i*3], lengths[i*3+1], lengths[i*3+2]); 

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

  
  //  printf(" len %f %f %f\n", len[0], len[1], len[2]); 
  glPushMatrix(); 
  //  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  glScalef(1.0, len[1]/len[0], len[2]/len[0]); 
  draw_cube(0,0,1);
  glPopMatrix(); 

  for (int i=0; i<nb; i++) {
      if (its_level[i]>=draw_level) 
      draw_bbx(i); 
  }
    
  //  glColor3f(1,1,1); 
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

void sort_list(float* list, int& cnt) 
{
  float* tmp_list = new float[cnt]; 
  for (int i=0; i<cnt; i++) tmp_list[i] = list[i]; 
  for (int i=0; i<cnt; i++) {
    int largest = -1; 
    int idx = i; 
    for (int j = i; j<cnt; j++) {
      if (tmp_list[j]>largest) {idx = j; largest = tmp_list[j]; }
    }
    int tmp = tmp_list[i]; 
    tmp_list[i] = tmp_list[idx]; 
    tmp_list[idx] =tmp; 
  }
  list[0] = tmp_list[cnt-1]; 
  int count = 1; 
  for (int i=cnt-2; i>=0; i--) {
    if (tmp_list[i]>list[count-1]) {
      list[count] = tmp_list[i]; 
      count++; 
    }
  }
  cnt = count; 
}


main(int argc, char* argv[]) {

  float *sizes; 
  FILE * in = fopen(argv[1], "r"); 
  float min[3], max[3]; 
  
  fread(&nb, sizeof(int), 1, in); 
  printf(" %d blocks.\n", nb); 

  fread(dims, sizeof(int), 3, in); 
  printf(" block dims %d %d %d\n", dims[0], dims[1], dims[2]); 

  int size = dims[0]*dims[1]*dims[2];  // number of elements 

  its_level = new int[nb]; 
  indices = new int[nb];      // nb of block index   
  centers = new float[nb*3];  // nb of center x y and z 
  lengths = new float[nb*3];  // nb of lengths in x y and z 
  mins = new float[nb*3];     // nb of min_x, min_y, min_z
  maxs = new float[nb*3];     // nb of max_x, max_y, max_z
  vectors = new float[size*3];   //block vector data 
  element_xsize = new float[nb]; 
  element_ysize = new float[nb]; 
  element_zsize = new float[nb]; 
  sizes = new float[nb]; 

  for (int i=0; i<nb; i++) {

    int idx = i*3; 
    fread(&(indices[i]), sizeof(int), 1, in); 
    printf(" read block id: %d\n", indices[i]); 
    fread(&(centers[idx]), sizeof(float), 3, in); 
    printf(" block %d center: [%f %f %f]\n", indices[i], centers[idx], centers[idx+1], centers[idx+2]); 

    fread(&(lengths[idx]), sizeof(float), 3, in); 
    printf(" block %d length: [%f %f %f]\n", indices[i], lengths[idx], lengths[idx+1], lengths[idx+2]); 

    float bounds[6]; 
    fread(bounds, sizeof(float), 6, in); 
    mins[idx]   = bounds[0];  maxs[idx]   = bounds[1]; 
    mins[idx+1] = bounds[2];  maxs[idx+1] = bounds[3]; 
    mins[idx+2] = bounds[4];  maxs[idx+2] = bounds[5]; 

    if (i==0) {
      min[0] = mins[0]; min[1] = mins[1]; min[2] = mins[2]; 
      max[0] = maxs[0]; max[1] = maxs[1]; max[2] = maxs[2]; 
    }
    else {
      if (mins[idx]<min[0]) min[0] = mins[idx]; 
      if (mins[idx+1]<min[1]) min[1] = mins[idx+1];      
      if (mins[idx+2]<min[2]) min[2] = mins[idx+2]; 
      if (maxs[idx]>max[0]) max[0] = maxs[idx]; 
      if (maxs[idx+1]>max[1]) max[1] = maxs[idx+1]; 
      if (maxs[idx+2]>max[2]) max[2] = maxs[idx+2]; 
    }
    //    printf(" block range: [%f %f %f]-[%f %f %f]\n", mins[i*3], mins[i*3+1], mins[i*3+2], 
    //	   maxs[i*3], maxs[i*3+1], maxs[i*3+2]); 

    float *xcomp = new float[size]; 
    float *ycomp = new float[size]; 
    float *zcomp = new float[size]; 

    fread(xcomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(ycomp, sizeof(float), size, in);  // read the data (vectors) 
    fread(zcomp, sizeof(float), size, in);  // read the data (vectors) 

    for (int c=0; c<size; c++) {
      vectors[c*3]   =xcomp[c]; 
      vectors[c*3+1] =ycomp[c]; 
      vectors[c*3+2] =zcomp[c]; 
    }

    delete [] xcomp; 
    delete [] ycomp; 
    delete [] zcomp; 

    printf("first vector component %f \n", vectors[1728]); 
    sizes[i] = element_xsize[i] = lengths[idx]/(float)dims[0]; 
    element_ysize[i] = lengths[idx+1]/(float)dims[1]; 
    element_zsize[i] = lengths[idx+2]/(float)dims[2]; 
  }

  int count = nb; 
  sort_list(sizes, count); 
  printf(" there are %d different sizes\n", count); 
  num_levels = count; 

  for (int i=0; i<nb; i++) {
    for (int j=0; j<count; j++) {
      if (element_xsize[i] == sizes[j]) its_level[i]= num_levels-1 - j; 
    }
  }

  for (int i=0; i<nb; i++) 
    printf(" bloak %d is at level %d \n", i, its_level[i]); 

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
    //  glutIdleFunc(idle); 

    //    glutTimerFunc(10, timer, 0); 
    glutMouseFunc(mymouse); 
    glutMotionFunc(mymotion);
    glutKeyboardFunc(mykey); 
    glutMainLoop(); 
}
