
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
#include "LatticeAMR.h" 
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

///////////////////////////////
//
// AMR book-keeping variables 
//

FlashAMR *amr; 

int nb;                // number of blocks
int num_levels;        // the total number of level 

float len[3];          // global physical x y z lengths
float center[3];       // center of the domain 

LatticeAMR * latticeAMR; 
volume_bounds_type_f* vb_list; 

OSUFlow **osuflow_list; 

list<vtListSeedTrace*> *sl_list; 
VECTOR3 **osuflow_seeds; 
int *osuflow_num_seeds; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 

int num_frames = 20000; 
int current_frame = 0; 

int draw_level=0; 

bool to_draw_bbx = false; 
bool to_draw_center = false; 
bool to_draw_vector = false; 

int total_seeds = 5000; 
int npart; 
int first_frame = 1; 
////////////////////////////////////////////////////////

void compute_streamlines() {
  
  float from[3], to[3]; 

  for (int i=0; i<npart; i++) {
    VECTOR3 minB, maxB; 
    from[0] = vb_list[i].xmin;  
    from[1] = vb_list[i].ymin;     
    from[2] = vb_list[i].zmin; 
    to[0] = vb_list[i].xmax;  
    to[1] = vb_list[i].ymax;     
    to[2] = vb_list[i].zmax; 

    osuflow_list[i]->SetRandomSeedPoints(from, to, 10); // set range for seed locations
    int num; 
    osuflow_seeds[i] = osuflow_list[i]->GetSeeds(num); 
    osuflow_num_seeds[i] = num; 
    sl_list[i].clear(); 
  }
  
  // Now begin to perform particle tracing in all subdomains
  bool has_seeds = true;      // initially we always have seeds
  int num_seeds_left = 20*npart; 

  int iter_num = 0; 
  //  while(has_seeds == true && num_seeds_left >0) {  // loop until all particles stop 
  while(iter_num < 10) {  // loop until all particles stop 

    printf(" ********* itern = %d \n", iter_num); 
    iter_num++; 
    latticeAMR->ResetSeedLists();    // clear up the lattice seed lists

    for (int i=0; i<npart; i++) {
      if (osuflow_num_seeds[i]==0) {  // npart is already done. 
	printf("skip domain %d \n", i); 
	continue; 
      }
      list<vtListSeedTrace*> list; 
      osuflow_list[i]->SetIntegrationParams(1, 5); 
            osuflow_list[i]->GenStreamLines(osuflow_seeds[i], FORWARD_DIR, 
      //      osuflow_list[i]->GenStreamLines(osuflow_seeds[i], BACKWARD_AND_FORWARD, 
      				      osuflow_num_seeds[i], 50, list); 
      printf("domain %d done integrations", i); 
      printf(" %d streamlines. \n", list.size()); 

      std::list<vtListSeedTrace*>::iterator pIter; 

      //looping through and store the trace points into sl_list[i] 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
        vtListSeedTrace *trace = *pIter; 
	sl_list[i].push_back(trace); 
      }
      // now redistributing the boundary streamline points to its neighbors. 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
	vtListSeedTrace *trace = *pIter; 
	if (trace->size() ==0) continue; 
	std::list<VECTOR3*>::iterator pnIter; 
	pnIter = trace->end(); 
	pnIter--; 
	VECTOR3 p3 = **pnIter; 
	VECTOR4 p(p3[0],p3[1],p3[2],0.0);

	//check p is in which neighbor's domain 
	int ei, ej, ek, et, el; 
	int neighbor = latticeAMR->GetNeighbor(i, p[0], p[1], p[2], p[3], ei, ej, ek, et, el); 
	if (neighbor!=-1 && neighbor!=i) latticeAMR->InsertSeed(ei, ej, ek, et, el, p); 
	//	printf(" insert a seed %f %f %f %f to [%d %d %d %d %d] rank %d \n",
	//	       p[0], p[1], p[2], p[3], ei, ej,ek,et,el,neighbor); 
	//	if (neighbor!=-1) printf("seedlists[%d].size() = %d\n", neighbor, 
	//     latticeAMR->seedlists[neighbor].size()); 
      }
    }
    // now create the seed arrays for the next run
    has_seeds = false;  
    num_seeds_left = 0; 
    for (int i=0; i<npart; i++) {
      //      if (osuflow_seeds[i]!=NULL) delete [] osuflow_seeds[i]; 
      osuflow_num_seeds[i] = latticeAMR->seedlists[i].size(); 
      num_seeds_left += osuflow_num_seeds[i]; 
      //      printf("seedlists[%d].size() = %d\n", i, osuflow_num_seeds[i]); 
      if (osuflow_num_seeds[i]!=0) has_seeds = true; 
      else continue; 
      osuflow_seeds[i] = new VECTOR3[osuflow_num_seeds[i]]; 
      std::list<VECTOR4>::iterator seedIter; 
      seedIter = latticeAMR->seedlists[i].begin(); 
      int cnt = 0; 
      for (; seedIter!=latticeAMR->seedlists[i].end(); seedIter++){
	VECTOR4 p4 = *seedIter; 
	VECTOR3 p(p4[0],p4[1],p4[2]); 
	osuflow_seeds[i][cnt++] = p; 
      }
    }
  }
}



void draw_streamlines() {

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  std::list<vtListSeedTrace*>::iterator pIter; 

  for (int i=0; i<npart; i++) {   // looping through all subdomains 

    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 

    /*
    glPointSize(5.0); 
    glColor3f(1,1,1); 
    glBegin(GL_POINTS); 
    for (int j=0; j<nSeeds; j++) {
      //      printf(" Seed: %f %f %f \n", seeds[j][0], seeds[j][1], seeds[j][2]); 
      glVertex3f(seeds[j][0], seeds[j][1], seeds[j][2]); 
    }
    glEnd(); 
    glPointSize(1.0); 
    */

    pIter = sl_list[i].begin(); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      vtListSeedTrace *trace = *pIter; 
      std::list<VECTOR3*>::iterator pnIter; 
      pnIter = trace->begin(); 
      glColor3f(1,1,0); 
      glBegin(GL_LINE_STRIP); 
	  //      glBegin(GL_POINTS); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
      //	printf(" %f %f %f \n", p[0], p[1], p[2]); 
	glVertex3f(p[0], p[1], p[2]); 
      }
      glEnd(); 
    }
  }
  glPopMatrix(); 
}

void animate_streamlines() {

  std::list<vtListSeedTrace*>::iterator pIter; 
  vtListSeedTrace *trace; 
  std::list<VECTOR3*>::iterator pnIter; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  float min_time = current_frame/(float) num_frames; 
  float max_time = (current_frame+1)/(float) num_frames; 

  glColor3f(1,1,0); 
  for (int i=0; i<npart; i++) {
    pIter = sl_list[i].begin(); 

    for (; pIter!=sl_list[i].end(); pIter++) {
      trace = *pIter; 
      pnIter = trace->begin(); 
      int cnt=0; 
      glBegin(GL_LINE_STRIP); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR3 p = **pnIter; 
	float ratio = cnt/(float)trace->size(); 
	//	if (ratio>= min_time && ratio < max_time) {
	if (cnt == current_frame % trace->size()) 
	  glColor3f(1,1,0); 
	else 
	  glColor3f(0.1,0.1,0.1); 

	float x = p[0]; 
	float y = p[1]; 
	float z = p[2]; 
	glVertex3f(x,y,z); 
	cnt++; 
      }
      glEnd(); 
    }
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

  int rank = latticeAMR->GetRank(center[0], center[1], center[2], 0);  

  if (rank !=-1) {
  float xmin = vb_list[rank].xmin; 
  float ymin = vb_list[rank].ymin; 
  float zmin = vb_list[rank].zmin; 
  float xmax = vb_list[rank].xmax; 
  float ymax = vb_list[rank].ymax; 
  float zmax = vb_list[rank].zmax; 
  
  glBegin(GL_LINE_STRIP); 
  glVertex3f(xmin, ymin, zmin); 
  glVertex3f(xmax, ymin, zmin); 
  glVertex3f(xmax, ymin, zmax); 
  glVertex3f(xmin, ymin, zmax); 
  glVertex3f(xmin, ymin, zmin); 
  glEnd(); 
    }
  glPopMatrix(); 

}

void draw_vector(int i) 
{
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glBegin(GL_LINES); 
  glColor3f(1,0,0); 


  float *center = amr->GetBlockCenter(i); 
  glVertex3f(center[0], center[1], center[2]); 
  float* vectors = amr->GetDataPtr(i); 
  glColor3f(0,0,1); 
  glVertex3f(center[0]+ vectors[0]*0.3, 
	     center[1]+ vectors[1]*0.3, 
	     center[2]+ vectors[2]*0.3); 
  glEnd(); 

  glPopMatrix(); 

}

/////////////////////////////////////////////////////
void draw_center(int i) 
{
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glPointSize(2); 

  glColor3f(1,1,1); 

  float *center = amr->GetBlockCenter(i); 

  int level = latticeAMR->GetFinestLevel(center[0], center[1], center[2],0);

  if (level %3 == 0) glColor3f(1,0,0); 
  if (level %3 == 1) glColor3f(0,1,0); 
  if (level %3 == 2) glColor3f(0,0,1); 

  glBegin(GL_POINTS); 
  glVertex3f(center[0], center[1], center[2]); 
  glEnd(); 

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

  if (toggle_draw_streamlines == true)
    draw_streamlines(); 
  else if (toggle_animate_streamlines == true)
    animate_streamlines(); 

  /*
  if (to_draw_bbx) 
    for (int i=0; i<npart; i++)
      draw_bbx(i); 
  */

  for (int i=0; i<nb; i++) {
    int its_level = amr->GetLevel(i); 
    if (its_level>=draw_level) {
      if (to_draw_bbx) draw_bbx(i); 
      if (to_draw_center) draw_center(i); 
      if (to_draw_vector) draw_vector(i); 
    }
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

void timer(int val) {
  if (toggle_animate_streamlines == true) {
    glutPostRedisplay(); 
  }
  glutTimerFunc(200, timer, 0); 
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
	case 'l': 
	  draw_level++; 
	  draw_level = draw_level % (num_levels); 
	  break; 
	case 'b':
	  to_draw_bbx = !to_draw_bbx; 
	  break; 
	case 'c': 
	  to_draw_center = !to_draw_center; 
	  break; 
	case 'v': 
	  to_draw_vector = !to_draw_vector; 
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
	glutPostRedisplay(); 
}


/////////////////////////////////////////////////////////////////////

main(int argc, char* argv[]) 
{

  float min[3], max[3]; 

  amr = new FlashAMR; 
  amr->LoadData(argv[1], min, max); 
  nb = amr->GetNumBlocks(); 
  num_levels = amr->GetNumLevels(); 

  int dims[3]; 
  amr->GetDims(dims); 

  len[0] = max[0]-min[0]; 
  len[1] = max[1]-min[1]; 
  len[2] = max[2]-min[2]; 
  center[0] = (min[0]+max[0])/2.0; 
  center[1] = (min[1]+max[1])/2.0; 
  center[2] = (min[2]+max[2])/2.0; 

  latticeAMR = new LatticeAMR(len[0], len[1], len[2], 1, num_levels); 

  for (int i=0; i<num_levels; i++)  {


    float blockSize[3]; 
    float levelMinB[3], levelMaxB[3]; 

    amr->GetLevelBlockSize(i, blockSize); 
    amr->GetLevelBounds(i, levelMinB, levelMaxB); 

    printf(" level size: [%f %f %f] min corner [%f %f %f] max corner [%f %f %f]\n", 
	   blockSize[0], blockSize[1], blockSize[2], 
	   levelMinB[0], levelMinB[1], levelMinB[2], 
	   levelMaxB[0], levelMaxB[1], levelMaxB[2]); 

    latticeAMR->CreateLevel(i, blockSize[0], blockSize[1], blockSize[2], 
			    dims[0], dims[1], dims[2], 
			    levelMinB[0], levelMaxB[0], levelMinB[1], levelMaxB[1], 
			    levelMinB[2], levelMaxB[2], 0, 0); 

  }
  for (int i=0; i<nb; i++) {

    float* center = amr->GetBlockCenter(i); 
    latticeAMR->CheckIn(amr->GetLevel(i), center[0], center[1], center[2], 0, 
			amr->GetDataPtr(i)); 
  }

  latticeAMR->MergeBlocks();   // this is what's new in gldrawFlash4: 
                               // Merge adjacent blocks that have the 
                               // same resolutions together into larger blocks 
                               // to minimize overhead 

  latticeAMR->MergeAndCompleteLevels(); 

  latticeAMR->InitSeedLists(); 

  vb_list = latticeAMR->GetBoundsList(npart); 
  printf("nb = %d npart = %d \n", nb, npart);  // they should be equal 

  sl_list = new list<vtListSeedTrace*>[npart]; //one streamlines list for each subdomain 

  // now initialize the osuflow objects, one per block 
  osuflow_list = new OSUFlow*[npart]; 
  osuflow_seeds = new VECTOR3*[npart]; 
  osuflow_num_seeds = new int[npart]; 

  for (int i=0; i<npart; i++) {
    osuflow_list[i] = new OSUFlow(); 
    float **data = latticeAMR->GetDataPtr(i); // i is the rank 
    if (data == NULL) {
      printf(" panic\n"); 
      exit(1); 
    }
    float minB[3], maxB[3]; 
    minB[0] = vb_list[i].xmin; maxB[0] = vb_list[i].xmax;     
    minB[1] = vb_list[i].ymin; maxB[1] = vb_list[i].ymax;
    minB[2] = vb_list[i].zmin; maxB[2] = vb_list[i].zmax;
    // dims[0/1/2] are the x y z dat dimensions 
    // minB and maxB are the physical space min/max corners 
    printf(" block %d =  dim %d %d %d ....\n", i,  vb_list[i].xdim, vb_list[i].ydim, vb_list[i].zdim); 

    osuflow_list[i]->CreateStaticFlowField(data[0], vb_list[i].xdim, vb_list[i].ydim, vb_list[i].zdim,
    					   minB, maxB); 
  }

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
