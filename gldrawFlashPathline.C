
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

bool toggle_draw_pathlines = false; 
bool toggle_animate_pathlines = false; 

//int num_timesteps; 
int num_frames; 
int start_time, end_time; 
int current_frame = 0; 
float time_incr; 

float center[3], len[3]; 
int first_frame = 1; 

int draw_level=0; 

bool to_draw_bbx = false; 
bool to_draw_center = false; 
bool to_draw_vector = false; 

int current_time = 0; 

///////////////////////////////
//
// AMR book-keeping variables 
//
int num_timesteps;     // number of time steps 

float domain_lengths[3];      // x y z lengths of the entire domain 

float domain_center[3];       // center of the entire domain 


//***********************************************
//    Per time step information 
//

TimeVaryingFlashAMR *tamr; 

int    total_levels;      // how many levels in total across all time steps 

//////////////////////////////////////////////////////
//  Lattice data structure built on top of the grid 

LatticeAMR * latticeAMR; 

OSUFlow **osuflow_list; 
volume_bounds_type_f* vb_list; 

list<vtListTimeSeedTrace*> *sl_list; 
VECTOR4 **osuflow_seeds; 
int *osuflow_num_seeds; 


int total_seeds = 50000; 
int npart; 

///////////////////////////////////////////////////////////////////////
//
//
void compute_pathlines() {

  float from[3], to[3]; 

  for (int i=0; i<npart; i++) {

    VECTOR3 minB, maxB; 
    VECTOR3 * seeds; 
    int num; 

    minB[0] = vb_list[i].xmin;  
    minB[1] = vb_list[i].ymin;     
    minB[2] = vb_list[i].zmin; 
    maxB[0] = vb_list[i].xmax;  
    maxB[1] = vb_list[i].ymax;     
    maxB[2] = vb_list[i].zmax; 

    from[0] = minB[0]; to[0] = maxB[0];
    from[1] = minB[1]; to[1] = maxB[1];
    from[2] = minB[2]; to[2] = maxB[2]; 

    //    osuflow_list[i]->SetRandomSeedPoints(from, to, total_seeds/npart); 
    osuflow_list[i]->SetRandomSeedPoints(from, to, 1); 
    seeds = osuflow_list[i]->GetSeeds(num); // only VECTOR3s are generated 
    osuflow_num_seeds[i] = num; 
    osuflow_seeds[i] = new VECTOR4[num];    // now copy and augment to 4D seeds 
    for (int j=0; j<num; j++)  {
      osuflow_seeds[i][j][0] = seeds[j][0]; 
      osuflow_seeds[i][j][1] = seeds[j][1]; 
      osuflow_seeds[i][j][2] = seeds[j][2]; 

      //      osuflow_seeds[i][j][3] = vb_list[i].tmin; 
      //      printf(" tmin = %d\n",vb_list[i].tmin); 
      osuflow_seeds[i][j][3] = vb_list[i].tmin + 0.5; 
      //            osuflow_seeds[i][j][3] =  0.5; 
      //      printf(" %f %f %f %f ", osuflow_seeds[i][j][0], 
      //	     osuflow_seeds[i][j][1], osuflow_seeds[i][j][2], osuflow_seeds[i][j][3]); 
    }
    sl_list[i].clear();   // clear the trace 
  }

  // Now begin to perform pathline tracing in all subdomains if seeds are in 
  bool has_seeds = true;      // initially we always have seeds
  int num_seeds_left = total_seeds; 

  while(has_seeds == true && num_seeds_left >50) {  // loop until all particles stop 

    latticeAMR->ResetSeedLists();    // clear up the lattice seed lists

     for (int np=0; np<npart; np++) {
       int i = np; 
       if (osuflow_num_seeds[i]==0) {  // domain i is already done. 
	 //	 printf("skip domain %d \n", i); 
	 continue; 
       }

       list<vtListTimeSeedTrace*> list; 
       osuflow_list[i]->SetIntegrationParams(1, 5); 
       osuflow_list[i]->GenPathLines(osuflow_seeds[i],list, FORWARD, 
				     osuflow_num_seeds[i], 50); 

       //       printf("domain %d done integrations", i); 
       printf(" %d pathlines. \n", list.size()); 

       std::list<vtListTimeSeedTrace*>::iterator pIter; 
       //------------------------------------------------
       // looping through the trace points
       pIter = list.begin(); 
       for (; pIter!=list.end(); pIter++) {
	 vtListTimeSeedTrace *trace = *pIter; 
	 sl_list[i].push_back(trace); 
       }
      // now redistributing the boundary pathline points to its neighbors. 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
	vtListTimeSeedTrace *trace = *pIter; 
	if (trace->size() ==0) continue; 
	std::list<VECTOR4*>::iterator pnIter; 
	pnIter = trace->end(); 
	pnIter--; 
	VECTOR4 p = **pnIter; 
	//check p is in which neighbor's domain 
	int ei, ej, ek, et, el; 
	int neighbor = latticeAMR->GetNeighbor(i, p[0], p[1], p[2], p[3], ei, 
					       ej, ek, et, el); 

	if (neighbor!=-1 && neighbor!=i) latticeAMR->InsertSeed(ei, ej, ek, et, el, p); 
	printf(" insert a seed %f %f %f %f to rank %d \n",
	       p[0], p[1], p[2], p[3], neighbor); 
      }
    }
    //------------------
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
      osuflow_seeds[i] = new VECTOR4[osuflow_num_seeds[i]]; 
      std::list<VECTOR4>::iterator seedIter; 
      seedIter = latticeAMR->seedlists[i].begin(); 
      int cnt = 0; 
      for (; seedIter!=latticeAMR->seedlists[i].end(); seedIter++){
	VECTOR4 p = *seedIter; 
	osuflow_seeds[i][cnt++] = p; 
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////
void draw_pathlines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
   glTranslatef(-center[0], -center[1], -center[2]); 

  for (int i=0; i<npart; i++) {

    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 
    glColor3f(1,1,1); 
    glBegin(GL_POINTS); 
    //    for (int j=0; j<nSeeds; j++) 
    //      glVertex3f(seeds[j][0], seeds[j][1], seeds[j][2]); 
    glEnd(); 

    glColor3f(1,1,0); 
    std::list<vtListTimeSeedTrace*>::iterator pIter; 

    pIter = sl_list[i].begin(); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      vtListTimeSeedTrace *trace = *pIter; 
      std::list<VECTOR4*>::iterator pnIter; 
      pnIter = trace->begin(); 
      glBegin(GL_LINE_STRIP); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR4 p = **pnIter; 
	float x = p[0]; 
	float y = p[1]; 
	float z = p[2]; 
	glVertex3f(x,y,z); 
      }
      glEnd(); 
    }
  }
  glPopMatrix(); 
}

void animate_pathlines() {
 
  std::list<vtListTimeSeedTrace*>::iterator pIter; 
  vtListTimeSeedTrace *trace; 
  std::list<VECTOR4*>::iterator pnIter; 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  float min_time = current_frame * time_incr; 
  float max_time = (current_frame+1) * time_incr; 


  for (int i=0; i<npart; i++) {
    pIter = sl_list[i].begin(); 
    glColor3f(1,1,0); 
    glBegin(GL_POINTS); 
    for (; pIter!=sl_list[i].end(); pIter++) {
      trace = *pIter; 
      pnIter = trace->begin(); 
      for (; pnIter!= trace->end(); pnIter++) {
	VECTOR4 p = **pnIter; 
	if (p[3]>= min_time && p[3] < max_time) {
	  float x = p[0]; 
	  float y = p[1]; 
	  float z = p[2]; 
	  glVertex3f(x,y,z); 
	}
      }
    }
    glEnd(); 
  }

  glPopMatrix(); 
  current_frame = (current_frame+1) % num_frames; 
}

void draw_vector(int t, int i) 
{

  FlashAMR *amr = tamr->GetTimeStep(t); 

  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glBegin(GL_LINES); 
  glColor3f(1,0,0); 

  float *center = amr->GetBlockCenter(i); 
  float *vectors = amr->GetDataPtr(i); 

  glVertex3f(center[0], center[1], center[2]); 
  glColor3f(0,0,1); 

  glVertex3f(center[0]+ vectors[0]*0.3, 
	     center[1]+ vectors[1]*0.3, 
	     center[2]+ vectors[2]*0.3); 

  glEnd(); 

  glPopMatrix(); 

}

////////////////////////////////////////////// 

void draw_cube(float r, float g, float b)
{
  glColor3f(r, g, b); 
  glutWireCube(1.0);   // draw a solid cube 
}

///////////////////////////////////////////////

void draw_bbx(int t, int i) {

  FlashAMR* amr = tamr->GetTimeStep(t); 

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

  // test the centers
  glPushMatrix(); 
  glTranslatef(center[0], center[1], center[2]); 
  glScalef(length[0], length[1], length[2]); 

  glutWireCube(1.0); 
  glPopMatrix(); 

  /*
  // test the volume bounds 

  int rank = latticeAMR->GetRank(center[0], center[1], center[2], t); 

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
  */
  glPopMatrix(); 
}

/////////////////////////////////////////////////////
void draw_center(int t, int i) 
{
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-center[0], -center[1], -center[2]); 

  glPointSize(2); 

  glColor3f(1,1,1); 

  FlashAMR *amr = tamr->GetTimeStep(t); 
  float *center = amr->GetBlockCenter(i); 

  int level = latticeAMR->GetFinestLevel(center[0], center[1], center[2], t); 

  if (level %3 == 0) glColor3f(1,0,0); 
  if (level %3 == 1) glColor3f(0,1,0); 
  if (level %3 == 2) glColor3f(0,0,1); 

  glBegin(GL_POINTS); 
  glVertex3f(center[0], center[1], center[2]); 
  glEnd(); 

  glPopMatrix(); 
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

  if (toggle_draw_pathlines == true)
    draw_pathlines(); 
  else if (toggle_animate_pathlines == true)
    animate_pathlines(); 

  FlashAMR *amr = tamr->GetTimeStep(current_time); 
  int nb = amr->GetNumBlocks(); 
  for (int i=0; i<nb; i++) {
    int blevel = amr->GetLevel(i); 
    if (blevel>=draw_level) {
      if (to_draw_bbx) draw_bbx(current_time, i); 
      if (to_draw_center) draw_center(current_time, i); 
      if (to_draw_vector) draw_vector(current_time, i); 
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
	case 's': 
	  compute_pathlines(); 
	  glutPostRedisplay(); 
	  break; 
	case 'l': 
	  draw_level++; 
	  draw_level = draw_level % (total_levels); 
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
	case 't': 
	  current_time +=1; 
	  current_time = current_time % num_timesteps; 
	  break; 
	case 'd': 
	  toggle_draw_pathlines = !toggle_draw_pathlines; 
	  toggle_animate_pathlines = false; 
	  break; 
	case'a': 
	  toggle_animate_pathlines = !toggle_animate_pathlines; 
	  toggle_draw_pathlines = false; 
	  first_frame = 1; 
	}
	glutPostRedisplay(); 
}

void timer(int val) {
  if (toggle_animate_pathlines == true) {
    glutPostRedisplay(); 
  }
  glutTimerFunc(200, timer, 0); 
}


/////////////////////////////////////////////////////////////////////

main(int argc, char* argv[]) 
{
  float min[3], max[3]; 
  int dims[3]; 

  tamr = new TimeVaryingFlashAMR; 
  tamr->LoadData(argv[1], min, max); 

  tamr->GetDims(dims); 
  num_timesteps = tamr->GetNumTimeSteps(); 
  total_levels = tamr->GetNumLevels(); 
  

  len[0] = domain_lengths[0] = max[0]-min[0]; 
  len[1] = domain_lengths[1] = max[1]-min[1]; 
  len[2] = domain_lengths[2] = max[2]-min[2]; 
  center[0] = domain_center[0] = (min[0]+max[0])/2.0; 
  center[1] = domain_center[1] = (min[1]+max[1])/2.0; 
  center[2] = domain_center[2] = (min[2]+max[2])/2.0; 

  latticeAMR = new LatticeAMR(domain_lengths[0], domain_lengths[1], domain_lengths[2], 
			      num_timesteps, total_levels); 

  for (int i=0; i<total_levels; i++)  {

    float blockSize[3]; 
    float levelMinB[3], levelMaxB[3]; 

    tamr->GetLevelBlockSize(i, blockSize); 
    tamr->GetLevelBounds(i, levelMinB, levelMaxB); 

    printf(" level size: [%f %f %f] min corner [%f %f %f] max corner [%f %f %f]\n", 
	   blockSize[0], blockSize[1], blockSize[2], 
	   levelMinB[0], levelMinB[1], levelMinB[2], 
	   levelMaxB[0], levelMaxB[1], levelMaxB[2]); 

    latticeAMR->CreateLevel(i, blockSize[0], blockSize[1], blockSize[2], 
			    dims[0], dims[1], dims[2], 
			    levelMinB[0], levelMaxB[0], levelMinB[1], levelMaxB[1], 
			    levelMinB[2], levelMaxB[2], 0, num_timesteps-1); 
  }

  for (int t=0; t<num_timesteps; t++) {
    FlashAMR * amr = tamr->GetTimeStep(t); 
    int nb = amr->GetNumBlocks(); 

    for (int i=0; i<nb; i++) {

      float *center = amr->GetBlockCenter(i); 

      latticeAMR->CheckIn(amr->GetLevel(i), center[0], center[1], center[2], t, 
			  amr->GetDataPtr(i)); 
    }
  }
  latticeAMR->CompleteLevels(2); // this call is very important. It finishes up all the 
                                 // book-keeping and complete the lattice setup
  latticeAMR->InitSeedLists(); 

  vb_list = latticeAMR->GetBoundsList(npart); 

  // now initialize one osuflow object for each partition 

  osuflow_list = new OSUFlow*[npart]; 

  sl_list = new list<vtListTimeSeedTrace*>[npart]; // pathline lists,one per partition

  osuflow_seeds = new VECTOR4*[npart]; 
  osuflow_num_seeds = new int[npart]; 

  for(int i=0; i<npart; i++) {

    osuflow_list[i] = new OSUFlow(); 
    float **data = latticeAMR->GetDataPtr(i); 
    if (data == NULL) {
      printf(" panic\n"); 
      exit(1); 
    }
    float minB[3], maxB[3]; 
    int min_t, max_t; 
    minB[0] = vb_list[i].xmin; maxB[0] = vb_list[i].xmax;     
    minB[1] = vb_list[i].ymin; maxB[1] = vb_list[i].ymax;
    minB[2] = vb_list[i].zmin; maxB[2] = vb_list[i].zmax;
    min_t = vb_list[i].tmin; max_t = vb_list[i].tmax; 

    // dims[0/1/2] are the x y z dat dimensions 
    // minB and maxB are the physical space min/max corners 

    osuflow_list[i]->CreateTimeVaryingFlowField(data, dims[0], dims[1], 
						dims[2], minB, maxB, 
						min_t, max_t);     
    osuflow_list[i]->ScaleField(0.1); 
  }

  // set up the animation frame time 
  num_frames = 5;     // number of frames to loop through the data time 
  time_incr = num_timesteps/(float) num_frames; 

  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 
  
  glutCreateWindow("Display"); 
  glutDisplayFunc(display); 

  glutTimerFunc(10, timer, 0); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 
  glutMainLoop(); 
}
