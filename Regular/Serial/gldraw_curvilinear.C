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
list<vtListSeedTrace*> sl_list; 

bool toggle_draw_streamlines = false; 
bool toggle_animate_streamlines = false; 
float center[3], len[3]; 
int first_frame = 1; 

//draw grid
int dim[3];
float * grid_fx, *grid_fy, *grid_fz;
list<vtListTimeSeedTrace*> sl_list_pathlines;
bool steady=true;
char flowname[250];
////////////////////////////////////////////////////////
void get_seed(int* ind, float* seed)
{
	int i,j,k;
	i=ind[0]; j=ind[1]; k=ind[2];
	int id=i+j*dim[0]+k*dim[0]*dim[1];

	seed[0]=grid_fx[id];
	seed[1]=grid_fy[id];
	seed[2]=grid_fz[id]; 


}
void compute_streamlines() {

	float from[3], to[3]; 

	from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
	to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

	from[0] = 0;   from[1] = 0;   from[2] = 0; 
	to[0] = 1;   to[1] = 1;   to[2] = 1; 

	printf("generating seeds...\n"); 

	//generate random seeds
	int seed_num=dim[0]*dim[1];
	VECTOR3* seed_list=new VECTOR3[seed_num];
	int ind[3];
	float seed[3];
	int count=0;
	int z=15;
//	for(int z=0; z<dim[2];z+=10)
	{
		for(int y=0; y<dim[1];y+=5)
		{
			for(int x=0; x<dim[0];x+=5)
			{
				ind[0]=x; ind[1]=y; ind[2]=z;

				get_seed(ind, seed);
				seed_list[count++].Set(seed[0],seed[1],seed[2]);

			}
		}
	}
		

	seed_num=count;
	osuflow->SetSeedPoints(seed_list, seed_num);
	// osuflow->SetRandomSeedPoints(from, to, 10); 
	int nSeeds; 
	VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
	//for (int i=0; i<nSeeds; i++) \
		printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], \
		seeds[i][1], seeds[i][2]); 

	sl_list.clear(); 

	osuflow->SetIntegrationParams(.1, .5); 
	osuflow->GenStreamLines(sl_list , FORWARD_DIR, 500, 0); 
	if(steady==true)
	  {
	printf("compute streamlines..\n"); 
	osuflow->SetIntegrationParams(.1, .5); 
	osuflow->GenStreamLines(sl_list , FORWARD_DIR, 500, 0); 
	  }
	else
	  {
	    printf("compute path lines\n");
	    osuflow->SetIntegrationParams(.1, .5);
	    int num_timesteps;
	    FILE* fp=fopen(flowname,"r");
	    fscanf(fp,"%d",&num_timesteps);
	    fclose(fp);
	    printf("time step=%d\n",num_timesteps);
	    float* tarray = new float[nSeeds]; 
	    for (int i=0;i<nSeeds; i++) 
	      tarray[i] = (float)(i % num_timesteps); 
	    osuflow->GenPathLines(seeds, sl_list_pathlines , FORWARD, nSeeds, 500, tarray); 
  	    printf(" done integrations\n");  
	    delete [] tarray;	    
	  }
	printf(" done integrations\n"); 
	printf("list size = %d\n", sl_list.size()); 
	/*
	std::list<vtListSeedTrace*>::iterator pIter; 
	pIter = sl_list.begin(); 
	for (; pIter!=sl_list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter; 
		printf("line len=%d\n",trace->size());
	}*/
}

void draw_grid()
{
	glPushMatrix(); 
	glColor3f(1,1,1); 
	glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
	glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 
		
	//for(int k=0; k<dim[2]-1;k+=dim[2]-2)
	//int k=0;
	for(int k=0; k<dim[2]-1;k+=10)  
	{
		for(int j=0; j<dim[1]-1;j++)  
		{
			for(int i=0; i<dim[0]-1;i++)
			{
				glBegin(GL_LINE_STRIP);
				//glBegin(GL_QUADS);
				int id[4];
				id[0]=i+j*dim[0]+k*dim[0]*dim[1];
				id[1]=i+1+j*dim[0]+k*dim[0]*dim[1];
				id[3]=i+(j+1)*dim[0]+k*dim[0]*dim[1];
				id[2]=i+1+(j+1)*dim[0]+k*dim[0]*dim[1];

				for(int l=0;l<4; l++)
				{
					VECTOR3 v1,v2,v3,v;
					v1.Set(grid_fx[id[l]], grid_fy[id[l]], grid_fz[id[l]]);
					v2.Set(grid_fx[id[(l+1)%4]], grid_fy[id[(l+1)%4]], grid_fz[id[(l+1)%4]]);
					if(l-1>=0)
						v3.Set(grid_fx[id[(l-1)]], grid_fy[id[(l-1)]], grid_fz[id[(l-1)]]);
					else
						v3.Set(grid_fx[id[3]], grid_fy[id[3]], grid_fz[id[3]]);
					v=cross(v3-v1,v2-v1);
					v.Normalize(); 
					
					glNormal3f(v.x(),v.y(),-v.z());
					glVertex3f(grid_fx[id[l]], grid_fy[id[l]], grid_fz[id[l]]); 
				}

				glEnd();
			}
		}
	}

	glPopMatrix(); 
}


void draw_streamlines() {

	glPushMatrix(); 

	glLineWidth(2);
	glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
	glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

	printf("draw streamlines here.\n"); 
	glColor3f(1,0,0); 
	std::list<vtListSeedTrace*>::iterator pIter; 
	if(steady==true)
	{
		pIter = sl_list.begin(); 
		for (; pIter!=sl_list.end(); pIter++) {
			vtListSeedTrace *trace = *pIter; 
			std::list<VECTOR3*>::iterator pnIter; 
			pnIter = trace->begin(); 
			if(trace->size()<=2)
				continue;
			glBegin(GL_LINE_STRIP); 
			for (; pnIter!= trace->end(); pnIter++) {
				VECTOR3 p = **pnIter; 
				//printf(" %f %f %f ", p[0], p[1], p[2]); 
				glVertex3f(p[0], p[1], p[2]); 
			}
			glEnd(); 
		}
	}
	glLineWidth(1);

	glPopMatrix(); 
}

void draw_pathlines() {
  
  glPushMatrix(); 

  glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
  glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

  printf("draw pathlines.\n"); 
  glColor3f(1,1,0); 
  std::list<vtListTimeSeedTrace*>::iterator pIter; 

  pIter = sl_list_pathlines.begin(); 
  for (; pIter!=sl_list_pathlines.end(); pIter++) {
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
 // sleep(.5); 
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

	GLfloat mat_specular[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat mat_diffuse[] = {0.3, 0.3, 0.3, 1.0};	
	GLfloat mat_ambient[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat mat_shininess[] = {5};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    
	GLfloat light_specular[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat light_diffuse[] = {0.5, 0.5, 0.5, 1.0};
	GLfloat light_ambient[] = {0.1, 0.1, 0.3, 1.0};	
	GLfloat light_position[] = {0.0, 0.0, 10.0, 0.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);    

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);		

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	//if (toggle_draw_streamlines == true)
	glPushMatrix(); 

	glRotatef(x_angle, 0, 1,0); 
	glRotatef(y_angle, 1,0,0); 
	glScalef(scale_size, scale_size, scale_size); 
	//draw_grid();
	glDisable(GL_LIGHTING);

	if(steady==true)
	draw_streamlines(); 
	else
	draw_pathlines();
	
	//if (toggle_animate_streamlines == true)
	 //animate_streamlines(); 
	//else
	//	draw_streamlines(); 

	glPopMatrix(); 

	printf(" len %f %f %f\n", len[0], len[1], len[2]); 
	
	glutSwapBuffers(); 
}

void timer(int val) {
	//printf("call idle....\n"); 
	if (toggle_animate_streamlines == true) {
		  animate_streamlines(); 
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


void readGrid(char* filename)
{
	char gird_file[255];
	memset(gird_file,0,255);
	sprintf(gird_file,"%s.grid",filename);
	printf("grid file=%s\n",gird_file);
	FILE* fIn = fopen(gird_file, "rb");
	if(fIn==0)
	  {
	    printf("can not find grid file, abort\n");
	    exit(0);
	  }

	fread(&dim[0], sizeof(int), 1, fIn);
	fread(&dim[1], sizeof(int), 1, fIn);
	fread(&dim[2], sizeof(int), 1, fIn);

	int nodeNum=dim[0]*dim[1]*dim[2];

	grid_fx = (float *)malloc(sizeof(float) * nodeNum);
	grid_fy = (float *)malloc(sizeof(float) * nodeNum);
	grid_fz = (float *)malloc(sizeof(float) * nodeNum);


	fread(grid_fx, sizeof(float), nodeNum, fIn);					// X for nodes

	fread(grid_fy, sizeof(float), nodeNum, fIn);					// Y for nodes

	fread(grid_fz, sizeof(float), nodeNum, fIn);					// Z for nodes
	
	fclose(fIn);
}

int main(int argc, char** argv) 
{
	VECTOR3 minB, maxB; 
	readGrid(argv[1]);
	osuflow = new OSUFlow(); 
	printf("read file %s\n", argv[1]); 
	if(argc>2)
	{
		int isteady=atoi(argv[2]);
		steady=(isteady==0);
	}
	
	memset(flowname,0,250);
	sprintf(flowname,"%s",argv[1]);
	minB[0] = 0; minB[1] = 0; minB[2] = 0; 
	maxB[0] = 20; maxB[1] = 20; maxB[2] = 20;  
	osuflow->LoadDataCurvilinear((const char*)argv[1], steady, minB, maxB); //true: a steady flow field 
	osuflow->Boundary(minLen, maxLen); // get the boundary 
	minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
	maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
	osuflow->SetBoundary(minB, maxB);  // set the boundary. just to test
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
	return 0;
}


