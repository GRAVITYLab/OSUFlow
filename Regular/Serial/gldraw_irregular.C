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

char graphbin[255], tetragrid[255],tetrasoln[255],tetratetra[255];

//draw grid
int dim[3];
float* ver;
int  * hexs;
int ver_num,hex_num;
GLuint grid_list;

float * grid_fx, *grid_fy, *grid_fz;

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

	char filename[255];

	int xdim,ydim,zdim;
	float refval[4];
	// get the data size
	//sprintf(slFn, "%s%d", m_strSlnFN, 0); 
	FILE * fIn = fopen(tetrasoln, "rb");

	fread(&xdim, sizeof(int), 1, fIn);
	fread(&ydim, sizeof(int), 1, fIn);
	fread(&zdim, sizeof(int), 1, fIn);
	fread(refval, sizeof(float), 4, fIn);
	int numToRead = xdim*ydim*zdim;
	int nodeNum = numToRead;

	float* pDensity = new float[numToRead];
	float* pMomentumU = new float[numToRead];
	float* pMomentumV = new float[numToRead];
	float* pMomentumW = new float[numToRead];
	float* pEnergy = new float[numToRead];

	// read data
	numToRead = xdim * ydim * zdim;
	
	// density, momentum
	fread(pDensity, sizeof(float), numToRead, fIn);
	fread(pMomentumU, sizeof(float), numToRead, fIn);
	fread(pMomentumV, sizeof(float), numToRead, fIn);
	fread(pMomentumW, sizeof(float), numToRead, fIn);
	fread(pEnergy, sizeof(float), numToRead, fIn);

	fclose(fIn);
//	printf("Load in flow data %s!\n", filename);

	//get seeds
	FILE* fp=fopen(tetragrid,"rb");	
	int nodenum,surf,tetra;
	fread(&nodenum,sizeof(int),1,fp);
	fread(&surf,sizeof(int),1,fp);
	fread(&tetra,sizeof(int),1,fp);
	float*x,*y,*z;
	x=new float[nodenum];
	y=new float[nodenum];
	z=new float[nodenum];
	fread(x,sizeof(float),nodenum,fp);
	fread(y,sizeof(float),nodenum,fp);
	fread(z,sizeof(float),nodenum,fp);
	fclose(fp);
	fp=fopen(tetratetra,"rb");
//	int nodenum,surf,tetra;
	fread(&nodenum,sizeof(int),1,fp);
	fread(&surf,sizeof(int),1,fp);
	fread(&tetra,sizeof(int),1,fp);
	int* tetraIds;
	tetraIds=new int[tetra*4];
	fread(tetraIds,sizeof(int),tetra*4,fp);
	fclose(fp);
	float from[3], to[3]; 

	printf("generating seeds...\n"); 

	//generate random seeds
	std::vector<VECTOR3> seedslist;
	int max_seed_num=10;
	int count=0;
	srand((unsigned)time(NULL));			// initialize random number generator
	while(count<max_seed_num)
	{
//	for(int i = 0; i < tetra; i++)
//	{
		int i = (float)rand()/(float)RAND_MAX*(nodenum-1);
		for(int j=0;j<4;j++)
		{	
			int id;
			id=tetraIds[i*4+j];

			VECTOR3 tmp;
			tmp.Set(pMomentumU[id],pMomentumV[id],pMomentumW[id]);
			float mag=tmp.GetMag();
		//	printf("mag=%f\n",mag);
			if((mag<0.17) && (mag>0.13))
			{
				//pick the index of node  0 for the tetra
				seedslist.push_back(VECTOR3(x[id],y[id],z[id]));	
				printf("seed=%f %f %f\n",x[id],y[id],z[id]);
				count++;
				break;
			}
		}

	
	}
	delete [] pDensity;
	delete [] pMomentumU;
	delete [] pMomentumV;
	delete [] pMomentumW;
	delete [] pEnergy;

	int seed_num=seedslist.size();
	VECTOR3* seed_list=new VECTOR3[seed_num];
	for(int i=0;i<seed_num;i++)
	{
		seed_list[i].Set(seedslist[i].x(),seedslist[i].y(),seedslist[i].z());
	}
	
	printf("finished generating seeds...\n"); 

	//seed_list[0].Set(0.770097,-0.562795,0.334426-0.2);
	//seed_list[1].Set(0.344697,-1.07905,0.382389-0.2);

	delete [] x;
	delete [] y;
	delete [] z;
	delete [] tetraIds;
//	seed_num=count;
	osuflow->SetSeedPoints(seed_list, seed_num);
	// osuflow->SetRandomSeedPoints(from, to, 10); 
	int nSeeds; 
	VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
	for (int i=0; i<nSeeds; i++) 
		printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
		seeds[i][1], seeds[i][2]); 

	sl_list.clear(); 
	printf("compute streamlines..\n"); 
	osuflow->SetIntegrationParams(1, 5); 
	osuflow->GenStreamLines(sl_list , BACKWARD_AND_FORWARD, 500, 0); 
	printf(" done integrations\n"); 
	printf("list size = %d\n", sl_list.size()); 

}

void draw_streamlines() {

	glPushMatrix(); 

	glLineWidth(2);
//	glScalef(1/(float)len[0], 1/(float)len[0], 1/(float)len[0]); 
//	glTranslatef(-len[0]/2.0, -len[1]/2.0, -len[2]/2.0); 

//	printf("draw streamlines.\n"); 
	glColor3f(1,0,0); 
	std::list<vtListSeedTrace*>::iterator pIter; 
	glLineWidth(2);

	pIter = sl_list.begin(); 
	for (; pIter!=sl_list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter; 
		std::list<VECTOR3*>::iterator pnIter; 
		pnIter = trace->begin(); 
	//	if(trace->size()<20)
	//		continue;
	//	printf(" len=%d\n ",trace->size()); 
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
//	draw_grid();
	glDisable(GL_LIGHTING);

	draw_streamlines(); 
	if(!glIsList(grid_list))
	{
	grid_list = glGenLists(1);
	glNewList(grid_list, GL_COMPILE);

	glEndList();
	}
	glCallList(grid_list);
	glPopMatrix(); 
	


	//printf(" len %f %f %f\n", len[0], len[1], len[2]); 
	glutSwapBuffers(); 
}


void mykey(unsigned char key, int x, int y)
{
	switch(key) {
	case 'q': exit(1);
		break; 
	case 's': compute_streamlines(); 
		glutPostRedisplay(); 
		break; 
	
	}
}
///////////////////////////////////////////////////////////////
void readGrid(char* filename)
{
	char gird_file[255];
	memset(gird_file,0,255);
	sprintf(gird_file,"%s.grid",filename);
	FILE* fIn = fopen(gird_file, "rb");


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
void read_ANL_grid(char* filename)
{
	char gird_file[255];
	memset(gird_file,0,255);
	sprintf(gird_file,"%s",filename);

	FILE* fIn = fopen(gird_file, "rb");

	fread(&ver_num, sizeof(int), 1, fIn);
	ver=new float[3*ver_num];
	fread(ver, sizeof(float), 3*ver_num, fIn);
	
	fread(&hex_num, sizeof(int), 1, fIn);
	hexs=new int[8*hex_num];
	fread(hexs, sizeof(int), 8*hex_num, fIn);

	fread(&ver_num, sizeof(int), 1, fIn);
	float* vel=new float[3*ver_num];
	fread(vel, sizeof(float), 3*ver_num, fIn);

	fclose(fIn);	
}
int main(int argc, char** argv) 
{
//argv: 1, grid, 2. tetra 3. solution, 4. tetra grid
	
	memset(graphbin,0,255);
	sprintf(graphbin,"%s",argv[1]);
	memset(tetragrid,0,255);
	sprintf(tetragrid,"%s.grid",argv[2]);
	memset(tetrasoln,0,255);
	sprintf(tetrasoln,"%s.soln",argv[2]);
	memset(tetratetra,0,255);
	sprintf(tetratetra,"%s.tetra",argv[2]);

	read_ANL_grid(graphbin);

	VECTOR3 minB, maxB; 

	osuflow = new OSUFlow(); 
	printf("read file %s\n", argv[1]); 
	minB[0] = 0; minB[1] = 0; minB[2] = 0; 
	maxB[0] = 20; maxB[1] = 20; maxB[2] = 20;  

	osuflow->LoadDataIrregular((const char*)argv[2], true, minB, maxB) ;
	int atphy;
	VECTOR3 pos,vecdata;
	PointInfo info;
	float from[3], to[3]; 

	printf("generating seeds...\n"); 
//  osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 
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
//	glutTimerFunc(10, timer, 0); 
//	glutMouseFunc(mymouse); 
//	glutMotionFunc(mymotion);
	glutKeyboardFunc(mykey); 
	glutMainLoop(); 
	return 0;
}


