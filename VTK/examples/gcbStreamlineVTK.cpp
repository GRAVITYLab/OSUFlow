/*

gcbStreamline:
The is a demo to show how to generate streamlines by OSUFlow and how to render the generated streamlines.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>
#include <stdlib.h>
#include <time.h>

#include "gcb.h"
#include "OSUFlow.h"
#include "LineRendererInOpenGL.h"

// VTK
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkSmartPointer.h"
#include "VectorFieldVTK.h"
// VTK


char *szVecFilePath;	// ADD-BY-LEETEN 09/29/2012
OSUFlow *osuflow;
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/07/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/07/2010-END
CLineRendererInOpenGL cLineRenderer;


void draw_cube(float r, float g, float b)
{
  glColor3f(r, g, b);
  glutWireCube(1.0);   // draw a solid cube
}

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
  srand (time(NULL));
	LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 400);
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(.001, .001);
  osuflow->GenStreamLines(sl_list , BACKWARD_AND_FORWARD, 2000, 0);
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	for(int i = 0; i < sl_list.size(); i++)
	{
		VECTOR4 v4Color;
#if 0
		switch((i/2)%7)
		{
		case 0: v4Color = VECTOR4(1.0f, 0.0f, 0.0f, 1.0f);	break;
		case 1: v4Color = VECTOR4(0.0f, 1.0f, 0.0f, 1.0f);	break;
		case 2: v4Color = VECTOR4(0.0f, 0.0f, 1.0f, 1.0f);	break;
		case 3: v4Color = VECTOR4(1.0f, 1.0f, 0.0f, 1.0f);	break;
		case 4: v4Color = VECTOR4(1.0f, 0.0f, 1.0f, 1.0f);	break;
		case 5: v4Color = VECTOR4(0.0f, 1.0f, 1.0f, 1.0f);	break;
		case 6: v4Color = VECTOR4(1.0f, 1.0f, 1.0f, 1.0f);	break;
		}
#else
		switch((i/2)%10)
		{
		case 0: v4Color = VECTOR4(164.f/255.f, 196.f/255.f, 0.0f, 1.0f);	break;
		case 1: v4Color = VECTOR4(96.f/255.f, 169.f/255.f, 23.f/255.f, 1.0f);	break;
		case 2: v4Color = VECTOR4(0, 138.f/255.f, 0, 1.0f);	break;
		case 3: v4Color = VECTOR4(0, 171.f/255.f, 169.f/255.f, 1.0f);	break;
		case 4: v4Color = VECTOR4(27.f/255.f, 161.f/255.f, 226.f/255.f, 1.0f);	break;
		case 5: v4Color = VECTOR4(0, 80.f/255.f, 239.f/255.f, 1.0f);	break;
		case 6: v4Color = VECTOR4(106.f/255.f, 0, 1.f, 1.0f);	break;
		case 7: v4Color = VECTOR4(170.f/255.f, 0, 1.f, 1.0f);	break;
		case 8: v4Color = VECTOR4(244.f/255.f, 114.f/255.f, 208.f/255.f, 1.0f);	break;
		case 9: v4Color = VECTOR4(216.f/255.f, 0, 115.f/255.f, 1.0f);	break;
		}
#endif
		liv4Colors.push_back(v4Color);
	}
	// ADD-BY-LEETEN 07/07/2010-END
	cLineRenderer._Update();
}

void draw_streamlines() 
{
	glPushAttrib(
		GL_LIGHTING_BIT |
		0
	);

	cLineRenderer._Draw();

	glPopAttrib();
}

///////////////////////////////////////////////////////////////////////////////
void
_KeyboardFunc(unsigned char ubKey, int iX, int iY)
{
	switch(ubKey)
	{
	case 's':
		compute_streamlines();
		glutPostRedisplay();
		break;

	case 'h':
		{
			int iHalo;
			cLineRenderer._GetInteger(CLineRenderer::ENABLE_HALO, &iHalo);
			iHalo = !iHalo;
			cLineRenderer._SetInteger(CLineRenderer::ENABLE_HALO, iHalo);
		}
		glutPostRedisplay();
		break;

	case 'l':
		{
			int iLighting;
			cLineRenderer._GetInteger(CLineRenderer::ENABLE_LIGHTING, &iLighting);
			iLighting = !iLighting;
			cLineRenderer._SetInteger(CLineRenderer::ENABLE_LIGHTING, iLighting);
		}

		glutPostRedisplay();
		break;

	// ADD-BY-LEETEN 09/29/2012-BEGIN
	case 'S':
		{
			VECTOR3 v3Min, v3Max;
			osuflow->Boundary(v3Min, v3Max);
			float pfDomainMin[4];
			float pfDomainMax[4];
			for(size_t d = 0; d < 3; d++)
			{
				pfDomainMin[d] = v3Min[d];
				pfDomainMax[d] = v3Max[d];
			}
			pfDomainMin[3] = 0.0f;
			pfDomainMax[3] = 0.0f;

			char szFilename[1024];
			strcpy(szFilename, szVecFilePath);
			strcat(szFilename, ".trace");

			OSUFlow::WriteFlowlines(
				pfDomainMin,
				pfDomainMax,
				&sl_list,
				NULL,
				szFilename);
			LOG(printf("Save the streamlines to %s", szFilename));
		}
		break;
	// ADD-BY-LEETEN 09/29/2012-END

	}
}

void
_DisplayFunc()
{
	//glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render the scene
    draw_streamlines();


    glPushMatrix();
    //glTranslatef(minLen[0], minLen[1], minLen[2]);
    //glScalef(len[0]/len[2], len[1]/len[2], 1.f);
    glScalef(len[0]/7, len[1]/7, len[2]/7);
    //glScalef(len[0], len[1], len[2]);
    draw_cube(1,1,1);
    glPopMatrix();


	// NOTE: Call glutSwapBuffers() at the end of your display function
	glutSwapBuffers();
}

void
init()
{
	LOG(printf("Initialize here."));
	glEnable(GL_DEPTH_TEST);

	// setup light 0
	static GLfloat pfLightAmbient[4] =	{0.1f, 0.1f, 0.1f, 1.0f};
	static GLfloat pfLightColor[4] =	{0.7f, 0.7f, 0.7f, 1.0f};
	glLightfv(GL_LIGHT0, GL_AMBIENT,		pfLightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,		pfLightColor);
	glLightfv(GL_LIGHT0, GL_SPECULAR,		pfLightColor);
	glLightf(GL_LIGHT0, GL_SPOT_EXPONENT,	4.0f);
	cLineRenderer._UpdateLighting();

	// ADD-BY-LEETEN 08/14/2010-BEGIN
	LOG(printf("The vector field is ready. Press key 's' to generate the primtives."));
	// ADD-BY-LEETEN 08/14/2010-END
}

void 
quit()
{
	LOG(printf("Clean up here."));
}

int
main(int argc, char* argv[])
{
	///////////////////////////////////////////////////////////////
	// when use GCB, it is still needed to initialize GLUT
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL );

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
	///////////////////////////////////////////////////////////////
	// initialize OSU flow
	osuflow = new OSUFlow(); 

	// load the scalar field
	LOG(printf("read file %s\n", argv[1])); 

	osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 

	szVecFilePath = argv[1];	// ADD-BY-LEETEN 09/29/2012

	// comptue the bounding box of the streamlines 
	VECTOR3 minB, maxB; 
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
	///////////////////////////////////////////////////////////////
	cLineRenderer._SetBoundingBox(
		minLen[0], minLen[1], minLen[2], 
		maxLen[0], maxLen[1], maxLen[2]);
	cLineRenderer._SetDataSource(&sl_list);
	// ADD-BY-LEETEN /2010-BEGIN
	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_TRACE);
	// ADD-BY-LEETEN 07/07/2010-END

	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Streamline");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbStreamline.cpp,v $
Revision 1.5  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/08/15 04:23:06  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.2  2010/07/07 18:00:33  leeten

[07/07/2010]
1. [ADD] Specify the color for all traces.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
