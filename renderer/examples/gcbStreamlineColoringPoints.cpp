/*

gcbStreamlineColoringPoints:
The is a demo to show how to generate streamlines by OSUFlow and how to color each point individually. There are 2 steps:

1. Specify the pointers to a list of 4D vector and coloring scheme as COLOR_PER_POINT:

	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_POINT);

2. Specify the color per point to the list: (this part is shown in the function compute_streamlines())

	liv4Colors.clear();
	int iT = 0;
	for(list<vtListSeedTrace*>::const_iterator
			pIter = sl_list.begin(); 
		pIter!=sl_list.end(); 
		pIter++, iT++) 
	{
		const vtListSeedTrace *trace = *pIter; 
		VECTOR4 v4Color(0.0f, 0.0f, 0.0f, 0.0f);
		int iP = 0;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			// specify the color
			...

			liv4Colors.push_back(v4Color);
		}
	}

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "LineRendererInOpenGL.h"

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/07/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/07/2010-END
CLineRendererInOpenGL cLineRenderer;

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
	LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 200); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreamLines(sl_list , FORWARD_DIR, 100, 0); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 

  // Specify the color per point to the list. 
  // The color is computed based on the dot product between 
  // the corresponding streamline segment and a reference vector, 
  // which is (0, 0, 1) in this example.
	liv4Colors.clear();
	int iT = 0;
	for(list<vtListSeedTrace*>::const_iterator
			pIter = sl_list.begin(); 
		pIter!=sl_list.end(); 
		pIter++, iT++) 
	{
		const vtListSeedTrace *trace = *pIter; 
		VECTOR4 v4Color(0.0f, 0.0f, 0.0f, 0.0f);
		int iP = 0;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			VECTOR3 p = **pnIter; 
			static VECTOR3 v3Prev; 
			const VECTOR3 v3RefDir(1.0f, 0.0f, 0.0f); 
			VECTOR3 v3Dir(0.0f, 0.0f, 0.0f);
			if(iP > 0)
			{
				v3Dir = p - v3Prev;
				v3Dir.Normalize();

				float fDot = 0.0f;
				for(int i = 0; i < 3; i++)
					fDot += v3RefDir[i] * v3Dir[i];

				float fAbsDot = fabsf(fDot);
				v4Color = VECTOR4(fAbsDot, 0.0f, 1.0f - fAbsDot, 1.0f);
			}

			v3Prev = p;

			liv4Colors.push_back(v4Color);
		}
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

	}
}

void
_DisplayFunc()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render the scene
    draw_streamlines(); 

	// NOTE: Call glutSwapBuffers() at the end of your display function
	glutSwapBuffers();
}

void
init()
{
	LOG(printf("Initialize here."));
	glEnable(GL_DEPTH_TEST);

	// setup light 0
	static GLfloat pfLightAmbient[4] =	{0.0f, 0.0f, 0.0f, 1.0f};
	static GLfloat pfLightColor[4] =	{0.5f, 0.5f, 0.5f, 1.0f};
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

	///////////////////////////////////////////////////////////////
	// initialize OSU flow
	osuflow = new OSUFlow(); 

	// load the scalar field
	LOG(printf("read file %s\n", argv[1])); 

	osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 

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

	///////////////////////////////////////////////////////////////
	cLineRenderer._SetBoundingBox(
		minLen[0], minLen[1], minLen[2], 
		maxLen[0], maxLen[1], maxLen[2]);
	cLineRenderer._SetDataSource(&sl_list);

	// Specify the pointers to a list of 4D vector and coloring scheme as COLOR_PER_POINT:
	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_POINT);

	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Streamline Coloring Points");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log$

*/
