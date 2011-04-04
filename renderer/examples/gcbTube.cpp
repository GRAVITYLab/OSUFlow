/*

gcbTube:
The is a demo to show how to generate streamlines as tubes by OSUFlow and how to render the generated tubes.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "TubeRendererInOpenGL.h"

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/07/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/07/2010-END

CTubeRendererInOpenGL cLineRenderer;

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
  osuflow->GenStreamLines(sl_list , BACKWARD_AND_FORWARD, 30, 0); 
  printf(" done integrations\n"); 
    printf("list size = %d\n", (int)sl_list.size()); 

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	for(int i = 0; i < sl_list.size(); i++)
	{
		VECTOR4 v4Color;
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

	int iLighting;
	cLineRenderer._GetInteger(CLineRenderer::ENABLE_LIGHTING, &iLighting);
	if( iLighting )
	{
		static GLfloat pfLightAmbient[4] =	{0.1f, 0.1f, 0.1f, 1.0f};
		static GLfloat pfLightDiffuse[4] =	{0.6f, 0.6f, 0.6f, 1.0f};
		static GLfloat pfLightSpecular[4] =	{0.3f, 0.3f, 0.3f, 1.0f};;
		static GLfloat fSpotExponent = 4.0f;
		glLightfv(GL_LIGHT0, GL_AMBIENT,	pfLightAmbient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE,	pfLightDiffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR,	pfLightSpecular);
		glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, fSpotExponent);

		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);

		glPushMatrix();
		glLoadIdentity();
		static GLfloat pfLightPos[4] =	{0.0f, 0.0f, 1.0f, 0.0f};
		glLightfv(GL_LIGHT0, GL_POSITION, pfLightPos);
		glPopMatrix();
	}

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

	case 'w':
		// MOD-BY-LEETEN 01/20/2011-FROM:
			// cLineRenderer._SetInteger(CTubeRenderer::DRAW_POLYGON, CTubeRenderer::DRAW_POLYGON_WIREDFRAME);
		// TO:
		{
			static int iWiredFrameOrFill = 0;
			iWiredFrameOrFill = 1 - iWiredFrameOrFill;
			switch(iWiredFrameOrFill)
			{
			case 0: cLineRenderer._SetInteger(CTubeRenderer::DRAW_POLYGON, CTubeRenderer::DRAW_POLYGON_WIREDFRAME);	break;
			case 1:	cLineRenderer._SetInteger(CTubeRenderer::DRAW_POLYGON, CTubeRenderer::DRAW_POLYGON_FILL);	break;
			}
		}
		// MOD-BY-LEETEN 01/20/2011-END
		glutPostRedisplay();
		break;

	#if	0	// DEL-BY-LEETEN 01/20/2011-BEGIN
	case 'f':
		cLineRenderer._SetInteger(CTubeRenderer::DRAW_POLYGON, CTubeRenderer::DRAW_POLYGON_FILL);
		glutPostRedisplay();
		break;
	#endif	// DEL-BY-LEETEN 01/20/2011-END

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
main(int argn, char* argv[])
{
	///////////////////////////////////////////////////////////////
	// when use GCB, it is still needed to initialize GLUT
	glutInit(&argn, argv);
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
	// ADD-BY-LEETEN /2010-BEGIN
	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_TRACE);
	// ADD-BY-LEETEN 07/07/2010-END

	///////////////////////////////////////////////////////////////
	// MOD-BY-LEETEN 10/01/2010-FROM:
		// glutCreateWindow("GCB Streamline");
	// TO:
	glutCreateWindow("GCB Tube");
	// MOD-BY-LEETEN 10/01/2010-END


	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbTube.cpp,v $
Revision 1.8  2011/01/21 20:50:59  leeten

[01/21/2011]
1. [MOD] Remove the hotkey 'f'. Now the rendering mode will be toggeled only by the hotkey 'w'.

Revision 1.7  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.5  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.4  2010/08/15 04:24:53  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.3  2010/08/15 04:23:06  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.2  2010/07/07 18:01:03  leeten

[07/07/2010]
1. [ADD] Specify the color for all traces.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
