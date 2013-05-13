/*

gcbStreamlineDynamicColor:
The is a demo to show how to generate streamlines by OSUFlow, how to render the generated streamlines, 
and how to change the streamline color on the fly.

To generate the streamline colors on the fly, the following steps should be followed:

1. Create a new class that inherited from CLineRendererInOpenGL and overload the following tow methods:

	void _CheckTrace(int iTrace, bool& bIsDrawingTrace);
	void _GetTraceColor(int iTrace, float& fR, float& fG, float& fB, float& fA);

Dependent on the applications, this new class can 
define other methods to control the coloring. In this sample code, two methods are defined to randomly 
decide whethere a streamline is drawn or not (_GenearateRandomFlags()) and to randomly assign streamline 
colors (_GenearateRandomColor()). Also, I overload the method _TraverseLinesBegin() in order to allocate 
the needed data structures for all streamlines. 

Then declare the a renderer of the new created class. 

2. Before calling cLineRenderer._Update(), set the coloring scheme to be "On the Fly" by executing the 
following statements:

	_SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_ON_THE_FLY);

3. Regenerate the colors according to the application's requirement. In this example, when the user presses 
the hot keys '0' and '1' to regenerate the colors and flags, respectively, of all streamlines.

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
list<VECTOR4> liv4Colors;
class CDynamicColorLineRenderer:
	public CLineRendererInOpenGL 
{
	vector<bool> vbTraceFlags;
	vector<VECTOR4> vv4TraceColors;
	virtual void _CheckTrace
	(
		int iTrace,
		bool& bIsDrawingTrace
	)
	{
		bIsDrawingTrace = vbTraceFlags[iTrace];
	}

	virtual void _GetTraceColor
	(
		int iTrace, 
		float& fR, 
		float& fG, 
		float& fB, 
		float& fA
	)
	{
		VECTOR4 v4Color = vv4TraceColors[iTrace];
		fR = v4Color[0];
		fG = v4Color[1];
		fB = v4Color[2];
		fA = v4Color[3];
	}

public:
	#if	0	// MOD-By-LEETEN 01/20/2011-FROM:
		void _SetNrOfTraces(int iNrOfTraces)
		{
			vbTraceFlags.resize(iNrOfTraces);
			vv4TraceColors.resize(iNrOfTraces);
		}
	#else	// MOD-By-LEETEN 01/20/2011-TO:
	void _TraverseLinesBegin(int iNrOfTraces)
	{
		vbTraceFlags.resize(iNrOfTraces);
		vv4TraceColors.resize(iNrOfTraces);
		CLineRendererInOpenGL::_TraverseLinesBegin(iNrOfTraces);
	}
	#endif	// MOD-By-LEETEN 01/20/2011-END
	
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	void _GenearateRandomFlags()
	{
		for(int t = 0; t < vbTraceFlags.size(); t++)
		{
			vbTraceFlags[t] = (rand()%2)?true:false;
		}
	}
	// ADD-BY-LEETEN 01/20/2011-END

	void _GenearateRandomColor()
	{
		for(int t = 0; t < vbTraceFlags.size(); t++)
		{
			// DEL-BY-LEETEN 01/20/2011-BEGIN
				// vbTraceFlags[t] = (rand()%2)?true:false;
			// DEL-BY-LEETEN 01/20/2011-END
			VECTOR4 v4Color;
			int iRand = rand() % 7;
			switch(iRand)
			{
			case 0: v4Color = VECTOR4(1.0f, 0.0f, 0.0f, 1.0f);	break;
			case 1: v4Color = VECTOR4(0.0f, 1.0f, 0.0f, 1.0f);	break;
			case 2: v4Color = VECTOR4(0.0f, 0.0f, 1.0f, 1.0f);	break;
			case 3: v4Color = VECTOR4(1.0f, 1.0f, 0.0f, 1.0f);	break;
			case 4: v4Color = VECTOR4(1.0f, 0.0f, 1.0f, 1.0f);	break;
			case 5: v4Color = VECTOR4(0.0f, 1.0f, 1.0f, 1.0f);	break;
			case 6: v4Color = VECTOR4(1.0f, 1.0f, 1.0f, 1.0f);	break;
			}
			vv4TraceColors[t] = v4Color;
		}
	}
};
CDynamicColorLineRenderer cLineRenderer;

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
	LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 100); 
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
  printf("list size = %d\n", (int)sl_list.size()); 

	#if	0	// DEL-BY-LEETEN 01/20/2011-BEGIN
	  cLineRenderer._SetNrOfTraces(sl_list.size());
	  cLineRenderer._GenearateRandomColor();
	#endif	// DEL-BY-LEETEN 01/20/2011-END

	cLineRenderer._Update();
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	cLineRenderer._GenearateRandomColor();
	// ADD-BY-LEETEN 01/20/2011-END
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	cLineRenderer._GenearateRandomFlags();
	// ADD-BY-LEETEN 01/20/2011-END
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
	case '0': 
		cLineRenderer._GenearateRandomColor();
		glutPostRedisplay();
		break;

	// ADD-BY-LEETEN 01/20/2011-BEGIN
	case '1': 
		cLineRenderer._GenearateRandomFlags();
		glutPostRedisplay();
		break;
	// ADD-BY-LEETEN 01/20/2011-END

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
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_ON_THE_FLY);

	///////////////////////////////////////////////////////////////
	// MOD-BY-LEETEN 01/20/2011-FROM:
		// glutCreateWindow("GCB Streamline");
	// TO:
	glutCreateWindow("GCB Streamline Dynamic Color");
	// MOD-BY-LEETEN 01/20/2011-END

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbStreamlineDynamicColor.cpp,v $
Revision 1.2  2011/01/21 20:51:41  leeten

[01/21/2011]
1. [ADD] Add explaination of the sample codes in the comment.
2. [ADD] Add a method _GenearateRandomFlags() to randonly generate the flags. This method will be invoked when the hotkey '1' is pressed.

Revision 1.1  2011/01/20 17:20:33  leeten

[01/19/2010]
1. [ADD] First time checkin.

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
