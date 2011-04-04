/*

gcbPathLine
The is a demo to show how to generate pathlines by OSUFlow and how to animate the generated pathlines.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "TimeLineRendererInOpenGL.h"

enum {ANIME_LINES};

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListTimeSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/09/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/09/2010-END

CTimeLineRendererInOpenGL cLineRenderer;

bool bIsAnimationOn;

int num_timesteps; 
int num_frames = 50; 
int current_frame = 0; 
float time_incr; 

////////////////////////////////////////////////////////////////////////////
void compute_pathlines() 
{

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 1000); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  float* tarray = new float[nSeeds]; 
  for (int i=0;i<nSeeds; i++) 
    tarray[i] = (float)(i % num_timesteps); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(0.01f, 0.5f); 
  //  osuflow->GenPathLines(seeds, sl_list , FORWARD, nSeeds, 5000); 
  osuflow->GenPathLines(seeds, sl_list , FORWARD, nSeeds, 5000, tarray); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 

	// ADD-BY-LEETEN 07/09/2010-BEGIN
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
	// ADD-BY-LEETEN 07/09/2010-END
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

void 
_MoveFrame()
{
	float min_time = float(current_frame);
	float max_time = float(current_frame+1); 
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MIN_TIME_STEP, min_time);
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MAX_TIME_STEP, max_time);
	cLineRenderer._Update();
	glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////////////
void 
_TimerCB(int iValue)
{
	switch(iValue)
	{
	case ANIME_LINES:
		if( bIsAnimationOn )
		{
			_MoveFrame();
			glutTimerFunc(66, _TimerCB, ANIME_LINES);
			current_frame = (current_frame+1) % num_timesteps; 
		}
		break;
	}
}

void
_KeyboardFunc(unsigned char ubKey, int iX, int iY)
{
	switch(ubKey)
	{
	case 's':
		compute_pathlines();
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

	case '0':
		current_frame = 0;
		_MoveFrame();
		break;

	case '>':
		current_frame = (current_frame + 1) % num_timesteps;
		_MoveFrame();
		break;

	case '<':
		current_frame--;
		if( current_frame < 0 )
			current_frame = num_timesteps - 1;
		_MoveFrame();
		break;

	case 'a':
		bIsAnimationOn = !bIsAnimationOn;
		if( false == bIsAnimationOn )
		{
			cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MIN_TIME_STEP, 0.0f);
			cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MAX_TIME_STEP, float(num_timesteps));
			cLineRenderer._Update();
			glutPostRedisplay();
		}
		else
			_TimerCB(ANIME_LINES);
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
_IdleFunc()
{
	glutPostRedisplay();
}

void
init()
{
	cLineRenderer._SetInteger(CTimeLineRendererInOpenGL::TIME_VARYING_ON, 1);
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MIN_TIME_STEP, 0.0f);
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MAX_TIME_STEP, float(num_timesteps));


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
	osuflow->LoadData((const char*)argv[1], false); //false : a time-varying flow field 

	osuflow->ScaleField(100.0); 

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

	num_timesteps = osuflow->NumTimeSteps(); 
	time_incr = num_timesteps/(float) num_frames; 
	printf(" reading in %d time steps.\n", num_timesteps); 

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
	// ADD-BY-LEETEN 07/09/2010-END

	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Line Animator");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbPathline.cpp,v $
Revision 1.2  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/08/15 04:23:06  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.2  2010/07/09 20:07:34  leeten

[07/09/2010]
1. [ADD] Support the coloring of streamlines.

Revision 1.1  2010/04/16 17:37:51  leeten

[04/12/2010]
1. [1ST] First Time Checkin.

Revision 1.1  2010/04/15 04:06:35  leeten

[04/12/2010]
1. [1ST] First Time Checkin.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
