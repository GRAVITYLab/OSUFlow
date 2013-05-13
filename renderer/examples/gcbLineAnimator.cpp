/*

gcbLineAnimator:
The is a demo to show how to generate streamlines by OSUFlow and how to animate the generated streamlines as moving particles.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "LineAnimatorInOpenGL.h"

enum {ANIME_LINES};

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 
// ADD-BY-LEETEN 07/07/2010-BEGIN
list<VECTOR4> liv4Colors;
// ADD-BY-LEETEN 07/07/2010-END

CLineAnimatorInOpenGL cLineRenderer;

bool bIsAnimationOn;

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
	LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 500); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
	#if	0	// MOD-BY-LEETEN 04/23/2010-FROM:
		  osuflow->SetIntegrationParams(1, 5); 
			osuflow->GenStreamLines(sl_list , BACKWARD_AND_FORWARD, 500, 0); 
	#else	// MOD-BY-LEETEN 04/23/2010-TO:
	osuflow->SetIntegrationParams(0.01, 0.5); 
	osuflow->GenStreamLines(sl_list , FORWARD_DIR, 500, 0); 
	#endif	// MOD-BY-LEETEN 04/23/2010-END
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
	// MOD-BY-LEETEN 04/15/2010-FROM:
		// cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_FRONT_ONE_PARTICLE, 0);
	// TO:
	cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_FRONT_ONE_PARTICLE, 0);
	cLineRenderer._SetInteger(CLineAnimatorInOpenGL::RESET_ALL_PARTICLES, 0);
	// MOD-BY-LEETEN 04/15/2010-END
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
_TimerCB(int iValue)
{
	switch(iValue)
	{
	case ANIME_LINES:
		if( bIsAnimationOn )
		{
			// MOD-BY-LEETEN 04/23/2010-FROM:
				// cLineRenderer._SetInteger(CLineAnimatorInOpenGL::MOVE_ALL_PARTICLES, 1);
			// TO:
			const int iParticleGap = 12;
			const int iParticleLength = 1;
			const int iMaxNrOfParticles = 1000;

			int iFrontParticleIndex; 
			cLineRenderer._GetInteger(CLineAnimatorInOpenGL::FRONT_PARTICLE_INDEX, &iFrontParticleIndex);

			if( iParticleGap == iFrontParticleIndex )
			{
				cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_FRONT_ONE_PARTICLE, iParticleGap);
				for(int i = 1; i < iParticleLength; i++)
					cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_FRONT_ONE_PARTICLE, 1);
			}


			int iBackParticleIndex = 0;
			while(true)
			{
				cLineRenderer._GetInteger(CLineAnimatorInOpenGL::BACK_PARTICLE_INDEX, &iBackParticleIndex);
				if( iBackParticleIndex < iMaxNrOfParticles )
					break;
				cLineRenderer._SetInteger(CLineAnimatorInOpenGL::REMOVE_BACK_PARTICLES, 1);
			} 

			cLineRenderer._SetInteger(CLineAnimatorInOpenGL::MOVE_ALL_PARTICLES, 1);
			// MOD-BY-LEETEN 04/23/2010-END

			cLineRenderer._Update();

			int iNrOfRenderedParticles;
			// MOD-BY-LEETEN 04/15/2010-FROM:
				// cLineRenderer._GetInteger(CLineAnimatorInOpenGL::NR_OF_RENDERED_PARTICLES, &iNrOfRenderedParticles);
			// TO:
			cLineRenderer._GetInteger(CLineRenderer::NR_OF_RENDERED_PARTICLES, &iNrOfRenderedParticles);
			// MOD-BY-LEETEN 04/15/2010-END
			if( 0 == iNrOfRenderedParticles )
			{
				cLineRenderer._SetInteger(CLineAnimatorInOpenGL::RESET_ALL_PARTICLES, 0);
				cLineRenderer._Update();
			}


			// MOD-BY-LEETEN 04/23/2010-TO:
				// glutTimerFunc(66, _TimerCB, ANIME_LINES);
			// TO:
			glutTimerFunc(33, _TimerCB, ANIME_LINES);
			// MOD-BY-LEETEN 04/23/2010-END

			glutPostRedisplay();
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

	// MOD-BY-LEETEN 04/15/2010-FROM:
		// case 'r':
	// TO:
	case 'c':
	// MOD-BY-LEETEN 04/15/2010-END
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::REMOVE_ALL_PARTICLES, 0);
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_FRONT_ONE_PARTICLE, 0);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;

	case '0':
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::RESET_ALL_PARTICLES, 0);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;

	case '+':
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::PUSH_BACK_ONE_PARTICLE, 1);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;

	#if	0	// MOD-BY-LEETEN 04/23/2010-FROM:
		case '-':
			cLineRenderer._SetInteger(CLineAnimatorInOpenGL::DEQUEUE_PARTICLES, 1);
			cLineRenderer._Update();
			glutPostRedisplay();
			break;
	#else	// MOD-BY-LEETEN 04/23/2010-TO:
	case '-':
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::REMOVE_BACK_PARTICLES, 1);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;
	#endif	// MOD-BY-LEETEN 04/23/2010-END

	case '>':
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::MOVE_ALL_PARTICLES, 1);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;

	case '<':
		cLineRenderer._SetInteger(CLineAnimatorInOpenGL::MOVE_ALL_PARTICLES, -1);
		cLineRenderer._Update();
		glutPostRedisplay();
		break;

	case 'a':
		bIsAnimationOn = !bIsAnimationOn;
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
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_TRACE);
	// ADD-BY-LEETEN 07/07/2010-END

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

$Log: gcbLineAnimator.cpp,v $
Revision 1.8  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.7  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.6  2010/08/15 04:24:53  leeten

[08/15/2010]
1. [ADD] Change the description of the application.
2. [ADD] Print a message a the end of _init() to indicate the user to press the hotkey 's'.

Revision 1.5  2010/07/07 17:59:24  leeten

[07/07/2010]
1. [ADD] Specify the color for all traces.

Revision 1.4  2010/04/23 20:31:33  leeten

[04/23/2010]
1. [DEL] Remove extra #endif.

Revision 1.3  2010/04/23 19:24:06  leeten

[04/23/2010]
1. [MOD] Change the names of the parameters.
2. [MOD] When the animation mode is on, particles will be continuously injected

Revision 1.2  2010/04/16 17:42:24  leeten

[04/12/2010]
1. [MOD] When initialize the streamlines, after one particle has been added, reset the index of the first particle to 0.
2. [MOD] Change the parameter's name from CLineAnimatorInOpenGL::NR_OF_RENDERED_PARTICLES to CLineRenderer::NR_OF_RENDERED_PARTICLES.

Revision 1.1  2010/04/15 04:06:35  leeten

[04/12/2010]
1. [1ST] First Time Checkin.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
