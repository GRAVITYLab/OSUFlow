/*

gcbFileRenderer:
This is a demo to show how to generated streamlines from the streamline trace files.

Creted by 
Kewei Lu (The Ohio State University)
Oct, 2012

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "TimeLineRendererInOpenGL.h"
#include "LineRendererInOpenGL.h"

char *szVecFilePath;
OSUFlow *osuflow; 
VECTOR4 minLen, maxLen; 
list<vtListTimeSeedTrace*> sl_list; 
float center[3], len[3]; 
list<VECTOR4> liv4Colors;
CTimeLineRendererInOpenGL cLineRenderer;

////////////////////////////////////////////////////////////////////////////
void assignColors() 
{
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
	cLineRenderer._Update();
}

void readTraceFile(char* fileName)
{
	FILE *fpPathline;
	fpPathline = fopen(fileName, "rb");

	assert(fpPathline);
	float pfMin[4];
	float pfMax[4];
	fread(&pfMin, sizeof(pfMin[0]), 4, fpPathline);
	fread(&pfMax, sizeof(pfMax[0]), 4, fpPathline);
	int iNrOfTimeSteps = 1 + (int)ceilf(pfMax[3] - pfMin[3]);
	minLen[0] = pfMin[0];minLen[1] = pfMin[1];minLen[2] = pfMin[2];minLen[3] = pfMin[3];
	maxLen[0] = pfMax[0];maxLen[1] = pfMax[1];maxLen[2] = pfMax[2];maxLen[3] = pfMax[3];

	size_t uFOffset = ftell(fpPathline);
	int iNrOfPathlines = 0;
	for(int iDelim = 0; iDelim >= 0; iNrOfPathlines++)
		fread(&iDelim, sizeof(iDelim), 1, fpPathline);
	iNrOfPathlines--;

	size_t uDataStart = ftell(fpPathline);

	fseek(fpPathline, uFOffset, SEEK_SET);

	int* piNrsOfCoords;
	piNrsOfCoords = (int*)malloc(iNrOfPathlines*sizeof(int));
	fread(piNrsOfCoords, sizeof(piNrsOfCoords[0]), iNrOfPathlines, fpPathline);

	// skip the delim
	fseek(fpPathline, uDataStart, SEEK_SET);

	int iTotalNrOfPoints = 0;
	for(int s = 0; s < iNrOfPathlines; s++)
	{
		int iNrOfPoints = piNrsOfCoords[s];
		if( 0 == iNrOfPoints )
		{
			LOG("Error");
			continue;
		}

		vtListTimeSeedTrace* vtNewListSeedTrace = new vtListTimeSeedTrace;
		vtNewListSeedTrace->clear();
		VECTOR4* pv4Coords;
		pv4Coords = (VECTOR4*) malloc(iNrOfPoints*sizeof(VECTOR4));
		fread(&pv4Coords[0], sizeof(pv4Coords[0]), iNrOfPoints, fpPathline);

		float T;
		float prevT;
		for(int c = 0; c < iNrOfPoints; c++)
		{
			T = pv4Coords[c][3];
			if(c!=0 && prevT == T)
				LOG(printf("WARNING:Two consequence points in the line list have the same value for time step, you should not use TimeLineRendererInOpenGL"));
			vtNewListSeedTrace->push_back(new VECTOR4(pv4Coords[c]));
			prevT = T;
		}

		sl_list.push_back(vtNewListSeedTrace);
		iTotalNrOfPoints += iNrOfPoints;
	}
	fclose(fpPathline);
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

	// load the scalar field
	LOG(printf("read file %s\n", argv[1]));
	readTraceFile(argv[1]);
	assignColors();

	szVecFilePath = argv[1];

	// comptue the bounding box of the streamlines 
	VECTOR3 minB, maxB; 
	minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
	maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
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
	cLineRenderer._SetColorSource(&liv4Colors);
	cLineRenderer._SetInteger(CLineRenderer::COLOR_SCHEME, CLineRenderer::CColorScheme::COLOR_PER_TRACE);

	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Line Rendering from Files");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MIN_TIME_STEP,0.0f);
	cLineRenderer._SetFloat(CTimeLineRendererInOpenGL::MAX_TIME_STEP,1.0f);
	cLineRenderer._Update();
	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

[02/10/2012]
1. [1ST] First Time Checkin.

*/
