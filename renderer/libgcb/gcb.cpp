/*

GCB (Glut CallBacks): this is a template for developing GLUT-vbased OpenGL programs. 
My goal is to use this file as a template to speed up the development of OpenGL/GLUT 
programs

Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May,2010

*/

#include <string.h>

#include "gcb.h"

///////////////////////////////////////////////////////////////
// 

// ADD-BY-LEETEN 10/01/2010-BEGIN
double dRotateCameraSpeed = 0.1;
double dRotateModelSpeed = 0.01;
double dZoomModelSpeed = 0.001;
// ADD-BY-LEETEN 10/01/2010-END

typedef GLdouble TMatrix[16];
TMatrix tViewMatrix, tInitViewMatrix;
TMatrix tModelMatrix, tInitModelMatrix;
TMatrix tProjectionMatrix;

int piViewport[4];

// ADD-BY-LEETEN 03/31/2008-BEGIN
// the different between the origin of the model and the cursor the on screen.
// it is used to adjust the shift of the model
double dXDiffOnScreen, dYDiffOnScreen;	

///////////////////////////////////////////////////////////////
// GLUT Callbacks 
void (*_DisplayFunc)();
void (*_ReshapeFunc)(int w, int h);
void (*_KeyboardFunc)(unsigned char key, int w, int h);
void (*_SpecialFunc)(int skey, int w, int h);
void (*_IdleFunc)();

void (*QuitFunc)();

///////////////////////////////////////////////////////////////
// The callbacks, functions and vairable for the mouse interfaces.
// The reference point of every kind of motion is the point where the mouse is just clicked.
static int iBeginX, iBeginY; // the 2D coordinate of the cursor
static int iCursorX, iCursorY; // the 2D coordinate of the cursor
static GLenum eMouseButton = 0;
static GLenum eModifier = 0;
static bool bMoving = false;

void
_Quit()
{
	if( QuitFunc )
		QuitFunc();
}

void
_AlignPan(double pdNewCoord[3], double pdOldCoord[3])
{
	if( 0 == (GLUT_ACTIVE_ALT & eModifier) )
		return;

	size_t uMaxDir = 0;
	for(size_t d = 1; d<3; d++)
		if( fabs(pdNewCoord[d] - pdOldCoord[d]) > 
			fabs(pdNewCoord[uMaxDir] - pdOldCoord[uMaxDir]) )
			uMaxDir = d;

	for(size_t d = 0; d<3; d++)
		if( d != uMaxDir )
			pdNewCoord[d] = pdOldCoord[d];
}

void
_AlignRotate(double pdAxis[3])
{
	if( 0 == (GLUT_ACTIVE_ALT & eModifier) )
		return;

	size_t uMaxDir = 0;	// X
	for(size_t d = 1; d < 3; d++)
		if( fabs(pdAxis[d]) > fabs(pdAxis[uMaxDir]) )
			uMaxDir = d;

	for(size_t d = 0; d < 3; d++)
		pdAxis[d] = (uMaxDir == d)?((pdAxis[d]>0.0)?1.0:-1.0):0.0;
}

// ADD-BY-LEETEN 03/31/2008-END

///////////////////////////////////////////////////////////////
// functions to manipulate the model and camera
void
trackball_ptov(int x, int y, int width, int height, double v[3])
{
    double d, a;
    /* project x,y onto a hemisphere centered within width, height , note z is up here*/
    v[0] = (double)(2*x - width) / (double)width;
    v[1] = (double)(2*y - height) / (double)height;    
    d = sqrt(v[0]*v[0] + v[1]*v[1]);
	v[2] = 0.0;
	if( d < 1.0 )
		v[2] = cos(M_PI_2 * d);
    a = 1.0f / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *=a;    v[1] *= a;    v[2] *= a;
}

static double pdOldCoord[3], pdNewCoord[3];
static double pdCurPos[3], pdBeginPos[3];
static double dX, dY, dZ;

void
_RotateCamera()
{

	// trackball_ptov(iBeginX, iBeginY, piViewport[2], piViewport[3], pdBeginPos);
	// trackball_ptov(iCursorX, iCursorY, piViewport[2], piViewport[3], pdCurPos);
	trackball_ptov(piViewport[2]/2, piViewport[3]/2, piViewport[2], piViewport[3], pdBeginPos); 
	trackball_ptov(piViewport[2]/2 + iCursorX - iBeginX, piViewport[3]/2 + iCursorY - iBeginY, piViewport[2], piViewport[3], pdCurPos);

	dX = pdCurPos[0] - pdBeginPos[0];
	dY = pdCurPos[1] - pdBeginPos[1];
	dZ = pdCurPos[2] - pdBeginPos[2];

	if (dX || dY || dZ) 
	{
		/* compute theta and cross product */
		double pdAxis[3];
		// MOD-BY-LEETEN 10/01/2010-FROM:
			// double dAngle = 0.1 * 90.0 * sqrt(dX*dX + dY*dY + dZ*dZ);
		// TO:
		double dAngle = dRotateCameraSpeed * 90.0 * sqrt(dX*dX + dY*dY + dZ*dZ);
		// MOD-BY-LEETEN 10/01/2010-END

		pdAxis[0] = pdBeginPos[1] * pdCurPos[2] - pdBeginPos[2] * pdCurPos[1];
		pdAxis[1] = pdBeginPos[2] * pdCurPos[0] - pdBeginPos[0] * pdCurPos[2];
		pdAxis[2] = pdBeginPos[0] * pdCurPos[1] - pdBeginPos[1] * pdCurPos[0];

		// ADD-BY-LEETEN 03/31/2008-BEGIN
		_AlignRotate(pdAxis);
		// ADD-BY-LEETEN 03/31/2008-END

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotated(dAngle, pdAxis[0], pdAxis[1], pdAxis[2]);
		glMultMatrixd(tViewMatrix);
		glGetDoublev(GL_MODELVIEW_MATRIX, tViewMatrix);

		glutPostRedisplay();
	}
}

void _PanCamera()
{
	gluUnProject(
		(double)iBeginX, (double)iBeginY, 0.0, 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdOldCoord[0], &pdOldCoord[1], &pdOldCoord[2]);

	gluUnProject(
		(double)iCursorX, (double)iCursorY, 0.0, 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdNewCoord[0], &pdNewCoord[1], &pdNewCoord[2]);

	// ADD-BY-LEETEN 03/31/2008-BEGIN
	_AlignPan(pdNewCoord, pdOldCoord);
	// ADD-BY-LEETEN 03/31/2008-END

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(tViewMatrix);
	glTranslated(pdOldCoord[0] - pdNewCoord[0], (pdOldCoord[1] - pdNewCoord[1]), pdOldCoord[2] - pdNewCoord[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, tViewMatrix);
	glutPostRedisplay();
}

void 
_ZoomCamera()
{
	gluUnProject(
		(double)piViewport[2]/2.0, (double)piViewport[3]/2.0, ((double)(iBeginY)/(double)piViewport[3]) * 2.0 - 1.0, 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdOldCoord[0], &pdOldCoord[1], &pdOldCoord[2]);

	gluUnProject(
		(double)piViewport[2]/2.0, (double)piViewport[3]/2.0, ((double)(iCursorY)/(double)piViewport[3]) * 2.0 - 1.0, 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdNewCoord[0], &pdNewCoord[1], &pdNewCoord[2]);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(tViewMatrix);
	glTranslated(pdOldCoord[0] - pdNewCoord[0], (pdOldCoord[1] - pdNewCoord[1]), pdOldCoord[2] - pdNewCoord[2]);
	glGetDoublev(GL_MODELVIEW_MATRIX, tViewMatrix);
	glutPostRedisplay();
}

void
_PanModel()
{
	gluProject(
		tModelMatrix[12], tModelMatrix[13], tModelMatrix[14], 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdOldCoord[0], &pdOldCoord[1], &pdOldCoord[2]);

	gluUnProject(
		// MOD-BY-LEETEN 03/31/2008-BEGIN
		// FROM: (double)iCursorX, (double)iCursorY, pdOldCoord[2],
		// TO:
		(double)iCursorX + dXDiffOnScreen, (double)iCursorY + dYDiffOnScreen, pdOldCoord[2],
		// MOD-BY-LEETEN 03/31/2008-END
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdNewCoord[0], &pdNewCoord[1], &pdNewCoord[2]);

	for(int i=0; i<3; i++)
		tModelMatrix[12 + i] = pdNewCoord[i];

	glutPostRedisplay();
}

void
_RotateModel()
{
	trackball_ptov(piViewport[2]/2, piViewport[3]/2, piViewport[2], piViewport[3], pdBeginPos); 
	trackball_ptov(piViewport[2]/2 + iCursorX - iBeginX, piViewport[3]/2 + iCursorY - iBeginY, piViewport[2], piViewport[3], pdCurPos);

	dX = pdCurPos[0] - pdBeginPos[0];
	dY = pdCurPos[1] - pdBeginPos[1];
	dZ = pdCurPos[2] - pdBeginPos[2];

	if (dX || dY || dZ) 
	{
		/* compute theta and cross product */
		double pdAxis[3];
		// MOD-BY-LEETEN 10/01/2010-FROM:
			// double dAngle = 0.1 * 90.0 * sqrt(dX*dX + dY*dY + dZ*dZ);
		// TO: 
		double dAngle = dRotateModelSpeed * 90.0 * sqrt(dX*dX + dY*dY + dZ*dZ);
		// MOD-BY-LEETEN 10/01/2010-END
		pdAxis[0] = pdBeginPos[1] * pdCurPos[2] - pdBeginPos[2] * pdCurPos[1];
		pdAxis[1] = pdBeginPos[2] * pdCurPos[0] - pdBeginPos[0] * pdCurPos[2];
		pdAxis[2] = pdBeginPos[0] * pdCurPos[1] - pdBeginPos[1] * pdCurPos[0];

		// ADD-BY-LEETEN 03/31/2008-BEGIN
		_AlignRotate(pdAxis);
		// ADD-BY-LEETEN 03/31/2008-END

		// Produce a new model matrix M' s.t.
		// the rotation R should be consistent to the final camera coordinate
		// i.e., M * V * R = M' * V 
		// => M' = M * V * R * V^-1 = M * V * R * V^T
		// Here assums that the transform of the view matrix is always a rigid transform.
		// Note the since only rotation is considered, the translation offsets are ignored 
		// during the multipilication. 
		static TMatrix tTempMatrix;

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		// reset the temp. matrix
		memset(tTempMatrix, 0, sizeof(tTempMatrix));

		// apply transpose
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				tTempMatrix[i * 4 + j] = tViewMatrix[j * 4 + i];

		for(int i=0; i<3; i++)
			tTempMatrix[4 * i + 3] = 0;

		glMultMatrixd(tTempMatrix);

		// now apply the rotation
		glRotated(dAngle, pdAxis[0], pdAxis[1], pdAxis[2]);

		// apply the view matrix
		memcpy(tTempMatrix, tViewMatrix, sizeof(tTempMatrix));
		memset(&tTempMatrix[12], 0, 3 * sizeof(tTempMatrix[0]));
		glMultMatrixd(tTempMatrix);

		// apply the original model matrix
		memcpy(tTempMatrix, tModelMatrix, sizeof(tTempMatrix));
		memset(&tTempMatrix[12], 0, 3 * sizeof(tTempMatrix[0]));
		glMultMatrixd(tTempMatrix);

		// get the final matrix
		glGetDoublev(GL_MODELVIEW_MATRIX, tTempMatrix);

		// restore the translation offset
		memcpy(tModelMatrix, tTempMatrix, 12 * sizeof(tTempMatrix[0]));

		glutPostRedisplay();
	}
}

void
_ZoomModel()
{
	gluProject(
		tModelMatrix[12], tModelMatrix[13], tModelMatrix[14], 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdOldCoord[0], &pdOldCoord[1], &pdOldCoord[2]);

	// MOD-BY-LEETEN 10/01/2010-FROM:
		// pdOldCoord[2] -= 0.01 * ((double)(iCursorY - iBeginY)/(double)piViewport[3]);
	// TO:
	pdOldCoord[2] -= dZoomModelSpeed * ((double)(iCursorY - iBeginY)/(double)piViewport[3]);
	// MOD-BY-LEETEN 10/01/2010-END

	gluUnProject(
		pdOldCoord[0], pdOldCoord[1], pdOldCoord[2],
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdNewCoord[0], &pdNewCoord[1], &pdNewCoord[2]);

	for(int i=0; i<3; i++)
		tModelMatrix[12 + i] = pdNewCoord[i];

	glutPostRedisplay();
}

// ADD-BY-LEETEN 10/01/2010-BEGIN
void 
gcbSetFloating(EGcbParameterName eName, float fValue)
{
	switch(eName)
 	{
 	case GCB_ROTATE_MODEL_SPEED:
 		dRotateCameraSpeed = (double)fValue;
 		break;
 
 	case GCB_ZOOM_MODEL_SPEED: 
 		dZoomModelSpeed = (double)fValue;
 		break;
 
 	case GCB_ROTATE_CAMERA_SPEED:
 		dRotateCameraSpeed  = (double)fValue;
 		break;
 	}
}
// ADD-BY-LEETEN 10/01/2010-END

/////////////////////////////////////////////////////////////////////////////////
void
_DisplayCB()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixd(tProjectionMatrix);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(tViewMatrix);

	glMultMatrixd(tModelMatrix);

	if( _DisplayFunc )
		_DisplayFunc();

	CHECK_OPENGL_ERROR("_DisplayCB()", true);
}

void
_ReshapeCB(int w, int h)
{
	if( w & h )
	{
		glViewport(0, 0, w, h);

		glGetIntegerv(GL_VIEWPORT, piViewport);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(30.0f, (float)piViewport[2]/(float)piViewport[3], 0.05f, 20.0f);
		glGetDoublev(GL_PROJECTION_MATRIX, tProjectionMatrix);
	
		if( _ReshapeFunc )
			_ReshapeFunc(w, h);
	}


	CHECK_OPENGL_ERROR("_ReshapeCB()", true);
}

// ADD-BY-LEETEN 01/20/2011-BEGIN
static void 
_OpenMatrix(char *szMatrixFilename)
{
	FILE *fpMatrix;
	fpMatrix = fopen(szMatrixFilename, "rt") ;

	if( !fpMatrix )
	{
		LOG(printf("the matrix file %s cannot be opened.", szMatrixFilename));
		return;
	}

	for(int i = 0,	r = 0; r < 4; r++)
		for(int		c = 0; c < 4; c++, i++)
		{
			float fTemp;
			fscanf(fpMatrix, "%f", &fTemp);
			tModelMatrix[i] = (double)fTemp;
		}

	for(int i = 0,	r = 0; r < 4; r++)
		for(int		c = 0; c < 4; c++, i++)
		{
			float fTemp;
			fscanf(fpMatrix, "%f", &fTemp);
			tViewMatrix[i] = (double)fTemp;
		}

	fclose(fpMatrix);
	LOG(printf("The model and view matrix are loaded from the file %s.", szMatrixFilename));
}

void 
_SaveMatrix(char *szMatrixFilename)
{
	FILE *fpMatrix;
	fpMatrix = fopen(szMatrixFilename, "wt");
	if( !fpMatrix )
	{
		LOG(printf("storing matrix failed."));
		return;
	}
	for(int i = 0,	r = 0; r < 4; r++)
	{
		for(int		c = 0; c < 4; c++, i++)
			fprintf(fpMatrix, "%f ", tModelMatrix[i]);
		fprintf(fpMatrix, "\n");
	}
	for(int i = 0,	r = 0; r < 4; r++)
	{
		for(int		c = 0; c < 4; c++, i++)
			fprintf(fpMatrix, "%f ", tViewMatrix[i]);
		fprintf(fpMatrix, "\n");
	}
	fclose(fpMatrix);
	LOG(printf("The model and view matrix are saved to the file %s.", szMatrixFilename));
}
// ADD-BY-LEETEN 01/20/2011-END

void 
_KeyCB(unsigned char key, int x, int y)
{
	switch(key) {
		// ADD-BY-LEETEN 08/06/2010-BEGIN
		case 'F': case 'f':
			{
				static int piPrevViewport[4]; // the viewport beforethe 
				static bool bIsFullScreened = false;
				bIsFullScreened = !bIsFullScreened;
				if(bIsFullScreened)
				{
					memcpy(piPrevViewport, piViewport, sizeof(piPrevViewport));
					piPrevViewport[0] = glutGet(GLUT_WINDOW_X);
					piPrevViewport[1] = glutGet(GLUT_WINDOW_Y);
					glutFullScreen();
				}
				else
				{
					glutPositionWindow(piPrevViewport[0], piPrevViewport[1]);
					glutReshapeWindow(piPrevViewport[2], piPrevViewport[3]);
				}
			}
			break;
		// ADD-BY-LEETEN 08/06/2010-END

		// ADD-BY-LEETEN 03/24/2008-BEGIN
		case 'x':	case 'X':
			memcpy(tViewMatrix, tInitViewMatrix, sizeof(tViewMatrix));;

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glMultMatrixd(tInitModelMatrix);
			glRotated(('x'==key)?90.0:-90, 0.0, 1.0, 0.0);
			glGetDoublev(GL_MODELVIEW_MATRIX, tModelMatrix);

			glutPostRedisplay();
			break;

		case 'y':	case 'Y':
			memcpy(tViewMatrix, tInitViewMatrix, sizeof(tViewMatrix));;

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glMultMatrixd(tInitModelMatrix);
			glRotated(('y'==key)?90.0:-90, 1.0, 0.0, 0.0);
			glGetDoublev(GL_MODELVIEW_MATRIX, tModelMatrix);

			glutPostRedisplay();
			break;

		case 'z':	case 'Z':
			memcpy(tViewMatrix, tInitViewMatrix, sizeof(tViewMatrix));;

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glMultMatrixd(tInitModelMatrix);
			glRotated(('z'==key)?0.0:180, 0.0, 1.0, 0.0);
			glGetDoublev(GL_MODELVIEW_MATRIX, tModelMatrix);

			glutPostRedisplay();
			break;

		case 'r': case 'R':
			memcpy(tViewMatrix, tInitViewMatrix, sizeof(tViewMatrix));;
			memcpy(tModelMatrix, tInitModelMatrix, sizeof(tModelMatrix));;

			glutPostRedisplay();
			break;

		case 'N': case 'n':
			glutPostRedisplay();
			break;

		// ADD-BY-LEETEN 01/20/2011-BEGIN
		case 'm':	// load the matrix
			_OpenMatrix("matrix.txt");
			glutPostRedisplay();
			break;

		case 'M':	// save the matrix
			_SaveMatrix("matrix.txt");
			break;
		// ADD-BY-LEETEN 01/20/2011-END

		case 27:
			exit(0);
			break;
	}
	if( _KeyboardFunc )
		_KeyboardFunc(key, x, y);
}

void 
_SpecialCB(int skey, int x, int y)
{
	switch(skey) {
		case GLUT_KEY_LEFT:
			break;

		case GLUT_KEY_RIGHT:
			break;

		case GLUT_KEY_UP:
			break;

		case GLUT_KEY_DOWN:
			break;
	}
	if( _SpecialFunc )
		_SpecialFunc(skey, x, y);
}

void
_MotionCB(int x, int y)
{
	// flip the y coordinate
	y = piViewport[3] - y;

	iCursorX = x;
	iCursorY = y;
}

void
_MouseCB(int button, int state, int x, int y)
{
	// flip the y coordinate
	y = piViewport[3] - y; 

	iCursorX = x;
	iCursorY = y;
	iBeginX = x;
	iBeginY = y;
	eMouseButton = button;
	eModifier = glutGetModifiers();
	bMoving = (state == GLUT_DOWN)?true:false;

	double pdCurrentModelCenterOnScreen[3];
	gluProject(
		tModelMatrix[12], tModelMatrix[13], tModelMatrix[14], 
		tViewMatrix, tProjectionMatrix, piViewport, 
		&pdCurrentModelCenterOnScreen[0], &pdCurrentModelCenterOnScreen[1], &pdCurrentModelCenterOnScreen[2]);

	dXDiffOnScreen = pdCurrentModelCenterOnScreen[0] - (double)iCursorX;
	dYDiffOnScreen = pdCurrentModelCenterOnScreen[1] - (double)iCursorY;
}

void
_IdleCB()
{
	if(bMoving)
	{
		switch(eMouseButton) 
		{
		case GLUT_LEFT_BUTTON: // pan
			switch( eModifier & ~GLUT_ACTIVE_ALT )
			{
			case 0:
				_PanModel();
				glutPostRedisplay();	
				break;

			case GLUT_ACTIVE_CTRL: // manipulate the object
				_RotateModel();
				glutPostRedisplay();	
				break;

			case GLUT_ACTIVE_SHIFT:
				_ZoomModel();
				glutPostRedisplay();	
				break;
			} // switch(eModifier)
			break;

		case GLUT_RIGHT_BUTTON:
			switch( eModifier & ~GLUT_ACTIVE_ALT )
			{
			case 0:
				_PanCamera();
				glutPostRedisplay();	
				break;

			case GLUT_ACTIVE_CTRL: // manipulate the object
				_RotateCamera();
				glutPostRedisplay();	
				break;

			case GLUT_ACTIVE_SHIFT:
				_ZoomCamera();
				glutPostRedisplay();	
				break;
			} // switch(eModifier)
			break;
		} // switch(eMouseButton) 
	} // if(bMoving)

	if( _IdleFunc )
		_IdleFunc();
}

void
gcbDisplayFunc(void (*_MyDisplayFunc)())
{
	_DisplayFunc = _MyDisplayFunc;
}

void
gcbReshapeFunc(void (*_MyReshapeFunc)(int, int))
{
	_ReshapeFunc = _MyReshapeFunc;
}

void
gcbKeyboardFunc(void (*_MyKeyboardFunc)(unsigned char, int, int))
{
	_KeyboardFunc = _MyKeyboardFunc;
}

void
gcbSpecialFunc(void (*_MySpecialFunc)(int, int, int))
{
	_SpecialFunc = _MySpecialFunc;
}

void
gcbIdleFunc(void (*_MyIdleFunc)())
{
	_IdleFunc = _MyIdleFunc;
}

// ADD-BY-LEETEN 08/14/2010-BEGIN
void
_PrintUsage()
{
	printf(	"\nUsage:\n"
			"\t[Hotkeys]\n"
			"\tF/f: toggle full screen.\n"
			"\tx/X: view from the -x/+x axis.\n"
			"\ty/Y: view from the -y/+y axis.\n"
			"\tz/Y: view from the -z/+z axis.\n"
			"\tr/R: reset the camera.\n"
			// ADD-BY-LEETEN 01/20/2011-BEGIN
			"\tm: Load the model & view matrices.\n"
			"\tM: Save the model & view matrices.\n"
			// ADD-BY-LEETEN 01/20/2011-END
			"\tESC: quit the application.\n"
			"\n"
			"\t[Mouse]\n"
			"\tLeft bottom: pan the object.\n"
			"\tLeft bottom + CTRL: rotate the object space.\n"
			"\tLeft bottom + SHIFT: pull/push the object from the camera.\n"
			"\tRight bottom: pan the camera.\n"
			"\tRight bottom + CTRL: rotate the view angle.\n"
			"\tRight bottom + SHIFT: zoom in/out the camera.\n"
			"\n"
			);
}
// ADD-BY-LEETEN 08/14/2010-END

void
gcbInit(void (*_InitFunc)(), void (*_QuitFunc)())
{
	// register GLUT callbacks
	glutDisplayFunc(_DisplayCB);
	glutReshapeFunc(_ReshapeCB);
	glutKeyboardFunc(_KeyCB);
	glutSpecialFunc(_SpecialCB);
	glutMotionFunc(_MotionCB);
	glutMouseFunc(_MouseCB);
	glutIdleFunc(_IdleCB);

	// the projection matrix is decided in the _ReshapeCB();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glGetDoublev(GL_MODELVIEW_MATRIX, tViewMatrix);
	memcpy(tInitViewMatrix, tViewMatrix, sizeof(tViewMatrix));
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glGetDoublev(GL_MODELVIEW_MATRIX, tModelMatrix);
	memcpy(tInitModelMatrix, tModelMatrix, sizeof(tModelMatrix));

	if( _InitFunc )
		_InitFunc();

	QuitFunc = _QuitFunc;
	atexit(_Quit);

	CHECK_OPENGL_ERROR("gcbInit()", true);

	// ADD-BY-LEETEN 08/14/2010-BEGIN
	LOG(_PrintUsage());
	// ADD-BY-LEETEN 08/14/2010-END
}

/*

$Log: gcb.cpp,v $
Revision 1.9  2011/01/21 20:54:08  leeten

[01/21/2011]
1. [ADD] Define two functions _OpenMatrix() and _SaveMatrix() to load/store the model/view matrices. The loading/storing can be triggered by the hotkeys 'm' and 'M.'

Revision 1.8  2010/10/01 20:42:51  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/08/15 04:20:44  leeten

[08/15/2010]
1. [ADD] Add a new function _PrintUsage() to show the basic usage of LIBGCB. The function will be called at the end of gcbInit().

Revision 1.2  2010/08/06 20:24:41  leeten

[08/05/2010]
1. [ADD] Add the hotkey' f' to toggle full screen mode.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
