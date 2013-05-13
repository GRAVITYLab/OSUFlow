/*

gcbTemplate:
The is a demo to show how to use GCB to quick create an OpenGL program with GLUT.

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "gcb.h"

void
_DisplayFunc()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render the scene
	glutWireCube(1.0);

	// NOTE: Call glutSwapBuffers() at the end of your display function
	glutSwapBuffers();
}

void
init()
{
	LOG(printf("Initialize OpenGL here."));
	glEnable(GL_DEPTH_TEST);
}

void 
quit()
{
	LOG(printf("Clean up here."));
}

int
main(int argn, char* argv[])
{
	// when use GCB, it is still needed to initialize GLUT
	glutInit(&argn, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL );
	glutCreateWindow("GCB Template");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbTemplate.cpp,v $
Revision 1.3  2010/10/01 20:38:23  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.2  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
