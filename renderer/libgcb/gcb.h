/*

LIBGCB: 
It is a library that provides a template to design a simple OpenGL program with a single display window and basic UI.

Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May,2010

*/

#ifndef __GLUT_CALLBACK_H__
#define __GLUT_CALLBACK_H__

#include <iostream>
using namespace std;

#if	0	// MOD-BY-LEETEN 08/06/2010-FROM:
	#define GLUT_BUILDING_LIB
	#include <GL/glut.h>
#else	// MOD-BY-LEETEN 08/06/2010-TO:

// DEL-BY-LEETEN 08/13/2010-BEGIN
	// #define USE_FREEGLUT
// DEL-BY-LEETEN 08/13/2010-END

#ifdef USE_FREEGLUT
	#define GLUT_DISABLE_ATEXIT_HACK
	#include <GL/freeglut.h>
	// MOD-BY-LEETEN 08/13/2010-FROM:
		// #pragma comment (lib, "freeglut.lib")      /* link with Windows MultiMedia lib */
	// TO:
	#ifdef  WIN32
		#pragma comment (lib, "freeglut.lib")      /* link with Windows MultiMedia lib */
	#endif
	// MOD-BY-LEETEN 08/13/2010-END
#else
	#define GLUT_BUILDING_LIB

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#endif

#ifdef LINUX
#include <GL/glut.h>
#include <GL/gl.h>
#endif

// ADD-BY-LEETEN 10/01/2010-BEGIN
#ifdef WIN32
#include <GL/glut.h>
#include <GL/gl.h>
#endif
// ADD-BY-LEETEN 10/01/2010-END

	// MOD-BY-LEETEN 08/13/2010-FROM:
		// #pragma comment (lib, "glut32.lib")      /* link with Windows MultiMedia lib */
	// TO:
	#if	0	// DEL-BY-LEETEN 08/23/2012-BEGIN
	#ifdef  WIN32
		#pragma comment (lib, "glut32.lib")      /* link with Windows MultiMedia lib */
	#endif
	#endif		// DEL-BY-LEETEN 08/23/2012-END
	// MOD-BY-LEETEN 08/13/2010-END

#endif

#endif	// MOD-BY-LEETEN 08/06/2010-END

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define	_USE_MATH_DEFINES
#include <math.h>

// my own useful macros
// ADD-BY-LEETEN 08/13/2010-BEGIN
#ifdef WIN32
// ADD-BY-LEETEN 08/13/2010-END

#define LOG(printf_stmt)	\
	{						\
		cout<<"["<<__FUNCTION__<<"() @ "<<strrchr(__FILE__, '\\') + 1<<"("<<__LINE__<<")]: ";		\
		(printf_stmt);		\
		cout<<endl;			\
	}						
// ADD-BY-LEETEN 08/13/2010-BEGIN
#else

#define LOG(printf_stmt)	\
	{						\
		cout<<"["<<__FUNCTION__<<"() @ "<<__FILE__<<"("<<__LINE__<<")]: ";		\
		(printf_stmt);		\
		cout<<endl;			\
	}						

#endif
// ADD-BY-LEETEN 08/13/2010-END

#define CHECK_OPENGL_ERROR(prefix, terminate)	\
		{\
			GLint iError = glGetError();	\
			if( iError )	\
			{\
				fprintf(stderr, "%s: %s\n", prefix, gluErrorString(iError));	\
				if( terminate )	exit(-iError);	\
			}\
		}


#ifdef NDEBUG
	#undef assert
	#define assert(expr) \
	{\
		if( !(expr) ) \
		{	\
		fprintf(stderr, "%s in %s (%d):", #expr, __FILE__, __LINE__);\
			perror("");	\
			exit(-1);	\
		}\
	}\

#else
	#include <assert.h>
#endif

// ADD-BY-LEETEN 10/01/2010-BEGIN
enum EGcbParameterName{
	GCB_PARAMETER_BASE = 0x0100,
	GCB_ROTATE_MODEL_SPEED, 
	GCB_ZOOM_MODEL_SPEED, 
	GCB_ROTATE_CAMERA_SPEED, 
	GCB_ZOOM_CAMERA_SPEED, 
};
void gcbSetFloating(EGcbParameterName eName, float fValue);
// ADD-BY-LEETEN 10/01/2010-END

//! Specify the content of the display callback
//!
void gcbDisplayFunc(void (*_MyDisplayFunc)());

//! Specify the content of the reshape callback
//!
void gcbReshapeFunc(void (*_MyReshapeFunc)(int, int));

//! Specify the content of the keyboard callback
//!
void gcbKeyboardFunc(void (*_MyKeyboardFunc)(unsigned char, int, int));

//! Specify the content of the special key callback
//!
void gcbSpecialFunc(void (*_MySpecialFunc)(int, int, int));

//! Specify the content of the idle callback
//!
void gcbIdleFunc(void (*_MyIdleFunc)());

//! Specify two callbacks _InitFunc and _QuitFunc. 
//! The function _InitFunc will be called once the OpenGL window has been created. 
//! The function _QuitFunc will be called when the application is going to terminiate.
//!
void gcbInit(void (*_InitFunc)() = NULL, void (*_QuitFunc)() = NULL);

#endif	// __GLUT_CALLBACK_H__

/*

$Log: gcb.h,v $
Revision 1.8  2010/10/01 20:42:51  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/08/14 22:58:15  leeten

[08/14/2010]
1. [DEL] Remove the definition of the preprocessor USE_FREEGULT.
2. [MOD] Change the definition of LOG for non-win32 environment.

Revision 1.2  2010/08/06 20:24:57  leeten

[08/06/2010]
1. [ADD] Use freeglut when USE_FREEGLUT is defined.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
