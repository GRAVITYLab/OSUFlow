/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#ifndef __OPENGL__H__
#define __OPENGL__H__

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

#ifdef LINUX
#include <GL/gl.h>
#include <GL/glu.h>
#endif

// ADD-BY-LEETEN 08/26/2010-BEGIN
#ifdef WIN32
#include <GL/glut.h>
#endif
// ADD-BY-LEETEN 08/26/2010-END

#define CHECK_OPENGL_ERROR(prefix, terminate)	\
		{\
			GLint iError = glGetError();	\
			if( iError )	\
			{\
				fprintf(stderr, "%s: %s\n", prefix, gluErrorString(iError));	\
				if( terminate )	exit(-iError);	\
			}\
		}

#define CREATE_1D_TEXTURE(TARGET, TID, FILTER, INT_FORMAT, WIDTH, FORMAT, TYPE, PTR)	\
		{\
		glGenTextures(1, &(TID));	\
		glBindTexture(TARGET, TID);	\
		glTexParameteri(TARGET, GL_TEXTURE_WRAP_S, GL_CLAMP);	\
		glTexParameteri(TARGET, GL_TEXTURE_MAG_FILTER, FILTER);	\
		glTexParameteri(TARGET, GL_TEXTURE_MIN_FILTER, FILTER);	\
		glTexImage1D(TARGET, 0, INT_FORMAT,	\
			(WIDTH), 0, (FORMAT), (TYPE), (PTR));	\
		}

#define CREATE_2D_TEXTURE(TARGET, TID, FILTER, INT_FORMAT, WIDTH, HEIGHT, FORMAT, TYPE, PTR)	\
		{\
		glGenTextures(1, &(TID));	\
		glBindTexture(TARGET, TID);	\
		glTexParameteri(TARGET, GL_TEXTURE_WRAP_S, GL_CLAMP);	\
		glTexParameteri(TARGET, GL_TEXTURE_WRAP_T, GL_CLAMP);	\
		glTexParameteri(TARGET, GL_TEXTURE_MAG_FILTER, FILTER);	\
		glTexParameteri(TARGET, GL_TEXTURE_MIN_FILTER, FILTER);	\
		glTexImage2D(TARGET, 0, INT_FORMAT,	\
			(WIDTH), (HEIGHT), 0, (FORMAT), (TYPE), (PTR));	\
		}


#endif	// #ifndef __OPENGL__H__

/*

$Log: opengl.h,v $
Revision 1.6  2010/12/06 15:26:33  leeten
no message

Revision 1.4  2010/08/26 20:44:39  leeten

[08/26/2010]
1. [ADD] Add the header <GL/glut.h> when the preprocessor WIN32 is defined.

Revision 1.3  2010/08/26 20:32:01  leeten

[08/26/2010]
1. [MOD] Modified by Tom Peterka.

Revision 1.2  2010/08/15 12:34:32  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.1  2010/04/15 03:55:12  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
