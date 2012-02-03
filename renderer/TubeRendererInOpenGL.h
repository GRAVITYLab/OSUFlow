/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once

#include <vector>

#define GLUT_BUILDING_LIB

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#endif

#ifdef LINUX
#include <GL/glut.h> 
#endif

// ADD-BY-LEETEN 08/26/2010-BEGIN
#ifdef WIN32
#include <GL/glut.h>
#endif
// ADD-BY-LEETEN 08/26/2010-END

#include "TubeRenderer.h"

//! The class to implement the rendering algorithms of 3D tubes in OpenGL.
/*!

This class implements the rendering of 3D tubes in OpenGL. 

*/

class CTubeRendererInOpenGL:
	public CTubeRenderer
{
protected:

	// ADD-BY-LEETEN 07/05/2010-BEGIN
	struct CVertexArray {
		float *pfCoords;
		float *pfNormals;
		// ADD-BY-LEETEN 07/07/2010-BEGIN
		float *pfColors;
		// ADD-BY-LEETEN 07/07/2010-END

		unsigned int *puIndices;

		CVertexArray()
		{
			pfCoords = NULL;
			pfNormals = NULL;
			// ADD-BY-LEETEN 07/07/2010-BEGIN
			pfColors = NULL;
			// ADD-BY-LEETEN 07/07/2010-END
			puIndices = NULL;
		}

		~CVertexArray()
		{
			if( pfCoords )		free(pfCoords);		pfCoords	= NULL;
			if( pfNormals )		free(pfNormals);	pfNormals	= NULL;
			// ADD-BY-LEETEN 07/07/2010-BEGIN
			if( pfColors )		free(pfColors);		pfColors	= NULL;
			// ADD-BY-LEETEN 07/07/2010-END
			if( puIndices  )	free(puIndices );	puIndices 	= NULL;
		} 
	} cVertexArray;
	// ADD-BY-LEETEN 07/05/2010-END

	// name of the display list
	GLuint uLid;

	// ADD-BY-LEETEN 02/03/2012-BEGIN
	GLuint lidLighting;
	// ADD-BY-LEETEN 02/03/2012-END

public:
	// ADD-BY-LEETEN 04/14/2010-BEGIN
	enum EParameter {
		PARAMETER_BASE = CTubeRenderer::MAX_NR_OF_PARAMETERS,
		MAX_NR_OF_PARAMETERS,
	} ;
	// ADD-BY-LEETEN 04/14/2010-END

	virtual void _Draw();
	virtual void _Update();
	// ADD-BY-LEETEN 02/03/2012-BEGIN
	virtual void _UpdateLighting();
	virtual void _TurnLightingOn();
	virtual void _TurnLightingOff();
	// ADD-BY-LEETEN 02/03/2012-END

	CTubeRendererInOpenGL(void);
	virtual ~CTubeRendererInOpenGL(void);
};

/*

$Log: TubeRendererInOpenGL.h,v $
Revision 1.9  2011-04-04 20:19:52  leeten

[04/04/2011]
1. [ADD] Declare the methods _TurnLightingOn(), _TurnLightingOff(), and _UpdateLighting() to control the lightingh.
2. [ADD] Add a display list lidLighting.

Revision 1.8  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.7  2010/08/26 20:44:49  leeten

[08/26/2010]
1. [ADD] Add the header <GL/glut.h> when the preprocessor WIN32 is defined.

Revision 1.6  2010/08/26 20:33:11  leeten

[08/26/2010]
1. [MOD] Modified by Tom Peterka.

Revision 1.5  2010/08/15 12:34:32  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.4  2010/07/07 17:56:44  leeten

[07/07/2010]
1. [ADD] Add a pointer pfColors as the array of the vertex color.

Revision 1.3  2010/07/05 14:27:19  leeten

[07/05/2010]
1. [ADD] Add structures and code segments to supprt the use of vertex array.

Revision 1.2  2010/04/15 04:04:01  leeten

[04/12/2010]
1. [MOD] Change the comment for doxygen.
2. [ADD] Declare the enum as EParameter s.t doxygen can generate cross-ref to the parameters' names.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
