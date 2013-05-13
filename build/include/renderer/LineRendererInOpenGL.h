/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once

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


#include "LineRenderer.h"

//! The class to implement the rendering algorithms of field lines in OpenGL.
/*!

This class implements the rendering of fields lines in OpenGL. Different effects, such as halo 
and lighting, can be added to the basic drawing of lines.

*/

class CLineRendererInOpenGL :
	public CLineRenderer
{
protected:
	// id of the lighting texture
	GLuint t2dLighting;

	// name of the display list
	GLuint uLid;

	// ADD-BY-LEETEN 07/05/2010-BEGIN
	// the temp. STL vectors that store the vertex attributes and indices
	vector<VECTOR4> pv4Coords;
	vector<VECTOR4> pv4TexCoords;
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	vector<VECTOR4> pv4Colors;
	// ADD-BY-LEETEN 07/07/2010-END

	// ADD-BY-LEETEN 01/18/2011-BEGIN
	vector<int> pviGlPrimitiveBases;
	vector<int> pviGlPrimitiveLengths;
	// ADD-BY-LEETEN 01/18/2011-END

	// the structure to store the vertex attributes for vertex array
	struct CVertexArray {
		// ADD-BY-LEETEN 08/26/2010-BEGIN
		int iNrOfVertices;
		// ADD-BY-LEETEN 08/26/2010-END

		float *pfCoords;
		float *pfTexCoords;
		// ADD-BY-LEETEN 07/07/2010-BEGIN
		float *pfColors;
		// ADD-BY-LEETEN 07/07/2010-END

		CVertexArray()
		{
			// ADD-BY-LEETEN 08/26/2010-BEGIN
			iNrOfVertices = 0;
			// ADD-BY-LEETEN 08/26/2010-END
			pfCoords = NULL;
			pfTexCoords = NULL;
			// ADD-BY-LEETEN 07/07/2010-BEGIN
			pfColors = NULL;
			// ADD-BY-LEETEN 07/07/2010-END
		}

		~CVertexArray()
		{
			if( pfCoords )
				free(pfCoords);
			pfCoords = NULL;

			if( pfTexCoords )
				free(pfTexCoords);
			pfTexCoords = NULL;

			// ADD-BY-LEETEN 07/07/2010-BEGIN
			if( pfColors )
				free(pfColors);
			pfColors = NULL;
			// ADD-BY-LEETEN 07/07/2010-END
		}
	} cVertexArray;
	// ADD-BY-LEETEN 07/05/2010-END

public:
	// ADD-BY-LEETEN 04/14/2010-BEGIN
	enum EParameter {
		PARAMETER_BASE = CLineRenderer::MAX_NR_OF_PARAMETERS,
		MAX_NR_OF_PARAMETERS,
	};
	// ADD-BY-LEETEN 04/14/2010-END

	//! Update the texture to mimic liine illumination based on the latest lighting parameters.
	/*!
	This function will get the current paramter of OpenGL Light 0 and then
	update a 2D texture, which mimics the illumination of streamlines.
	*/
	virtual void _UpdateLighting();

	// ADD-BY-LEETEN 04/14/2010-BEGIN
	//! Enable the texture to mimic liine illumination.
	/*!
	This function will bind the 2D texture to mimic the illumination of streamlines.
	*/
	virtual void _TurnLightingOn();
	// ADD-BY-LEETEN 04/14/2010-END

	// ADD-BY-LEETEN 09/14/2010-BEGIN
	//! Disable the texture to mimic liine illumination.
	/*!
	This function will reset the texture matrix to disable the illumination of streamlines.
	*/
	virtual void _TurnLightingOff();
	// ADD-BY-LEETEN 09/14/2010-END

	virtual void _Draw();
	virtual void _Update();

	virtual void _TraverseLinesBegin(int iNrOfTraces);
	virtual void _TraverseLinesEnd();

	virtual void _TraverseTraceBegin(int iTraceIndex, int iNrOfPoints);
	// MOD-By-LEETEN 01/20/2011-FROM:
		// virtual void _TraverseTraceEnd();
	// TO:
	virtual void _TraverseTraceEnd(int iTraceIndex);
	// MOD-By-LEETEN 01/20/2011-END

	virtual void _TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT);

	CLineRendererInOpenGL(void);
	virtual ~CLineRendererInOpenGL(void);
};

/*

$Log: LineRendererInOpenGL.h,v $
Revision 1.12  2011/01/20 17:18:16  leeten

[01/19/2010]
1. [ADD] Add a parameter iTraceindex to the method _TraverseTraceEnd().

Revision 1.11  2011/01/19 19:26:10  leeten

[01/19/2010]
1. [ADD] Add two vectors pviGlPrimitiveBases and pviGlPrimitiveLengths to store the base and length of each streamline.

Revision 1.10  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.8  2010/08/26 20:44:00  leeten

[08/26/2010]
1. [ADD] Add a field iNrOfVertices to the struct CVertexArray to record the #points for OpenGL (not the #particles).

Revision 1.7  2010/08/26 20:31:41  leeten

[08/26/2010]
1. [MOD] Modified by Tom Peterka.

Revision 1.6  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.5  2010/07/07 17:35:37  leeten
[07/07/2010]
1. [ADD] Add a new filed pv4Colors to store the colors.

Revision 1.4  2010/07/05 14:27:19  leeten

[07/05/2010]
1. [ADD] Add structures and code segments to supprt the use of vertex array.

Revision 1.3  2010/04/16 17:32:05  leeten

[04/16/2010]
1. [ADD] Refine the comments for doxygen.

Revision 1.2  2010/04/15 04:01:56  leeten

[04/12/2010]
1. [MOD] Change the comment for doxygen.
2. [ADD] Declare the enum as EParameter s.t doxygen can generate cross-ref to the parameters' names.
3. [MOD] Change the value of PARAMETER_BASE.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
