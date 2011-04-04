/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once

#include <vector>
#include "LineRenderer.h"

//! The class generates 3D tubes for streamlines. 
/*!

This abstract class converts the input streamlines into 3D tubes. The tubes, however, are not converted into 
primitives that can be immediately rendered yet. To render the tubes, a new class should be dervied based on 
this class by defining the method _Draw() to render the tubes.

*/

class CTubeRenderer:
	public CLineRenderer
{
protected:
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	vector<int> vviTracePatchBases;
	vector<int> vviTracePatchLengths;
	// ADD-BY-LEETEN 01/20/2011-END

	//! Number of quad patches per line segment
	/*!
	*/
	int iNrOfPatches;

	//! A vector of 3D vector to store the coordinate of the vertices on the patches
	/*!
	*/
	vector<VECTOR3> vv3Coords;

	//! A vector of 3D vector to store the normal at the vertices on the patches
	/*!
	*/
	vector<VECTOR3> vv3Normals;

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	//! A vector of 3D vector to store the colors at the vertices on the patches
	/*!
	*/
	vector<VECTOR4> vv4Colors;
	// ADD-BY-LEETEN 07/07/2010-END

	//! A vector of integers to store the trace index of each point
	/*!
	*/
	vector<int>		viPointTraceIndices;

	//! A vector of integers to store the index to the point of each patch
	/*!
	*/
	vector<int>		viPatchPointIndices;

	// ADD-BY-LEETEN 10/01/2010-BEGIN
	//! Total #patches from all streamlines
	/*!
	*/
	int iTotalNrOfPatches;
	// ADD-BY-LEETEN 10/01/2010-END

	//! A integer as the drawing mode.
	/*!
	*/
	int iDrawPolygon;
public:
	#if	0	// MOD-BY-LEETEN 04/14/2010-FROM:
		enum {
			PARAMETER_BASE = CLineRenderer::PARAMETER_BASE + 0x10,
	#else	// MOD-BY-LEETEN 04/14/2010-TO:
	enum EParameter {
		PARAMETER_BASE = CLineRenderer::MAX_NR_OF_PARAMETERS,
	#endif	// MOD-BY-LEETEN 04/14/2010-END

		//! Polygon drawing method. 
		/*! 
		One integer from the enum EDrawPolygon should as the drawing styles of the patches. 
		\sa EDrawPolygon
		*/
		DRAW_POLYGON,

		//! The number of quad patches on each line segment
		/*!
		One integer as the number of quad patches per line segment on the tube.
		*/
		NR_OF_PATCHES,	

		// ADD-BY-LEETEN 04/14/2010-BEGIN
		MAX_NR_OF_PARAMETERS,
		// ADD-BY-LEETEN 04/14/2010-END
	};

	// MOD-BY-LEETEN 04/15/2010-FROM:
		// enum{
	// TO:
	enum EDrawPolygon 
	{
	// MOD-BY-LEETEN 04/15/2010-END
		//! Draw the patches as filled polygones
		DRAW_POLYGON_FILL, 

		//! Draw the patches as polygons in wiredframe 
		DRAW_POLYGON_WIREDFRAME
	};

	virtual void _SetInteger(int iParameter,	int iValue);
	virtual void _GetInteger(int iParameter,	int* piValue);

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

	CTubeRenderer(void);
	virtual ~CTubeRenderer(void);
};

/*

$Log: TubeRenderer.h,v $
Revision 1.8  2011/01/20 17:18:46  leeten

[01/19/2010]
1. [ADD] Add a parameter iTraceindex to the method _TraverseTraceEnd().

Revision 1.7  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.5  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.4  2010/07/07 17:47:19  leeten

[07/07/2010]
1. [ADD] Add a new field vv4Colors to store the vertex color during the traversal of streamlines.

Revision 1.3  2010/04/16 17:33:25  leeten

[04/16/2010]
1. [ADD] Refine the comments for doxygen.

Revision 1.2  2010/04/15 04:03:18  leeten

[04/12/2010]
1. [MOD] Change the comment for doxygen.
2. [ADD] Declare the enum as EParameter s.t doxygen can generate cross-ref to the parameters' names.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
