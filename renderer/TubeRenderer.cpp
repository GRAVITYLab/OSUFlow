/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "TubeRenderer.h"

void 
CTubeRenderer::_SetInteger(int iParameter,	int iValue)
{
	switch(iParameter)
	{
	case DRAW_POLYGON:
		this->iDrawPolygon = iValue;
		break;

	case NR_OF_PATCHES:
		this->iNrOfPatches = iNrOfPatches;
		break;
	}
	CLineRenderer::_SetInteger(iParameter, iValue);
}

void 
CTubeRenderer::_GetInteger(int iParameter,	int* piValue)
{
	switch(iParameter)
	{
	case DRAW_POLYGON:
		*piValue = iDrawPolygon;
		break;

	case NR_OF_PATCHES:
		*piValue = this->iNrOfPatches;
		break;

	default:
		CLineRenderer::_GetInteger(iParameter, piValue);
	}
}

void 
CTubeRenderer::_TraverseLinesBegin(int iNrOfTraces)
{
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	vviTracePatchBases.resize(iNrOfTraces);
	vviTracePatchLengths.resize(iNrOfTraces);
	// ADD-BY-LEETEN 01/20/2011-END

	vv3Coords.clear();
	vv3Normals.clear();
	viPatchPointIndices.clear();
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	vv4Colors.clear();
	cColorScheme._Reset();
	// ADD-BY-LEETEN 07/07/2010-END
}

void 
CTubeRenderer::_TraverseLinesEnd()
{
}

void 
CTubeRenderer::_TraverseTraceBegin(int iTraceIndex, int iNrOfPoints)
{
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	vviTracePatchBases[iTraceIndex] = viPatchPointIndices.size();
	// ADD-BY-LEETEN 01/20/2011-END
}

void 
// MOD-By-LEETEN 01/20/2011-FROM:
	// CTubeRenderer::_TraverseTraceEnd()
// TO:
CTubeRenderer::_TraverseTraceEnd(int iTraceIndex)
// MOD-By-LEETEN 01/20/2011-END
{
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	vviTracePatchLengths[iTraceIndex] = viPatchPointIndices.size() - vviTracePatchBases[iTraceIndex];
	// ADD-BY-LEETEN 01/20/2011-END

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	cColorScheme._MoveToNextTrace();
	// ADD-BY-LEETEN 07/07/2010-END
}

void 
CTubeRenderer::_TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT)
{
	static VECTOR3 v3PrevPoint;
	VECTOR3 v3Point;
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	static VECTOR4 v4PrevColor;
	VECTOR4 v4Color = cColorScheme.V4GetColor();
	// ADD-BY-LEETEN 07/07/2010-END
	v3Point.Set(fX, fY, fZ);
	if( 0 < iPointIndex )
	{
		// compute the tangent
		VECTOR3 v3Tangent = v3Point - v3PrevPoint;
		v3Tangent.Normalize();

		// compute an up vector that is orthogonal to the tangent
		int iMaxDir = 0;
		for(int i = 1; i < 3; i++)
			if( fabsf(v3Tangent[iMaxDir]) > fabsf(v3Tangent[i]) )
				iMaxDir = i;

		VECTOR3 v3Normal = VECTOR3(-v3Tangent[1], v3Tangent[0], 0.0f);
		v3Normal.Normalize();

		// find another vector which is orthogonal to both the tangent and the up vector
		VECTOR3 v3Up = VECTOR3(
			v3Normal[1] * v3Tangent[2] - v3Normal[2] * v3Tangent[1],
			v3Normal[2] * v3Tangent[0] - v3Normal[0] * v3Tangent[2],
			v3Normal[0] * v3Tangent[1] - v3Normal[1] * v3Tangent[0]
		);
		v3Up.Normalize();

		int iNrOfIterations = (1 == iPointIndex)?2:1;
		for(int i = 0; i < iNrOfIterations; i++)
		{
			for(int p = 0; p < iNrOfPatches; p++)
			{
				float fAngle = (float) p * 2.0f * (float)M_PI / (float)iNrOfPatches;

				VECTOR3 v3PointNormal = 
					VECTOR3(v3Normal * cosf(fAngle) + v3Up * sinf(fAngle));
				v3PointNormal.Normalize();

				VECTOR3 v3PointOffset = v3PointNormal;
				v3PointOffset.scale(cLine.fWidth/2.0f);

				VECTOR3 v3Coord;
				if(1 == iPointIndex && i == 0 )
					v3Coord = v3PrevPoint + v3PointOffset;
				else
					v3Coord = v3Point + v3PointOffset;

				// store the cooridnate to the list
				vv3Coords.push_back(v3Coord);

				// store the normal
				vv3Normals.push_back(v3PointNormal);

				// ADD-BY-LEETEN 07/07/2010-BEGIN
				if(1 == iPointIndex && i == 0 )
					vv4Colors.push_back(v4PrevColor);
				else
					vv4Colors.push_back(v4Color);
				// ADD-BY-LEETEN 07/07/2010-END

				// store the line index of this point
				viPointTraceIndices.push_back(iTraceIndex);
			} // for p
		} // for i

		for(int p = 0; p < iNrOfPatches; p++)
		{
			int iNextP = (p + 1) % iNrOfPatches;
			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches + p);
			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches - iNrOfPatches + iNextP);
			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches - iNrOfPatches + p );

			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches + p);
			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches + iNextP);
			viPatchPointIndices.push_back(vv3Coords.size()-iNrOfPatches - iNrOfPatches + iNextP );
		}
	}
	v3PrevPoint = v3Point;
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	v4PrevColor = v4Color;
	cColorScheme._MoveToNextPoint();
	// ADD-BY-LEETEN 07/07/2010-END
}

void
CTubeRenderer::_Update()
{
	if( !this->pDataSource )
		return;

	_TraverseLines();
}

CTubeRenderer::CTubeRenderer(void)
{
	this->iNrOfPatches = 6;
}

CTubeRenderer::~CTubeRenderer(void)
{
}

/*

$Log: TubeRenderer.cpp,v $
Revision 1.5  2011/01/20 17:14:29  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [ADD] Compute the point base per streamline in the method _TraverseTraceBegin().
3. [MOD] Compute the #points per streamlines in the new _TraverseTraceEnd().

Revision 1.4  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.3  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.2  2010/07/07 17:45:15  leeten
[07/07/2010]
1. [ADD] Get the user-specified color through the methods of CColorScheme.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
