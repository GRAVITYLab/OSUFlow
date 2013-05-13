/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "TimeLineRendererInOpenGL.h"

void 
CTimeLineRendererInOpenGL::_SetInteger(int iParameter,	int iValue)
{
	switch(iParameter)
	{
	case TIME_VARYING_ON:
		bIsTimeVaryingOn = bool(iValue);
		break;
	}
	CLineRendererInOpenGL::_SetInteger(iParameter, iValue);
}

void 
CTimeLineRendererInOpenGL::_GetInteger(int iParameter,	int* piValue)
{
	switch(iParameter)
	{
	case TIME_VARYING_ON:
		*piValue = int(bIsTimeVaryingOn);
		break;
	default:
		CLineRendererInOpenGL::_GetInteger(iParameter, piValue);
	}
}

void 
CTimeLineRendererInOpenGL::_SetFloat(int iParameter,	float fValue)
{
	switch(iParameter)
	{
	case MIN_TIME_STEP:
		fMinTimeStep = fValue;
		break;
	case MAX_TIME_STEP:
		fMaxTimeStep = fValue;
		break;
	}
	CLineRendererInOpenGL::_SetFloat(iParameter, fValue);
}

void 
CTimeLineRendererInOpenGL::_GetFloat(int iParameter,	float* pfValue)
{
	switch(iParameter)
	{
	case MIN_TIME_STEP:
		*pfValue = fMinTimeStep;
		break;
	case MAX_TIME_STEP:
		*pfValue = fMaxTimeStep;
		break;
	default:
		CLineRendererInOpenGL::_GetFloat(iParameter, pfValue);
	}
}

void 
CTimeLineRendererInOpenGL::_TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT)
{
	static list<int>::iterator liiParticleIterator;
	static VECTOR3 v3PrevPoint;
	static VECTOR3 v3PrevTangent;
	// ADD-BY-LEETEN 04/16/2010-BEGIN
	static float fPrevT;
	// ADD-BY-LEETEN 04/16/2010-END
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	static VECTOR4 v4PrevColor;
	VECTOR4 v4Color = cColorScheme.V4GetColor();
	// ADD-BY-LEETEN 07/07/2010-END
	VECTOR3 v3Point(fX, fY, fZ);

	VECTOR3 v3Tangent = v3Point - v3PrevPoint;
	v3Tangent.Normalize();

	if( 0 < iPointIndex )
	{
		float fT0 = max(fPrevT, fMinTimeStep);
		float fT1 = min(fT, fMaxTimeStep);
		if( fT0 < fT1 )
		{
			float fMinCoeff = (fT0 - fPrevT)/(fT - fPrevT);
			float fMaxCoeff = (fT1 - fPrevT)/(fT - fPrevT);
			if( iPointIndex > 0 )
			{
				if( 1 == iPointIndex )
					pv4TexCoords.push_back(VECTOR4(v3Tangent[0], v3Tangent[1], v3Tangent[2], 1.0));
				else
					pv4TexCoords.push_back(VECTOR4(v3PrevTangent[0], v3PrevTangent[1], v3PrevTangent[2], 1.0));
				VECTOR3 v3PrevP = v3PrevPoint + fMinCoeff * (v3Point - v3PrevPoint);
				pv4Coords.push_back(VECTOR4(v3PrevP[0], v3PrevP[1], v3PrevP[2], 1.0));
				pv4Colors.push_back(v4PrevColor);

				pv4TexCoords.push_back(VECTOR4(v3Tangent[0], v3Tangent[1], v3Tangent[2], 1.0));
				VECTOR3 v3NextP = v3PrevPoint + fMaxCoeff * (v3Point - v3PrevPoint);
				pv4Coords.push_back(VECTOR4(v3NextP[0], v3NextP[1], v3NextP[2], 1.0));
				pv4Colors.push_back(v4Color);
			}
			iNrOfRenderedParticles++;
		}
	}
	fPrevT = fT;

	v3PrevPoint = v3Point;
	v3PrevTangent = v3Tangent;
	v4PrevColor = v4Color;
	cColorScheme._MoveToNextPoint();
}

void 
CTimeLineRendererInOpenGL::_TraverseLines()
{
	iNrOfRenderedParticles = 0;

	const list<vtListTimeSeedTrace*>* sl_list = (const list<vtListTimeSeedTrace*>*)this->pDataSource;

	_TraverseLinesBegin(sl_list->size());

	int iT = 0;
	for(list<vtListTimeSeedTrace*>::const_iterator
			pIter = sl_list->begin(); 
		pIter!=sl_list->end(); 
		pIter++, iT++) 
	{
		const vtListTimeSeedTrace *trace = *pIter; 

		_TraverseTraceBegin(iT, trace->size());

		int iP = 0;
		for(list<VECTOR4*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			VECTOR4 p = **pnIter; 
			_TraversePoint(iP, iT, p[0], p[1], p[2], p[3]); 
		}
		// MOD-By-LEETEN 01/20/2011-FROM:
			// _TraverseTraceEnd();
		// TO:
		_TraverseTraceEnd(iT);
		// MOD-By-LEETEN 01/20/2011-END
	}
	_TraverseLinesEnd();
}

CTimeLineRendererInOpenGL::CTimeLineRendererInOpenGL(void)
{
}

CTimeLineRendererInOpenGL::~CTimeLineRendererInOpenGL(void)
{
}

/*

$Log: TimeLineRendererInOpenGL.cpp,v $
Revision 1.11  2011/01/20 17:13:15  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [MOD] Pass the trace index to the new _TraverseTraceEnd().

Revision 1.10  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.9  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.8  2010/07/09 05:41:09  leeten

[07/09/2010]
1. [ADD] Change the computation of the vertex coordiantes.

Revision 1.7  2010/07/07 17:50:15  leeten

[07/07/2010]
1. [ADD] Add the section for CVS log.


*/
