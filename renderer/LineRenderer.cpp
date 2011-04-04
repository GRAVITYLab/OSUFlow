/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "LineRenderer.h"

// ADD-BY-LEETEN 01/19/2011-BEGIN
void 
CLineRenderer::_CheckTrace
(
	int iTrace,
	bool& bIsDrawingTrace
)
{
	bIsDrawingTrace = true;
}

void 
CLineRenderer::_GetTraceColor
(
	int iTrace, 
	float& fR, 
	float& fG, 
	float& fB, 
	float& fA
)
{
	if( 0 == iTrace )
		this->cColorScheme._Reset();
	VECTOR4 v4Color = this->cColorScheme.V4GetColor();
	fR = v4Color[0];
	fG = v4Color[1];
	fB = v4Color[2];
	fA = v4Color[3];
	this->cColorScheme._MoveToNextTrace();
}
// ADD-BY-LEETEN 01/19/2011-END

void 
CLineRenderer::_SetInteger(int iParameter,	int iValue)
{
	switch(iParameter)
	{
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	case COLOR_SCHEME:
		this->cColorScheme.iScheme = iValue;
		break;
	// ADD-BY-LEETEN 07/07/2010-END

	case ENABLE_LIGHTING:
		this->bIsLightingEnabled = (bool)iValue;
		break;

	case ENABLE_HALO:
		this->bIsHaloEnabled = (bool)iValue;
		break;
	}
	CRenderer::_SetInteger(iParameter, iValue);
}

void 
CLineRenderer::_GetInteger(int iParameter,	int* piValue)
{
	switch(iParameter)
	{
	case ENABLE_LIGHTING:
		*piValue = int(this->bIsLightingEnabled);
		break;

	case ENABLE_HALO:
		*piValue = int(this->bIsHaloEnabled);
		break;

	default:
		CRenderer::_GetInteger(iParameter, piValue);
	}
}

void 
CLineRenderer::_SetFloat(int iParameter,		float fValue)
{
	switch(iParameter)
	{
	case LINE_WIDTH:
		this->cLine.fWidth = fValue;
		break;

	case HALO_WIDTH:
		this->cHalo.fWidth = fValue;
		break;
	}
	CRenderer::_SetFloat(iParameter, fValue);
}

void 
CLineRenderer::_GetFloat(int iParameter,		float* pfValue)
{
	switch(iParameter)
	{
	case LINE_WIDTH:
		*pfValue = this->cLine.fWidth;
		break;

	case HALO_WIDTH:
		*pfValue = this->cHalo.fWidth;
		break;

	default:
		CRenderer::_GetFloat(iParameter, pfValue);
	}
}


void 
CLineRenderer::_SetIntegerv(int iParameter,	int iNrOfValues, int piValues[])
{
	switch(iParameter)
	{
	default:
		CRenderer::_SetIntegerv(iParameter, iNrOfValues, piValues);
	}
}

void 
CLineRenderer::_GetIntegerv(int iParameter,	int iNrOfValues, int piValues[])
{
	switch(iParameter)
	{
	default:
		CRenderer::_GetIntegerv(iParameter, iNrOfValues, piValues);
	}
}


void 
CLineRenderer::_SetFloatv(int iParameter,		int iNrOfValues, float pfValues[])
{
	switch(iParameter)
	{
	case LINE_COLOR:
		switch(iNrOfValues)
		{
		case 3:
		case 4:	
			for(int i = 0; i < 4; i++)
				if( i < 3 )
					this->cLine.v4Color[i] = pfValues[i];
				else
					this->cLine.v4Color[i] = ( 4 == iNrOfValues)?pfValues[i]:1.0;

			break;
		default:
			perror("Invalid parameter");
			exit(EXIT_FAILURE);
		}
		break;

	case HALO_COLOR:
		switch(iNrOfValues)
		{
		case 3:	
		case 4:	
			for(int i = 0; i < 4; i++)
				if( i < 3 )
					this->cHalo.v4Color[i] = pfValues[i];
				else
					this->cHalo.v4Color[i] = ( 4 == iNrOfValues)?pfValues[i]:1.0;

			break;

		default:
			assert("Invalid parameter");
			exit(EXIT_FAILURE);
		}
		break;
	}
	CRenderer::_SetFloatv(iParameter, iNrOfValues, pfValues);
}

void 
CLineRenderer::_GetFloatv(int iParameter,		int iNrOfValues, float pfValues[])
{
	switch(iParameter)
	{
	case LINE_COLOR:
		for(int i = 0; i < 4; i++)
			pfValues[i] = this->cLine.v4Color[i];
		break;

	case HALO_COLOR:
		for(int i = 0; i < 4; i++)
			pfValues[i] = this->cHalo.v4Color[i];
		break;

	default:
		CRenderer::_GetFloatv(iParameter, iNrOfValues, pfValues);
	}
}

void 
CLineRenderer::_TraverseLines()
{
	// ADD-BY-LEETEN 04/15/2010-BEGIN
	iNrOfRenderedParticles = 0;
	// ADD-BY-LEETEN 04/15/2010-END

	const list<vtListSeedTrace*>* sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;

	_TraverseLinesBegin(sl_list->size());

	int iT = 0;
	for(list<vtListSeedTrace*>::const_iterator
			pIter = sl_list->begin(); 
		pIter!=sl_list->end(); 
		pIter++, iT++) 
	{
	    const vtListSeedTrace *trace = *pIter; 

		_TraverseTraceBegin(iT, trace->size());

		int iP = 0;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			VECTOR3 p = **pnIter; 
			_TraversePoint(iP, iT, p[0], p[1], p[2], 0.0f); 
		}
		// MOD-By-LEETEN 01/20/2011-FROM:
			// _TraverseTraceEnd();
		// TO:
		_TraverseTraceEnd(iT);
		// MOD-By-LEETEN 01/20/2011-END
	}

	_TraverseLinesEnd();
}

// ADD-BY-LEETEN 07/07/2010-BEGIN
void 
CLineRenderer::_SetColorSource(const void *pColorSource)
{
	CRenderer::_SetColorSource(pColorSource);
	cColorScheme.plv4Colors = (const list<VECTOR4>*)this->pColorSource;
}
// ADD-BY-LEETEN 07/07/2010-END

CLineRenderer::CLineRenderer(void)
{
	cLine = CLine(1.0f, 1.0f, 1.0f, 1.0f, 2.0f);
	cHalo = CLine(0.0f, 0.0f, 0.0f, 1.0f, 4.0f);
}

CLineRenderer::~CLineRenderer(void)
{
}

/*

$Log: LineRenderer.cpp,v $
Revision 1.9  2011/01/20 17:09:04  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [MOD] Pass the trace index to the new _TraverseTraceEnd().

Revision 1.8  2011/01/19 19:28:16  leeten

[01/19/2010]
1. [ADD] Define a new method _CheckTrace() to decide whether a streamline is drawn or not.
2. [ADD] Define a new method _GetTraceColor() to get the color of the specfied trace.

Revision 1.7  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.5  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.4  2010/07/07 17:28:00  leeten
[07/07/2010]
1. [ADD] Add the processing for the new parameter COLOR_SCHEME in _SetInteger().
2. [ADD] Overload the function _SetColorSource().

Revision 1.3  2010/04/16 17:25:22  leeten

[04/16/2010]
1. [ADD] Set the varialbe iNrOfRenderedParticles to 0 at the beginning of _TraverseLines().

Revision 1.2  2010/04/15 03:57:06  leeten

[04/12/2010]
1. [MOD] Call CRenderer::_SetFloat() other than CRenderer::_SetInteger() in CLineRenderer::_SetFloat().

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
