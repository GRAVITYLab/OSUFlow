/*

Creted by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "Renderer.h"

void 
CRenderer::_SetInteger(int iParameter,	int iValue)
{
}

void 
CRenderer::_GetInteger(int iParameter,	int* piValue)
{
}

void 
CRenderer::_SetFloat(int iParameter,		float fValue)
{
}

void 
CRenderer::_GetFloat(int iParameter,		float *pfValue)
{
}

void 
CRenderer::_SetIntegerv(int iParameter,	int iNrOfValues, int piValues[])
{
}

void 
CRenderer::_GetIntegerv(int iParameter,	int iNrOfValues, int piValues[])
{
}

void 
CRenderer::_SetFloatv(int iParameter,		int iNrOfValues, float pfValues[])
{
}

void 
CRenderer::_GetFloatv(int iParameter,		int iNrOfValues, float pfValues[])
{
}

void 
CRenderer::_SetPointer(int iParameter, const void* _pPtr)
{
	switch(iParameter) 
	{
	case DATA_SOURCE_POINTER:
		this->pDataSource = _pPtr;
		break;

	case COLOR_SOURCE_POINTER:
		this->pColorSource = _pPtr;
		break;
	}
}

// ADD-BY-LEETEN 07/07/2010-BEGIN
//! Specify the colors of all traces/vertices
void CRenderer::_SetColorSource
(
	const void* pColorSource
)
{
	this->pColorSource = pColorSource;
}
// ADD-BY-LEETEN 07/07/2010-END

void 
CRenderer::_SetDataSource(const void * pDataSource)
{
	this->pDataSource = pDataSource;
}

void 
CRenderer::_SetBoundingBox(
	float fLeft,  float fBottom, float fFront,
	float fRight, float fTop,	 float fBack)
{
	cBoundingBox.pv3Corners[0][0] = fLeft;
	cBoundingBox.pv3Corners[0][1] = fBottom;
	cBoundingBox.pv3Corners[0][2] = fFront;
	cBoundingBox.pv3Corners[1][0] = fRight;
	cBoundingBox.pv3Corners[1][1] = fTop;
	cBoundingBox.pv3Corners[1][2] = fBack;
}

CRenderer::CRenderer(void)
{
}

CRenderer::~CRenderer(void)
{
}

/*

$Log: Renderer.cpp,v $
Revision 1.5  2011/01/19 19:29:35  leeten

[01/19/2010]
1. [ADD] Define the method _SetPointer().

Revision 1.4  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.3  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.2  2010/07/07 17:43:45  leeten

[07/07/2010]
1. [ADD] Define a new method _SetColorSource() to specify the buffer of the colors.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
