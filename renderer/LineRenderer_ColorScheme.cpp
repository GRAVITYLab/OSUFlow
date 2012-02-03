/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "LineRenderer.h"

// ADD-BY-LEETEN 02/03/2012-BEGIN
//! Specify the default line color
void
CLineRenderer::CColorScheme::_SetColor(const VECTOR4& v4Color)
{
	this->v4Color = v4Color;
}
// ADD-BY-LEETEN 02/03/2012-END

VECTOR4 
CLineRenderer::CColorScheme::V4GetColor()
{
	VECTOR4 v4Color;
	switch(iScheme )
	{
	case COLOR_PER_POINT:
	case COLOR_PER_TRACE:
		if( NULL != plv4Colors )
		{
			v4Color = *liv4Colors;
			break;
		}
		else
		{
			// warning
		}
		// ADD-BY-LEETEN 02/03/2012-BEGIN
		break;
		// ADD-BY-LEETEN 02/03/2012-END

	case COLOR_ALL_WHITE:
		// MOD-BY-LEETEN 02/03/2012-FROM:
			// v4Color = VECTOR4(1.0f, 1.0f, 1.0f, 1.0f);
		// TO:
		v4Color = this->v4Color;
		// MOD-BY-LEETEN 02/03/2012-END
		break;
	}
	return v4Color;
}

void
CLineRenderer::CColorScheme::_Reset()
{
	switch(iScheme)
	{
	case COLOR_PER_TRACE:
	case COLOR_PER_POINT:
		if( NULL == plv4Colors )
		{
			// warning
			return;
		}
	
		liv4Colors = plv4Colors->begin();
		break;
	}
}

void
CLineRenderer::CColorScheme::_MoveToNextPoint()
{
	switch(iScheme)
	{
	case COLOR_PER_POINT:
		if( NULL == plv4Colors )
		{
			// warning
			return;
		}

		if( liv4Colors == plv4Colors->end() )
		{
			// warning
			return;
		}

		liv4Colors++;
		break;
	}
}

void
CLineRenderer::CColorScheme::_MoveToNextTrace()
{
	switch(iScheme)
	{
	case COLOR_PER_TRACE:
		if( NULL == plv4Colors )
		{
			// warning
			return;
		}

		if( liv4Colors == plv4Colors->end() )
		{
			// warning
			return;
		}

		liv4Colors++;
		break;
	}
}

CLineRenderer::CColorScheme::CColorScheme()
{
	iScheme = COLOR_ALL_WHITE;
	plv4Colors = NULL;
	// ADD-BY-LEETEN 02/03/2012-BEGIN
	v4Color.Set(1.0f, 1.0f, 1.0f, 1.0f);
	// ADD-BY-LEETEN 02/03/2012-END
}

/*

$Log: LineRenderer_ColorScheme.cpp,v $
Revision 1.4  2011-02-07 02:57:41  leeten

[02/06/2011]
1. [ADD] When the coloring scheme is COLOR_ALL_WHITE, use the color specified in the vector v4Color. Define the method to _SetColor to specify this color.

Revision 1.3  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.2  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.1  2010/07/07 17:30:15  leeten
[07/07/2010]
1. [1ST] First time checkin. This file defines the methods of CLineRenderer::CColorScheme.


*/
