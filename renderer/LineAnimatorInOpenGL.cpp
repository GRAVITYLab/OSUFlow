/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include "LineAnimatorInOpenGL.h"
#include "opengl.h"

void 
CLineAnimatorInOpenGL::_SetInteger(int iParameter, int iValue)
{
	switch(iParameter)
	{
	case PUSH_FRONT_ONE_PARTICLE:
		{
			int iNewIndex = 0;
			if( !liParticles.empty() )
			{
				int iFirstIndex = liParticles.front();
				iNewIndex = iFirstIndex  - iValue;
			}
			liParticles.push_front(iNewIndex);
		}
		break;

	case PUSH_BACK_ONE_PARTICLE:
		{
			int iNewIndex = 0;
			if( !liParticles.empty() )
			{
				int iLastIndex = liParticles.back();
				iNewIndex = iLastIndex + iValue;
			}
			liParticles.push_back(iNewIndex);
		}
		break;

	case MOVE_ALL_PARTICLES:
		for(list<int>::iterator 
				liiParticles = liParticles.begin(); 
			liiParticles != liParticles.end(); 
			liiParticles ++)
				*liiParticles += iValue;
		break;

	#if	0	// MOD-BY-LEETEN 04/16/2010-FROM:
		case DEQUEUE_PARTICLES:
			{
				int iNrOfParticles = liParticles.size();
				for(int i = 0; i < iValue && i < iNrOfParticles; i++)
					liParticles.pop_front();
			}
			break;
	 
		case POP_PARTICLES:
			{
				int iNrOfParticles = liParticles.size();
				for(int i = 0; i < iValue && i < iNrOfParticles; i++)
					liParticles.pop_back();
			}
			break;
	#else	// MOD-BY-LEETEN 04/16/2010-TO:
	case REMOVE_FRONT_PARTICLES:
		{
			int iNrOfParticles = liParticles.size();
			for(int i = 0; i < iValue && i < iNrOfParticles; i++)
				liParticles.pop_front();
		}
		break;
 
	case REMOVE_BACK_PARTICLES:
		{
			int iNrOfParticles = liParticles.size();
			for(int i = 0; i < iValue && i < iNrOfParticles; i++)
				liParticles.pop_back();
		}
		break;
	#endif	// MOD-BY-LEETEN 04/16/2010-END

	case REMOVE_ALL_PARTICLES:
		liParticles.clear();
		break;

	case RESET_ALL_PARTICLES:
		{
			if( !liParticles.empty() )
			{
				int iFirstIndex = liParticles.front();
				for(list<int>::iterator 
						liiParticles = liParticles.begin(); 
					liiParticles != liParticles.end(); 
					liiParticles ++)
					*liiParticles -= iFirstIndex;
			}
		}
		break;
	}

	CLineRendererInOpenGL::_SetInteger(iParameter, iValue);
}

void 
CLineAnimatorInOpenGL::_GetInteger(int iParameter, int* piValue)
{
	switch(iParameter)
	{
	case NR_OF_RENDERED_PARTICLES:
		*piValue = iNrOfRenderedParticles;
		break;

	// ADD-BY-LEETEN 04/16/2010-BEGIN
	case NR_OF_PARTICLES:
		*piValue = liParticles.size();
		break;

	case FRONT_PARTICLE_INDEX:
		*piValue = ( liParticles.empty() )?(-1):liParticles.front();
		break;

	case BACK_PARTICLE_INDEX:
		*piValue = ( liParticles.empty() )?(-1):liParticles.back();
		break;
	// ADD-BY-LEETEN 04/16/2010-END

	default:
		CLineRendererInOpenGL::_GetInteger(iParameter, piValue);
	}
}

#if	0	// DEL-BY-LEETEN 04/15/2010-BEGIN
	void 
	CLineAnimatorInOpenGL::_TraverseLinesBegin(int iNrOfTraces)
	{
		CLineRendererInOpenGL::_TraverseLinesBegin(iNrOfTraces);
		iNrOfRenderedParticles = 0;
	}

	void 
	CLineAnimatorInOpenGL::_TraverseTraceBegin(int iTraceIndex, int iNrOfPoints)
	{
		glBegin(GL_LINES); 
	}
#endif	// DEL-BY-LEETEN 04/15/2010-END

void 
CLineAnimatorInOpenGL::_TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT)
{
	static list<int>::iterator liiParticleIterator;
	static VECTOR3 v3PrevPoint;
	static VECTOR3 v3PrevTangent;
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	static VECTOR4 v4PrevColor;
	VECTOR4 v4Color = cColorScheme.V4GetColor();
	// ADD-BY-LEETEN 07/07/2010-END
	VECTOR3 v3Point(fX, fY, fZ);

	if( 0 == iPointIndex )
	{
		for(liiParticleIterator  = liParticles.begin();		
			liiParticleIterator != liParticles.end() && *liiParticleIterator < 0;
			liiParticleIterator++)
			;
	}

	VECTOR3 v3Tangent = v3Point - v3PrevPoint;
	v3Tangent.Normalize();

	if( liiParticleIterator != liParticles.end() && iPointIndex == *liiParticleIterator )
	{
		#if	0	// MOD-BY-LEETEN 07/05/2010-FROM:
			if( iPointIndex > 0 )
			{
				if( 1 == iPointIndex )
					glTexCoord3fv(&v3Tangent[0]);
				else
					glTexCoord3fv(&v3PrevTangent[0]);
				glVertex3fv(&v3PrevPoint[0]);

				glTexCoord3fv(&v3Tangent[0]);
				glVertex3fv(&v3Point[0]);
			}
		#else	// MOD-BY-LEETEN 07/05/2010-TO:
		if( 1 == iPointIndex )
			pv4TexCoords.push_back(VECTOR4(v3Tangent[0], v3Tangent[1], v3Tangent[2], 1.0));
		else
			pv4TexCoords.push_back(VECTOR4(v3PrevTangent[0], v3PrevTangent[1], v3PrevTangent[2], 1.0));
		pv4Coords.push_back(VECTOR4(v3PrevPoint[0], v3PrevPoint[1], v3PrevPoint[2], 1.0));
		// ADD-BY-LEETEN 07/07/2010-BEGIN
		pv4Colors.push_back(v4PrevColor);
		// ADD-BY-LEETEN 07/07/2010-END

		pv4TexCoords.push_back(VECTOR4(v3Tangent[0], v3Tangent[1], v3Tangent[2], 1.0));
		pv4Coords.push_back(VECTOR4(v3Point[0], v3Point[1], v3Point[2], 1.0));
		// ADD-BY-LEETEN 07/07/2010-BEGIN
		pv4Colors.push_back(v4Color);
		// ADD-BY-LEETEN 07/07/2010-END
		#endif	// MOD-BY-LEETEN 07/05/2010-END
		liiParticleIterator++;
		iNrOfRenderedParticles++;
	}

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	cColorScheme._MoveToNextPoint();

	v4PrevColor = v4Color;
	// ADD-BY-LEETEN 07/07/2010-END
	v3PrevPoint = v3Point;
	v3PrevTangent = v3Tangent;
}

CLineAnimatorInOpenGL::CLineAnimatorInOpenGL(void)
{
}

CLineAnimatorInOpenGL::~CLineAnimatorInOpenGL(void)
{
}

/*

$Log: LineAnimatorInOpenGL.cpp,v $
Revision 1.8  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.7  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.6  2010/07/07 17:04:57  leeten

[07/07/2010]
1. [ADD] Get the color from the new structure CColorScheme and specify the color for each point on all traces.

Revision 1.5  2010/07/05 14:27:19  leeten

[07/05/2010]
1. [ADD] Add structures and code segments to supprt the use of vertex array.

Revision 1.4  2010/04/19 19:37:22  leeten

[04/19/2010]
1. [MOD] Change the names of the preprocessors DEQUEUE_PARTICLES and POP_PARTICLES to REMOVE_FRONT_PARTICLES and REMOVE_BACK_PARTICLES, respecitvely.
2. [ADD] Add the code to return #particles, index of the front particle, and index of the last particle by calling the method _GetInteger with parameter NR_OF_PARTICLES, FRONT_PARTICLE_INDEX and BACK_PARTICLE_INDEX.

Revision 1.3  2010/04/16 17:33:46  leeten

[04/16/2010]
1. [ADD] Add the log.


*/
