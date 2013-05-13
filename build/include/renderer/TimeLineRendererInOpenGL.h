/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once
#include "LineRendererInOpenGL.h"

//! The class to render and animate the path lines
/*!
*/

class CTimeLineRendererInOpenGL :
	public CLineRendererInOpenGL
{
protected:
	//! A flag whether the time varying mode is on.
	/*!
	*/
	bool bIsTimeVaryingOn;

	//! An floating point as the min. time step of the time interval of interest. 
	/*!
	*/
	float fMinTimeStep;

	//! An floating point as the max. time step of the time interval of interest. 
	/*!
	*/
	float fMaxTimeStep;

public:
	enum EParameter {
		PARAMETER_BASE = CLineRendererInOpenGL::MAX_NR_OF_PARAMETERS,

		//! The ON/OFF of the time varying mode.
		/*!
		An integer as the flag whether the time varying mode is on/off.
		\sa bIsTimeVaryingOn
		*/
		TIME_VARYING_ON,

		//! The min. time step of the time interval of interest.
		/*!
		A floating point number as the minimal time step of interest.
		\sa fMinTimeStep
		*/
		MIN_TIME_STEP,

		//! The max. time step of the time interval of interest.
		/*!
		A floating point number as the maximal time step of interest.
		\sa fMaxTimeStep
		*/
		MAX_TIME_STEP,

		MAX_NR_OF_PARAMETERS,
	} ;

	virtual void _SetInteger(int iParameter,	int iValue);
	virtual void _GetInteger(int iParameter,	int *piValue);
	virtual void _SetFloat(int iParameter,		float fValue);
	virtual void _GetFloat(int iParameter,		float *pfValue);

	virtual void _TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT);
	virtual void _TraverseLines();

	CTimeLineRendererInOpenGL(void);
	virtual ~CTimeLineRendererInOpenGL(void);
};

/*

$Log: TimeLineRendererInOpenGL.h,v $
Revision 1.5  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/07/07 17:50:15  leeten

[07/07/2010]
1. [ADD] Add the section for CVS log.


*/
