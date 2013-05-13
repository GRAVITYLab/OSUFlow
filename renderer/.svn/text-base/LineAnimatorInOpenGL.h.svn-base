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

#include "LineRendererInOpenGL.h"

//! The class to animate the line segments
/*!
*/

class CLineAnimatorInOpenGL :
	public CLineRendererInOpenGL
{
protected:
	//! This list stores the indices of all particles 
	/*!
	*/
	list<int> liParticles;

	// DEL-BY-LEETEN 04/15/2010-BEGIN
		// int iNrOfRenderedParticles;
	// DEL-BY-LEETEN 04/15/2010-END

public:
	enum EParameter
	{
		PARAMETER_BASE = CLineRendererInOpenGL::MAX_NR_OF_PARAMETERS,

		//! Add one particle to the head of the list
		/*! 
		One integer as the index difference from the first particle. 
		The index of the first paticle will be substracted by the specified value
		as the index of the new-added particles.
		*/
		PUSH_FRONT_ONE_PARTICLE,	

		//! Add one particle to the end of the list
		/*! 
		One integer as the index difference from the latest particle 
		*/
		PUSH_BACK_ONE_PARTICLE,	

		//! Add an offset to the indcies fo all particles
		/*! 
		One integer as the offset
		*/
		MOVE_ALL_PARTICLES,	

		#if	0	// MOD-BY-LEETEN 04/16/2010-FROM:
			//! Remove particles from the beginning of the list
			/*! 
			One integer as the number of particles to be removed
			*/
			DEQUEUE_PARTICLES,	
	 
			//! Remove particles from the end of the list
			/*! 
			One integer as the number of particles to be removed
			*/
			POP_PARTICLES,	
		#else	// MOD-BY-LEETEN 04/16/2010-TO:
		//! Remove particles from the beginning of the list
		/*! 
		One integer as the number of particles to be removed
		*/
		REMOVE_FRONT_PARTICLES,	
 
		//! Remove particles from the end of the list
		/*! 
		One integer as the number of particles to be removed
		*/
		REMOVE_BACK_PARTICLES,	
		#endif	// MOD-BY-LEETEN 04/16/2010-END

		//! Remove all particles
		/*!
		*/
		REMOVE_ALL_PARTICLES,

		//! Reset the indices of all particles s.t. the index of the first particle is 0
		/*!
		*/
		RESET_ALL_PARTICLES,

		#if	0	// DEL-BY-LEETEN 04/15/2010-BEGIN
			//! Get the number of rendered particles
			/*!
			One pointer to an integer should be specified to return the number of rendered particles.
			*/
			NR_OF_RENDERED_PARTICLES,
		#endif	// DEL-BY-LEETEN 04/15/2010-END

		// ADD-BY-LEETEN 04/16/2010-BEGIN
		//! Number of particles
		/*!
		One integer as the current number of particles.
		*/
		NR_OF_PARTICLES,

		//! Index of the frontmost (earliest) particle
		/*!
		One integer as the index of the frontmost (earliest) particle
		*/
		FRONT_PARTICLE_INDEX,

		//! Index of the backmost (latest) particle
		/*!
		One integer as the index of the backmost (latest) particle
		*/
		BACK_PARTICLE_INDEX,

		// ADD-BY-LEETEN 04/16/2010-END

		MAX_NR_OF_PARAMETERS
	};

	virtual void _SetInteger(int iParameter, int iValue);
	virtual void _GetInteger(int iParameter, int* piValue);

	#if	0	// DEL-BY-LEETEN 04/15/2010-BEGIN
		virtual void _TraverseLinesBegin(int iNrOfTraces);
		virtual void _TraverseTraceBegin(int iTraceIndex, int iNrOfPoints);
	#endif	// DEL-BY-LEETEN 04/15/2010-END
	virtual void _TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT);

	CLineAnimatorInOpenGL(void);
	virtual ~CLineAnimatorInOpenGL(void);
};

/*

$Log: LineAnimatorInOpenGL.h,v $
Revision 1.6  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.5  2010/08/26 20:30:54  leeten

[08/26/2010]
1. [MOD] Modified by Tom Peterka.

Revision 1.4  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.3  2010/04/19 19:39:11  leeten

[04/19/2010]
1. [MOD] Change the names of the preprocessors DEQUEUE_PARTICLES and POP_PARTICLES to REMOVE_FRONT_PARTICLES and REMOVE_BACK_PARTICLES, respecitvely.
2. [ADD] Add new paramters: NR_OF_PARTICLES, FRONT_PARTICLE_INDEX and BACK_PARTICLE_INDEX.

Revision 1.2  2010/04/16 17:29:47  leeten

[04/16/2010]
1. [DEL] Remove the methods _TraverseLinesBegin and _TraverseTraceBegin.


*/
