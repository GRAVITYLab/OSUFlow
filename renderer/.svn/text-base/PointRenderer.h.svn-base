#pragma once

#define GLUT_BUILDING_LIB

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#endif

#ifdef LINUX
#include <GL/glut.h> 
#endif

#include "Renderer.h"

class CPointRenderer :
	public CRenderer
{
protected:
	// name of the display list
	GLuint uPid;

	// Data Structure to store all vertices and their attributes
	float *pointArray;
	float *colorArray;
	float *normalArray;

	// Number of particles in the system
	int particlecount;

	// Maximum timestep of a particle
	int maxT;

	// Number of traces
	int numT;

	MATRIX4 modelview;
	const list<vtListSeedTrace*>* sl_list;
	int currentTrace, currentTraceOffset, particlecountCurrentTrace;
	int singletracemodeEnabled;

	// Data Structure to store quads - one for each glyph

public:
	virtual void _Draw();

	virtual void _Update();
	virtual void _UpdateTraceMode();
	virtual void _UpdateSelectedTrace();

	virtual void _InitHeadlight();

	virtual void _TraversePoints_Point();

	virtual void _RenderPoints_Point();

	virtual void _TraversePointsBegin();
	virtual void _TraversePointsEnd();

	virtual void _ZMinMaxEyeSpace( float *z );
	virtual int _CountParticles( int *maxT, int *numT );
	virtual void _CurrentModelview();

	CPointRenderer(void);
	~CPointRenderer(void);
};
