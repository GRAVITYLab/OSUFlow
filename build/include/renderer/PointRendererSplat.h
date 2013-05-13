#pragma once

#define GLUT_BUILDING_LIB

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#endif

#ifdef LINUX
#include <GL/glut.h> 
#endif

#include "PointRenderer.h"

class CPointRendererSplat :
	public CPointRenderer
{
protected:
	// name of the display list
	GLuint uPid;

	// Data Structure to store all point
	float *triangleArray;

	// Current status of lighting
	int lightEnabled;

	// Current status of sorting
	int sortEnabled;

	// Current color scheme
	int currentColorScheme;

	// Number of particles in the system
	int particlecount;

	// Maximum timestep of a particle
	int maxT;

public:
	virtual void _Draw();
	virtual void _Update();

	virtual void _UpdateColorScheme( int inColorScheme );
	virtual void _UpdateSorting( int sortFlag );
	virtual void _UpdateLighting( int lightFlag );

	virtual void _TraversePoints_Splat();

	virtual void _RenderPoints_Splat();

	CPointRendererSplat(void);
	~CPointRendererSplat(void);
};
