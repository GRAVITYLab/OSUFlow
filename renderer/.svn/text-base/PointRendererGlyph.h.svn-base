#pragma once

#define GLUT_BUILDING_LIB

#if defined(MAC_OSX_OMPI) || defined(MAC_OSX_MPICH)
#include <GLUT/glut.h> 
#endif

#ifdef LINUX
#include <GL/glut.h> 
#endif

#include "PointRenderer.h"

#define GLYPH_ARROW 1
#define GLYPH_SPHERE 2
#define GLYPH_CONE 3

class CPointRendererGlyph :
	public CPointRenderer
{
protected:
	// name of the display list
	GLuint uPid;

	// Data Structure to store all point
	float *triangleArray;
	float *angleArray;
	float *heightArray;

	// type of the glyph specified by the user
	int currentGlyphType;

	// Current status of lighting
	int lightEnabled;

	// Current status of sorting
	int sortEnabled;

	// Current color scheme
	int currentColorScheme;

public:
	virtual void _Draw();
	virtual void _Update();

	virtual void _UpdateGlyphType( int glyphType );
	virtual void _UpdateColorScheme( int inColorScheme );
	virtual void _UpdateSorting( int sortFlag );
	virtual void _UpdateLighting( int lightFlag );

	virtual void _TraversePoints_Glyph_Arrow();
	virtual void _TraversePoints_Glyph_gluSphere();
	virtual void _TraversePoints_Glyph_gluCone();

	virtual void _RenderPoints_Glyph_Arrow();
	virtual void _RenderPoints_Glyph_gluSphere();
	virtual void _RenderPoints_Glyph_gluCone();

	CPointRendererGlyph(void);
	~CPointRendererGlyph(void);
};
