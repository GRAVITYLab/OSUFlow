/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once

#include "Renderer.h"

//! The base class of the rendering algorithms of field lines for OSUFlow.
/*!

This class include the colors/widths of the lines and the halo lines, 
and the material of the lines under the lighting.

*/

class CLineRenderer :
	public CRenderer
{
protected:
	// ADD-BY-LEETEN 04/15/2010-BEGIN
	//! The number of particles are actually rendered
	/*!
	*/
	int iNrOfRenderedParticles;
	// ADD-BY-LEETEN 04/15/2010-END

	//! The flag whether the lighting is on
	/*!
	*/
	bool bIsLightingEnabled;

	//! The flag whether the halo effect is on
	/*!
	*/
	bool bIsHaloEnabled;

	//! The strcut to hold the rendering properties of lines, including color and width
	/*!
	*/
	struct CLine {
		VECTOR4 v4Color;
		float fWidth;

		CLine()
		{
			v4Color[0] = v4Color[1] = v4Color[2] = v4Color[3] = 0.0f;
			fWidth = 1.0f;
		}

		CLine(float fR, float fG, float fB, float fA, float fW)
		{
			v4Color[0] = fR;
			v4Color[1] = fG;
			v4Color[2] = fB;
			v4Color[3] = fA;
			fWidth = fW;
		}
	} cHalo, cLine;

	//! The material perperties of lines for illumination
	/*!
	The strcut includdes the ambient, diffuse, and specula components of the Phong illumination model.
	*/
	struct CMaterial {
		VECTOR4 v4Ambient;
		VECTOR4 v4Diffuse;
		VECTOR4 v4Specular;
		float fShininess;
		
		CMaterial()
		{
			v4Ambient = VECTOR4(0.1f, 0.1f, 0.1f, 1.0f);
			v4Diffuse = VECTOR4(0.3f, 0.3f, 0.3f, 1.0f);
			v4Specular = VECTOR4(0.6f, 0.6f, 0.6f, 1.0f);
			fShininess = 4.0f;
		}
	} cMaterial;

public:
	enum  EParameter{
		PARAMETER_BASE = CRenderer::MAX_NR_OF_PARAMETERS,

		//! Lighting on/off
		/*!
		One integer should be specified to control the on (1) or off (0) of the lighting.
		*/
		ENABLE_LIGHTING,	

		//! Halo on/off	
		/*!
		One integer should be specified to control the on (1) or off (0) of the halo effect.
		*/
		ENABLE_HALO,		

		//! Line color
		/*!
		An array of four floating point numbers should be specified as the color of the lines.
		*/
		LINE_COLOR,			

		//! Line width
		/*!
		One floating point number should be specified as the width of the lines.
		*/
		LINE_WIDTH,			

		//! Halo color
		/*!
		An array of four floating point numbers should be specified as the color of the halo lines.
		*/
		HALO_COLOR,			

		//! Halo width
		/*!
		One floating point number should be specified as the width of the halo lines. 
		The width of the actual halo is (HALO_WIDTH - LINE_WIDTH)/2.
		*/
		HALO_WIDTH,

		// ADD-BY-LEETEN 04/15/2010-BEGIN
		//! The number of rendered particles
		/*!
		One pointer to an integer should be specified to return the number of rendered particles.
		*/
		NR_OF_RENDERED_PARTICLES,
		// ADD-BY-LEETEN 04/15/2010-END

		// ADD-BY-LEETEN 07/07/2010-BEGIN
		//! The color scheme
		/*!
		One integer should be specified as the color scheme
		*/
		COLOR_SCHEME,
		// ADD-BY-LEETEN 07/07/2010-END

		// ADD-BY-LEETEN 04/14/2010-BEGIN
		MAX_NR_OF_PARAMETERS,
		// ADD-BY-LEETEN 04/14/2010-END
	};

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	//! The coloring scheme of the streamlines
	/*!
	The strcut includdes the parameter to control the coloring scheme of the streamlines. 
	The streamlines can be specified per point or per streamline (trace). 
	*/
	struct CColorScheme
	{
		enum{
			//! Assign white color to all points on all traces
			COLOR_ALL_WHITE,

			//! Manually assign colors to all points on all traces
			COLOR_PER_POINT,

			//! Manually assign colors to all traces
			COLOR_PER_TRACE,

			// ADD-BY-LEETEN 01/20/2011-BEGIN
			//! Assign the color on the fly 
			COLOR_ON_THE_FLY,
			// ADD-BY-LEETEN 01/20/2011-END
		};

		//! The color scheme
		int iScheme;

		// ADD-BY-LEETEN 02/03/2012-BEGIN
		//! The default line color, which will be used when the color scheme is ALL WHITE
		VECTOR4 v4Color;

		//! Specify the default line color
		void _SetColor(const VECTOR4& v4Color);
		// ADD-BY-LEETEN 02/03/2012-END

		//! The color buffer
		const list<VECTOR4> *plv4Colors;

		//! A poiunter to the color buffer
		list<VECTOR4>::const_iterator liv4Colors;

		//! Reset the pointer to the beginning of the buffer
		void _Reset();

		//! Reset the pointer to the color for next point on the current trace
		void _MoveToNextPoint();

		//! Reset the pointer to the color for next trace
		void _MoveToNextTrace();

		//! Get the color from the pointer
		VECTOR4 V4GetColor();

		//! Constructor
		CColorScheme();
	} cColorScheme;

	// ADD-BY-LEETEN 01/19/2011-BEGIN
	virtual void _CheckTrace(int iTrace, bool& bIsDrawingTrace);
	virtual void _GetTraceColor(int iTrace, float& fR, float& fG, float& fB, float& fA); 
	// ADD-BY-LEETEN 01/19/2011-END

	virtual void _SetColorSource(const void *pColorSource);
	// ADD-BY-LEETEN 07/07/2010-END

	virtual void _SetInteger(int iParameter,	int iValue);
	virtual void _GetInteger(int iParameter,	int *piValue);
	virtual void _SetFloat(int iParameter,		float fValue);
	virtual void _GetFloat(int iParameter,		float *pfValue);
	virtual void _SetIntegerv(int iParameter,	int iNrOfValues, int piValues[]);
	virtual void _GetIntegerv(int iParameter,	int iNrOfValues, int piValues[]);
	virtual void _SetFloatv(int iParameter,		int iNrOfValues, float pfValues[]);
	virtual void _GetFloatv(int iParameter,	int iNrOfValues, float pfValues[]);

	//! Traverse throught the inptu data in order to create graphics primtitives for rendering. 
	//!
	virtual void _TraverseLines();

	//! This interface will be called at the beginning of _TraverseLines()
	/*!
	\param iNrOfTraces The total number of traces in the data struture
	*/
	virtual void _TraverseLinesBegin(int iNrOfTraces) = 0;

	//! This interface will be called at the end of _TraverseLines()
	//!
	virtual void _TraverseLinesEnd() = 0;

	//! This interface will be called when begin to traverse a trace
	/*!
	\param iTraceIndex Index of the current traversing trace. The first trace has index 0
	\param iNrOfPoints The total number of points on the current traversing trace
	*/
	virtual void _TraverseTraceBegin(int iTraceIndex, int iNrOfPoints) = 0;

	//! This interface will be called when a trace has been traversed
	//!
	// MOD-By-LEETEN 01/20/2011-FROM:
		// virtual void _TraverseTraceEnd() = 0;
	// TO:
	virtual void _TraverseTraceEnd(int iTraceIndex) = 0;
	// MOD-By-LEETEN 01/20/2011-END

	//! This interface will be called when traverse a point 
	/*!
	\param iPointIndex index of the current traversing point. The first point has index 0
	\param iTraceIndex index of the current traversing traces. The first point has index 0
	\param fX	X coordinate of the traversing point
	\param fY	Y coordinate of the traversing point
	\param fZ	Z coordinate of the traversing point
	\param fT	Timestamp of the traversing point
	*/
	virtual void _TraversePoint(int iPointIndex, int iTraceIndex, float fX, float fY, float fZ, float fT) = 0;

	virtual void _Draw() = 0;
	virtual void _Update() = 0;

	CLineRenderer(void);
	~CLineRenderer(void);
};

/*

$Log: LineRenderer.h,v $
Revision 1.10  2011-02-07 02:54:55  leeten

[02/06/2011]
1. [ADD] Adda a vector v4Color to CColorScheme and a method _SetColor to specify this color. This color will be used when the coloring scheme is COLOR_ALL_WHITE.

Revision 1.9  2011/01/20 17:17:50  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [ADD] Add a new enum COLOR_ON_THE_FLY as the new coloring scheme.
3. [MOD] Add a parameter iTraceindex to the method _TraverseTraceEnd().

Revision 1.8  2011/01/19 19:28:48  leeten

[01/19/2010]
1. [ADD] Declare a new method _CheckTrace() to decide whether a streamline is drawn or not.
2. [ADD] Declare a new method _GetTraceColor() to get the color of the specfied trace.

Revision 1.7  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.6  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.5  2010/07/07 17:26:45  leeten
[07/07/2010]
1. [ADD] Define a new parameter COLOR_SCHEME.
2. [ADD] Define a new sub structure CColorScheme to control the coloring scheme for the points/traces.

Revision 1.4  2010/04/23 19:11:06  leeten

[2010/04/23]
1. [MOD] Change the filename of the headers.

Revision 1.3  2010/04/16 17:31:28  leeten

[04/16/2010]
1. [ADD] Add a varialbe iNrOfRenderedParticles as the #particles that have been rendered. The value of iNrOfRenderedParticles is associated with NR_OF_RENDERED_PARTICLES
2. [ADD] Refine the comments for doxygen.

Revision 1.2  2010/04/15 04:00:53  leeten

[04/12/2010]
1. [MOD] Change the comment for doxygen.
2. [ADD] Declare the enum as EParameter s.t doxygen can generate cross-ref to the parameters' names.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
