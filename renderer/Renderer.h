/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#pragma once

#define	_USE_MATH_DEFINES
#include <math.h>

#include "OSUFlow.h"

//! The base class of the rendering algorithms for OSUFlow.
/*!

It includes methods to set and get the parameter's values, 
methods to specify the input data structure a method to update 
and to draw the geometric primtives for rendering.

*/

class CRenderer
{
protected:
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	const void *pColorSource;
	// ADD-BY-LEETEN 07/07/2010-END

	const void *pDataSource;

	//! The bounding box of the underline vector field.
	/*!
	Contains two 3D vectors to define left/right/top/bottom/front/back coordinates of the bounding box.
	*/
	struct CBoundingBox 
	{
		VECTOR3 pv3Corners[2];
	} cBoundingBox;

	// ADD-BY-LEETEN 08/23/2012-BEGIN
	//! A flag to indicate whether the bounding box is given
	int iIsWithBoundingBox;
	// ADD-BY-LEETEN 08/23/2012-END
	
public:
	// ADD-BY-LEETEN 04/14/2010-BEGIN
	//! Names of the parameters
	/*!
	For the following classes, the first symbol should be PARAMETER_BASE and the last one should be MAX_NR_OF_PARAMETERS.
	*/
	enum EParameter {
		PARAMETER_BASE = 0x00,
		// ADD-BY-LEETEN 01/18/2011-BEGIN
		DATA_SOURCE_POINTER,
		COLOR_SOURCE_POINTER,
		// ADD-BY-LEETEN 01/18/2011-END
		WITH_BOUNDING_BOX,	// ADD-BY-LEETEN 08/23/2012
		MAX_NR_OF_PARAMETERS,
	} ;
	// ADD-BY-LEETEN 04/14/2010-END

	//! Specify a parameter whose values contain one intergeer
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param iValue value of the parameter

	\sa EParameter
	*/
	virtual void _SetInteger(int iParameter,	int iValue);

	//! Get the values of the specified a parameter, whose type is one intergeer
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param *piValue pointer to store the value of the parameter

	\sa EParameter
	*/
	virtual void _GetInteger(int iParameter,	int *piValue);

	//! Specify a parameter whose values contain one floating-point number
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param fValue value of the parameter

	\sa EParameter
	*/
	virtual void _SetFloat(int iParameter,		float fValue);

	//! Get the values of the specified a parameter, whose type is one floating-point number
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param *pfValue pointer to store the value of the parameter

	\sa EParameter
	*/
	virtual void _GetFloat(int iParameter,		float *pfValue);

	//! Specify a parameter whose values contain an array of integers
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param piValue an array of integers as the values of the parameter

	\sa EParameter
	*/
	virtual void _SetIntegerv(int iParameter,	int iNrOfValues, int piValues[]);

	//! Get the values of the specified a parameter, whose type is an array of integers
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param piValue an array of integers as the values of the parameter

	\sa EParameter
	*/
	virtual void _GetIntegerv(int iParameter,	int iNrOfValues, int piValues[]);

	//! Specify a parameter whose values contain an array of floating-point numbers
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param pfValue an array of integers as the values of the parameter

	\sa EParameter
	*/
	virtual void _SetFloatv(int iParameter,		int iNrOfValues, float pfValues[]);

	//! Get the values of the specified a parameter, whose type is an array of floating-point number
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param pfValue an array of integers as the values of the parameter

	\sa EParameter
	*/
	virtual void _GetFloatv(int iParameter,	int iNrOfValues, float pfValues[]);

	// ADD-BY-LEETEN 01/18/2011-BEGIN
	/*! 
	\param iParameter name of the parameter. The names are defined in the enum EParameter.
	\param _pPtr the input pointer
	*/
	virtual void _SetPointer(int iParameter, const void* _pPtr);
	// ADD-BY-LEETEN 01/18/2011-END

	//! Specify a pointer to the input data
	/*! 
	\param pDataSource a read-only pointer to the input data strcture
	*/
	virtual void _SetDataSource(const void * pDataSource);

	//! Specify the bounding box of the vector field
	/*! 
	\param fLeft	the min. X value of the bounding box
	\param fBottom	the min. Y value of the bounding box
	\param fFront	the min. Z value of the bounding box
	\param fRight	the max. X value of the bounding box
	\param fTop		the max. Y value of the bounding box
	\param fBack	the max. Z value of the bounding box
	*/
	virtual void _SetBoundingBox(
		float fLeft,  float fBottom, float fFront,
		float fRight, float fTop,	 float fBack);

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	//! Specify the colors of all traces/vertices
	virtual void _SetColorSource(
		const void *pColorSource);
	// ADD-BY-LEETEN 07/07/2010-END

	//! Enforce the update of the geometry of the rendering primitives. This method should be called when the input data is changed.
	//!
	virtual void _Update() = 0;	

	//! Draw the rendering primitives, which should be called in the display callbacks.
	//!
	virtual void _Draw() = 0;

	CRenderer(void);
	virtual ~CRenderer(void);
};

/*

$Log: Renderer.h,v $
Revision 1.8  2011/01/19 19:30:15  leeten

[01/19/2010]
1. [ADD] Define the enum DATA_SOURCE_POINTER and COLOR_SOURCE_POINTER.

Revision 1.7  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.6  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.5  2010/07/07 17:53:15  leeten

[07/07/2010]
1. [ADD] Declare a new method _SetColorSource() to specify the buffer of the colors.

Revision 1.4  2010/04/23 19:11:06  leeten

[2010/04/23]
1. [MOD] Change the filename of the headers.

Revision 1.3  2010/04/16 17:32:31  leeten

[04/16/2010]
1. [ADD] Refine the comments for doxygen.

Revision 1.2  2010/04/15 04:02:47  leeten

[04/12/2010]
1. [MOD] Change the comment for doxygen.
2. [ADD] Declare the enum as EParameter s.t doxygen can generate cross-ref to the parameters' names.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
