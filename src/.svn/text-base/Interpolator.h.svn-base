/*************************************************************************
*						OSU Flow Vector Field							 *
*																		 *
*																		 *
*	Created:	Han-Wei Shen, Liya Li									 *
*				The Ohio State University								 *
*	Date:		06/2005													 *
*																		 *
*	Interpolator														 *
*************************************************************************/

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "header.h"

enum LerpType {LINEAR_LERP, NEAREST_LERP};

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

// barycentric interpolation
float BaryInterp(float dataValue[4], float coeff[3]);

// trilinear interpolation
float TriLerp(float lll, float hll, float lhl, float hhl, float llh, float hlh, float lhh, float hhh, float coeff[3]);

// bilinear interpolation
float BiLerp(float ll, float hl, float lh, float hh, float coeff[2]);

// linear interpolation
float Lerp(float x, float y, float ratio);

// Gaussian smoothing filter
void operateGaussianLPF(int width, int height, int element, float *pData);

#endif
