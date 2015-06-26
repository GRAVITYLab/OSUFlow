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

// linear interpolation
inline float Lerp(float x, float y, float ratio)
{
    return (x * (1 - ratio) + y * ratio);
}

// bilinear interpolation
inline float BiLerp(float ll, float hl, float lh, float hh, float coeff[2])
{
    return (Lerp(Lerp(ll, hl, coeff[0]), Lerp(lh, hh, coeff[0]), coeff[1]));
}

// trilinear interpolation
inline float TriLerp(float lll, float hll, float lhl, float hhl,
              float llh, float hlh, float lhh, float hhh,
              float coeff[3])
{
    return (Lerp(BiLerp(lll, hll, lhl, hhl, coeff),
                 BiLerp(llh, hlh, lhh, hhh, coeff),
                 coeff[2]));
}

// barycentric interpolation
inline float BaryInterp(float dataValue[4], float coeff[3])
{
    return (dataValue[0] +
            (dataValue[1]-dataValue[0])*coeff[0] +
            (dataValue[2]-dataValue[0])*coeff[1] +
            (dataValue[3]-dataValue[0])*coeff[2]);
}

// Gaussian smoothing filter
void operateGaussianLPF(int width, int height, int element, float *pData);

#endif
