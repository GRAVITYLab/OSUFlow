/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Interpolator
//
///////////////////////////////////////////////////////////////////////////////

#include "Interpolator.h"

#pragma warning(disable : 4251 4100 4244 4101)

/////////////////////////////////////////////////////////////////////////////////////////////////
// Gaussian low-pass filter
//////////////////////////////////////////////////////////////////////////////////////////////////
void operateGaussianLPF(int width, int height, int element, float *pData)
{
	int idx0, idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, iFor, jFor;
	int count;
	float val;

	int gauss_fact[9] = {1, 2, 1, 2, 4, 2, 1, 2, 1};

	for(jFor = 0; jFor < height; jFor++)
	{
		for(iFor = 0; iFor < width; iFor++)
		{
			idx4 = jFor * width + iFor;

			count = 0;
			val = 0;

			//idx0
			if((jFor == (height-1)) || (iFor == 0))
				idx0 = -1;							// invalid value
			else
			{
				idx0 = (jFor+1) * width + (iFor - 1);
				val += gauss_fact[0] * pData[idx0*element];
				count += gauss_fact[0];
			}

			//idx1
			if(jFor == (height-1))
				idx1 = -1;							// invalid value
			else
			{
				idx1 = (jFor+1) * width + iFor;
				val += gauss_fact[1] * pData[idx1*element];
				count += gauss_fact[1];
			}

			//idx2
			if((jFor == (height-1)) || (iFor == (width-1)))
				idx2 = -1;							// invalid value
			else
			{
				idx2 = (jFor+1) * width + (iFor+1);
				val += gauss_fact[2] * pData[idx2*element];
				count += gauss_fact[2];
			}

			//idx3
			if(iFor == 0)
				idx3 = -1;							// invalid value
			else
			{
				idx3 = jFor * width + (iFor - 1);
				val += gauss_fact[3] * pData[idx3*element];
				count += gauss_fact[3];
			}

			val += gauss_fact[4] * pData[idx4*element];
			count += gauss_fact[4];

			//idx5
			if(iFor == (width-1))
				idx5 = -1;							// invalid value
			else
			{
				idx5 = jFor * width + (iFor + 1);
				val += gauss_fact[5] * pData[idx5*element];
				count += gauss_fact[5];
			}

			//idx6
			if((jFor == 0) || (iFor == 0))
				idx6 = -1;							// invalid value
			else
			{
				idx6 = (jFor-1) * width + (iFor - 1);
				val += gauss_fact[6] * pData[idx6*element];
				count += gauss_fact[6];
			}

			//idx7
			if(jFor == 0)
				idx7 = -1;							// invalid value
			else
			{
				idx7 = (jFor-1) * width + iFor;
				val += gauss_fact[7] * pData[idx7*element];
				count += gauss_fact[7];
			}

			//idx8
			if((jFor == 0) || (iFor == (width-1)))
				idx8 = -1;							// invalid value
			else
			{
				idx8 = (jFor-1) * width + (iFor + 1);
				val += gauss_fact[8] * pData[idx8*element];
				count += gauss_fact[8];
			}

			pData[idx4*element] = val/(float)count;
		} // end of width for
	} // end of height for
}

