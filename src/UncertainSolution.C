/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li
//                 The Ohio State University
//                 Date:		06/2005
//                 FieldLine
//
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "UncertainSolution.h"

//#pragma warning(disable : 4251 4100 4244 4101)


/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation

                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
              this software for any application provided this
              copyright notice is preserved.

*/
template <class T>
T box_muller(T m, T s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
    T x1, x2, w, y1;
    static T y2;
    static char use_last = 0;

    if (use_last)		        /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * (T)rand()/RAND_MAX - 1.0;
            x2 = 2.0 * (T)rand()/RAND_MAX - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );
}


//////////////////////////////////////////////////////////////////////////
// input
//		pData:		this is a 2D array for storing static or
//					time-varying data
//      pStd:       standard deviation corresponding to pData
//		nodeNum:	number of total nodes
//		timeSteps:	how many time steps, for static data, this is 1
//////////////////////////////////////////////////////////////////////////
// initialzie data is the data is not that large
GaussianSolution::GaussianSolution(VECTOR3** pData, VECTOR3** pStd, int nodeNum, int timeSteps)
{
    assert((pData != NULL) && (pStd != NULL) && (nodeNum > 0) && (timeSteps > 0));

	m_nNodeNum = nodeNum;
	m_nTimeSteps = timeSteps;
	m_MinT = 0; m_MaxT = timeSteps-1;
	m_pDataArray = pData;
    // uncertainty
    m_pStdArray = pStd;
    m_bNormalize = false;
}

GaussianSolution::GaussianSolution(VECTOR3** pData, VECTOR3** pStd, int nodeNum, int timeSteps, int min_t,
		   int max_t)
{
	assert((pData != NULL) && (nodeNum > 0) && (timeSteps > 0));

	m_nNodeNum = nodeNum;
	m_nTimeSteps = timeSteps;
	m_MinT = min_t; m_MaxT = max_t;
	m_pDataArray = pData;
    // uncertainty
    m_pStdArray = pStd;
    m_bNormalize = false;
}


GaussianSolution::~GaussianSolution()
{
	int iFor;

	if(m_pDataArray != NULL)
	{
		for(iFor = 0; iFor < m_nTimeSteps; iFor++)
		{
			delete[] m_pDataArray[iFor];
		}
		delete[] m_pDataArray;
	}
    // uncertainty
    if (m_pStdArray != NULL)
    {
        for (iFor = 0; iFor < m_nTimeSteps; iFor++)
        {
            delete[] m_pStdArray[iFor];
        }
        delete[] m_pStdArray;
    }
}

void GaussianSolution::Reset()
{
	m_pDataArray = NULL;
	m_nNodeNum = 0;
	m_nTimeSteps = 1;
	m_MinT = 0; m_MaxT = 0;
    // uncertainty
    m_pStdArray = NULL;
}

//////////////////////////////////////////////////////////////////////////
// change vector data on-the-fly
//////////////////////////////////////////////////////////////////////////
void GaussianSolution::SetValue(int t, VECTOR3* pData, int nodeNum)
{
    printf("To prevent memory leak, it is suggested not to set values after solution is constructed.\n");
    Solution::SetValue(t, pData, nodeNum);
}


void GaussianSolution::sample_vector(int id, int t, VECTOR3 &nodeData)
{
    // ensure input is correct
    assert( id >= 0 && id < this->m_nNodeNum && t >= 0 && t < this->m_nTimeSteps  );

    nodeData.Set( box_muller( m_pDataArray[t][id][0], m_pStdArray[t][id][0]),
                  box_muller( m_pDataArray[t][id][1], m_pStdArray[t][id][1]),
                  box_muller( m_pDataArray[t][id][2], m_pStdArray[t][id][2])
                );
}

//////////////////////////////////////////////////////////////////////////
// get value of node id at time step t
// input
//		id:			node Id
//		t:			time step in check
// output
//		nodeData:	vector value at this node
// return
//		1:			operation successful
//		-1:			invalid id
//////////////////////////////////////////////////////////////////////////
int GaussianSolution::GetValue(int id, float t, VECTOR3& nodeData)
{
    float adjusted_t = t - m_MinT;
	if((id < 0) || (id >= m_nNodeNum) || (adjusted_t < 0.0) || (adjusted_t > (float)(m_nTimeSteps-1)))
		return -1;

    if (!isTimeVarying())
        sample_vector( id, (int)adjusted_t, nodeData ) ;
	else
	{
        int lowT = (int)floor(adjusted_t);
        int highT = lowT + 1;
        float ratio = adjusted_t - floorf(adjusted_t);
        // clamp
		if(lowT >= (m_nTimeSteps-1))
		{
			highT = lowT;
			ratio = 0.0;
		}
        VECTOR3 nodeDataLowT, nodeDataHighT;
        sample_vector( id, lowT, nodeDataLowT );
        sample_vector( id, highT, nodeDataHighT );
        nodeData.Set(Lerp(nodeDataLowT[0], nodeDataHighT[0], ratio),
                     Lerp(nodeDataLowT[1], nodeDataHighT[1], ratio),
                     Lerp(nodeDataLowT[2], nodeDataHighT[2], ratio));
	}

    if (m_bNormalize) {
        nodeData.Normalize();
    }

	return 1;
}

//////////////////////////////////////////////////////////////////////////
// to normalize the vector field
// input
// bLocal: whether to normalize in each timestep or through all timesteps.
//		   if locally, then divide its magnitude; if globally, then divide
//		   by the maximal magnitude through the whole field
//////////////////////////////////////////////////////////////////////////
void GaussianSolution::Normalize(bool bLocal)
{
    m_bNormalize = true;
}

void GaussianSolution::Scale(float scaleF)
{
	int iFor, jFor;
	float u, v, w;

	for(iFor = 0; iFor < m_nTimeSteps; iFor++)
	  {
	    for(jFor = 0; jFor < m_nNodeNum; jFor++)
          {
            m_pDataArray[iFor][jFor] = m_pDataArray[iFor][jFor] * scaleF;
            m_pStdArray[iFor][jFor] = m_pStdArray[iFor][jFor] * scaleF;
	      }
	  }
}

void GaussianSolution::Translate(VECTOR3& translate) {
	int iFor, jFor;
	float u, v, w;

	for (iFor = 0; iFor < m_nTimeSteps; iFor++) {
        for (jFor = 0; jFor < m_nNodeNum; jFor++) {
            m_pDataArray[iFor][jFor] = m_pDataArray[iFor][jFor] + translate;
		}
	}
}
