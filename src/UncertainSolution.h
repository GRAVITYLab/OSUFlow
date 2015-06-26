/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li
//                 The Ohio State University
//                 Date:		06/2005
//                 Vector Field Data
//
///////////////////////////////////////////////////////////////////////////////



#ifndef _UNCERTAIN_SOLUTION_H_
#define _UNCERTAIN_SOLUTION_H_

#include "Solution.h"

//Silence an annoying and unnecessary compiler warning
//#pragma warning(disable : 4251 4100 4244)

template <class T>
T box_muller(T m, T s);	/* normal random variate generator */


class GaussianSolution: public Solution
{
protected:
    VECTOR3** m_pStdArray;				// data value
    bool m_bNormalize ;                   // normalize the vector after sampling

    // what we have in Solution:
//	VECTOR3** m_pDataArray;				// data value
//	VECTOR3* m_pMinValue;				// value with min magnitude for each time step
//	VECTOR3* m_pMaxValue;				// value with max magnitude for each time step
//	int m_nNodeNum;						// how many nodes each time step
//	int m_nTimeSteps;					// how many time steps
//	float m_fMinMag;					// minimal magnitude
//	float m_fMaxMag;					// maximum magnitude
//	int m_MinT, m_MaxT;

    void sample_vector( int id, int t, VECTOR3& nodeData );

public:
	// constructor
    GaussianSolution(VECTOR3** pData, VECTOR3** pStd, int nodeNum, int timeSteps);
    GaussianSolution(VECTOR3** pData, VECTOR3** pStd, int nodeNum, int timeSteps, int min_t, int max_t);
    ~GaussianSolution();

    virtual void Reset();

    // solution functions
    virtual void SetValue(int t, VECTOR3* pData, int nodeNum);
    virtual int GetValue(int id, const float t, VECTOR3& nodeData);
    virtual void Normalize(bool bLocal);
    virtual void Scale(float scaleF);
    virtual void Translate(VECTOR3& translate);

};

#endif
