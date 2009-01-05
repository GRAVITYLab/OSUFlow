/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Vector Field Data
//
///////////////////////////////////////////////////////////////////////////////



#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "header.h"
#include "VectorMatrix.h"
#include "Interpolator.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

class Solution
{
private:
	VECTOR3** m_pDataArray;				// data value
	VECTOR3* m_pMinValue;				// value with min magnitude for each time step
	VECTOR3* m_pMaxValue;				// value with max magnitude for each time step
	int m_nNodeNum;						// how many nodes each time step
	int m_nTimeSteps;					// how many time steps
	float m_fMinMag;					// minimal magnitude
	float m_fMaxMag;					// maximum magnitude
	int m_MinT, m_MaxT; 

public:
	// constructor
	Solution();
	Solution(VECTOR3** pData, int nodeNum, int timeSteps);
	Solution(VECTOR3** pData, int nodeNum, int timeSteps, int min_t, int max_t);
	Solution(int nodeNum, int timeSteps);
	~Solution();

	void Reset();

	// solution functions
	void SetMinMaxTime(int min_t, int max_t) {m_MinT = min_t; m_MaxT = max_t;}
	void SetValue(int t, VECTOR3* pData, int nodeNum);
	int GetMinMaxValueAll(VECTOR3& minVal, VECTOR3& maxVal);
	int GetMinMaxValue(int t, VECTOR3& minVal, VECTOR3& maxVal); 
	void ComputeMinMaxValue(void);
	bool isTimeVarying(void);
	int GetValue(int id, const float t, VECTOR3& nodeData);
	void Normalize(bool bLocal);
	void Scale(float scaleF); 
};

#endif
