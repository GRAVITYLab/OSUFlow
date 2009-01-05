/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 FieldLine
//
///////////////////////////////////////////////////////////////////////////////

#include "Solution.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
//
// definition of CSolution class
//
//////////////////////////////////////////////////////////////////////////

Solution::Solution()
{
	Reset();
}

//////////////////////////////////////////////////////////////////////////
// input
//		pData:		this is a 2D array for storing static or 
//					time-varying data
//		nodeNum:	number of total nodes
//		timeSteps:	how many time steps, for static data, this is 1
//////////////////////////////////////////////////////////////////////////
// initialzie data is the data is not that large
Solution::Solution(VECTOR3** pData, int nodeNum, int timeSteps)
{
	assert((pData != NULL) && (nodeNum > 0) && (timeSteps > 0));

	m_nNodeNum = nodeNum;
	m_nTimeSteps = timeSteps;	
	m_MinT = 0; m_MaxT = timeSteps-1; 
	m_pDataArray = new VECTOR3*[timeSteps];

	for(int iFor = 0; iFor < timeSteps; iFor++)
	{
		m_pDataArray[iFor] = new VECTOR3[nodeNum];
		assert(m_pDataArray[iFor] != NULL);
		for(int jFor = 0; jFor < nodeNum; jFor++)
			m_pDataArray[iFor][jFor] = pData[iFor][jFor];
	}	
}

Solution::Solution(VECTOR3** pData, int nodeNum, int timeSteps, int min_t, 
		   int max_t)
{
	assert((pData != NULL) && (nodeNum > 0) && (timeSteps > 0));

	m_nNodeNum = nodeNum;
	m_nTimeSteps = timeSteps;	
	m_MinT = min_t; m_MaxT = max_t; 
	m_pDataArray = new VECTOR3*[timeSteps];

	for(int iFor = 0; iFor < timeSteps; iFor++)
	{
		m_pDataArray[iFor] = new VECTOR3[nodeNum];
		assert(m_pDataArray[iFor] != NULL);
		for(int jFor = 0; jFor < nodeNum; jFor++)
			m_pDataArray[iFor][jFor] = pData[iFor][jFor];
	}	
}


Solution::Solution(int nodeNum, int timeSteps)
{
	m_nNodeNum = nodeNum;
	m_nTimeSteps = timeSteps;		
	m_MinT = 0; m_MaxT = timeSteps-1; 	
	m_pDataArray = new VECTOR3*[timeSteps];

	for(int iFor = 0; iFor < timeSteps; iFor++)
		m_pDataArray[iFor] = NULL;
}

Solution::~Solution()
{
	int iFor;

	if(m_pDataArray != NULL)
	{
		for(iFor = 0; iFor < m_nTimeSteps; iFor++)
			delete[] m_pDataArray[iFor];
		delete[] m_pDataArray;
	}
}

void Solution::Reset()
{
	m_pDataArray = NULL;
	m_nNodeNum = 0;
	m_nTimeSteps = 1;			
	m_MinT = 0; m_MaxT = 0; 
}

//////////////////////////////////////////////////////////////////////////
// change vector data on-the-fly
//////////////////////////////////////////////////////////////////////////
void Solution::SetValue(int t, VECTOR3* pData, int nodeNum)
{
  //	if((t >= 0) && (t < m_nTimeSteps))
  	if((t >= m_MinT) && (t <= m_MaxT))
	{
	  t = t-m_MinT; 
	  m_pDataArray[t] = new VECTOR3[nodeNum];
	  assert(m_pDataArray[t] != NULL);
	  for(int jFor = 0; jFor < nodeNum; jFor++)
	    m_pDataArray[t][jFor] = pData[jFor];
	}
}

//////////////////////////////////////////////////////////////////////////
// whether field is time varying
//////////////////////////////////////////////////////////////////////////
bool Solution::isTimeVarying(void)
{
	return (m_nTimeSteps > 1);
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
int Solution::GetValue(int id, float t, VECTOR3& nodeData)
{
  float adjusted_t = t - m_MinT; 
	if((id < 0) || (id >= m_nNodeNum) || (adjusted_t < 0.0) || (adjusted_t > (float)(m_nTimeSteps-1)))
		return -1;

	if(!isTimeVarying())
		nodeData = m_pDataArray[(int)adjusted_t][id];
	else
	{
		int lowT, highT;
		float ratio;
		lowT = (int)floor(adjusted_t);
		ratio = adjusted_t - (float)floor(adjusted_t);
		highT = lowT + 1;
		if(lowT >= (m_nTimeSteps-1))
		{
			highT = lowT;
			ratio = 0.0;
		}
		nodeData.Set(Lerp(m_pDataArray[lowT][id][0], m_pDataArray[highT][id][0], ratio), 
					 Lerp(m_pDataArray[lowT][id][1], m_pDataArray[highT][id][1], ratio),
					 Lerp(m_pDataArray[lowT][id][2], m_pDataArray[highT][id][2], ratio));
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
void Solution::Normalize(bool bLocal)
{
	int iFor, jFor;
	float mag, u, v, w;

	m_fMinMag = FLT_MAX;
	m_fMaxMag = -FLT_MAX;

	if(bLocal)
	{
		for(iFor = 0; iFor < m_nTimeSteps; iFor++)
		{
			for(jFor = 0; jFor < m_nNodeNum; jFor++)
			{
				mag = m_pDataArray[iFor][jFor].GetMag();
				if(mag != 0.0)
				{
					u = m_pDataArray[iFor][jFor][0]/mag;
					v = m_pDataArray[iFor][jFor][1]/mag;
					w = m_pDataArray[iFor][jFor][2]/mag;
					m_pDataArray[iFor][jFor].Set(u, v, w);
				}
				
				if(mag < m_fMinMag)
					m_fMinMag = mag;
				if(mag > m_fMaxMag)
					m_fMaxMag = mag;
			}
		}
	}
	else
	{
		for(iFor = 0; iFor < m_nTimeSteps; iFor++)
		{
			for(jFor = 0; jFor < m_nNodeNum; jFor++)
			{
				mag = m_pDataArray[iFor][jFor].GetMag();
				if(mag < m_fMinMag)
					m_fMinMag = mag;
				if(mag > m_fMaxMag)
					m_fMaxMag = mag;
			}
		}

		for(iFor = 0; iFor < m_nTimeSteps; iFor++)
		{
			for(jFor = 0; jFor < m_nNodeNum; jFor++)
			{
				u = m_pDataArray[iFor][jFor][0]/m_fMaxMag;
				v = m_pDataArray[iFor][jFor][1]/m_fMaxMag;
				w = m_pDataArray[iFor][jFor][2]/m_fMaxMag;
				m_pDataArray[iFor][jFor].Set(u, v, w);
			}
		}
	}
}




void Solution::Scale(float scaleF)
{
	int iFor, jFor;
	float u, v, w;

	for(iFor = 0; iFor < m_nTimeSteps; iFor++)
	  {
	    for(jFor = 0; jFor < m_nNodeNum; jFor++)
	      {
		u = m_pDataArray[iFor][jFor][0]*scaleF; 
		v = m_pDataArray[iFor][jFor][1]*scaleF; 
		w = m_pDataArray[iFor][jFor][2]*scaleF; 
		m_pDataArray[iFor][jFor].Set(u, v, w);
	      }
	  }
}



// compute the min and max value with minimal and maximal magnitude
void Solution::ComputeMinMaxValue(void)
{
	int indexMin, indexMax, iFor, jFor;
	float minMag, maxMag, mag;
	
	m_pMinValue = new VECTOR3[m_nTimeSteps];
	m_pMaxValue = new VECTOR3[m_nTimeSteps];

	for(iFor = 0; iFor < m_nTimeSteps; iFor++)
	{
		minMag = FLT_MAX;
		maxMag = -FLT_MAX;
		for(jFor = 0; jFor < m_nNodeNum; jFor++)
		{
			mag = m_pDataArray[iFor][jFor].GetMag();
			if(mag < minMag)
			{
				minMag = mag;
				indexMin = jFor;
			}
			if(mag > maxMag)
			{
				maxMag = mag;
				indexMax = jFor;
			}
		}

		m_pMinValue[iFor] = m_pDataArray[iFor][indexMin];
		m_pMaxValue[iFor] = m_pDataArray[iFor][indexMax];
	}
}

// get the min and max value for all time steps
int Solution::GetMinMaxValueAll(VECTOR3& minVal, VECTOR3& maxVal)
{
	minVal = m_pMinValue[0];
	maxVal = m_pMaxValue[0];

	for(int tFor = 1; tFor < m_nTimeSteps; tFor++)
	{
		if(minVal.GetMag() > m_pMinValue[tFor].GetMag())
			minVal = m_pMinValue[tFor];

		if(maxVal.GetMag() < m_pMaxValue[tFor].GetMag())
			maxVal = m_pMaxValue[tFor];
	}

	return 1;
}

// get the min and max value for timestep t
int Solution::GetMinMaxValue(int t, VECTOR3& minVal, VECTOR3& maxVal)
{
        t = t - m_MinT; 
	if((t >= 0) && (t < m_nTimeSteps))
	{
		minVal = m_pMinValue[t];
		maxVal = m_pMaxValue[t];
		return 1;
	}
	else
		return -1;
}

