/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Vector Field: 3D Static or Time-Varying
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _VECTOR_FIELD_H_
#define _VECTOR_FIELD_H_

#include "header.h"
#include "VectorMatrix.h"
#include "Grid.h"
#include "Solution.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

//////////////////////////////////////////////////////////////////////////
// vector field class
//////////////////////////////////////////////////////////////////////////

class CVectorField
{
private:
	Grid* m_pGrid;						// grid
	Solution* m_pSolution;				// vector data
	int m_nTimeSteps;
	bool m_bIsNormalized;				// whether the solution is normalized or not
	int m_MinT, m_MaxT; // the min and max time step of the data field 

public:
	// constructor and destructor
	CVectorField();
	CVectorField(Grid* pGrid, Solution* pSolution, int timesteps, int min_t=0);
	~CVectorField();

	int lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos);
	int at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData);
	int at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData);
	int at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue);
	int at_phys(VECTOR3 pos, float t, VECTOR3& vecData);
	int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData);
	int at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue);
	float volume_of_cell(int cellId);
	void NormalizeField(bool bLocal);
	void ScaleField(float scale); 


	bool IsNormalized(void);
	void getDimension(int& xdim, int& ydim, int& zdim);
	CellType GetCellType(void) { return m_pGrid->GetCellType(); }
	int GetTimeSteps(void) {return m_nTimeSteps;}
	int GetMinTimeStep(void) {return m_MinT;}
	int GetMaxTimeStep(void) {return m_MaxT;}
	void GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t);
	void GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t);
	void GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t);
	void GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType);
	void GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType);
	void GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType);
	void Boundary(VECTOR3& minB, VECTOR3& maxB) { m_pGrid->Boundary(minB, maxB); };
	void SetBoundary(VECTOR3 minB, VECTOR3 maxB) {m_pGrid->SetBoundary(minB, maxB); 
	}
	void at_curl(int, VECTOR3&, VECTOR3&);
	void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize) 
	{ m_pGrid->BoundaryIntersection(intersectP, startP, endP, stepSize, oldStepSize); }
	void GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort);
	void GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap);

protected:
	// reset
	void Reset(void);
	// field functions
	bool isTimeVarying(void);
	// curl
	void curl(float, float, float, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&);
};

#endif
