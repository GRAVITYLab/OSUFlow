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
	virtual ~CVectorField();

	virtual int lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos);
	virtual int at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData);
	virtual int at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData);
	virtual int at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue);
	virtual int at_phys(VECTOR3 pos, float t, VECTOR3& vecData);
	virtual int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData);
	virtual int at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue);
	virtual float volume_of_cell(int cellId);
	virtual void NormalizeField(bool bLocal);
	virtual void ScaleField(float scale);
	virtual void TranslateField(VECTOR3& translate);

	// ADD-BY-LEETEN 02/02/2012-BEGIN
	virtual void Scan
	  (
	   void (*func)(int iLocalT, int iNode, VECTOR3 *pv3)
	   );
	// ADD-BY-LEETEN 02/02/2012-END

	virtual bool IsNormalized(void);
	virtual void getDimension(int& xdim, int& ydim, int& zdim);
	virtual CellType GetCellType(void) { return m_pGrid->GetCellType(); }
	virtual int GetTimeSteps(void) {return m_nTimeSteps;}
	virtual int GetMinTimeStep(void) {return m_MinT;}
	virtual int GetMaxTimeStep(void) {return m_MaxT;}
	virtual void GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t);
	virtual void GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t);
	virtual void GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t);
	virtual void GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType);
	virtual void GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType);
	virtual void GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType);
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) { m_pGrid->Boundary(minB, maxB); };
	virtual void SetBoundary(VECTOR3 minB, VECTOR3 maxB) {m_pGrid->SetBoundary(minB, maxB);
	}
	virtual void at_curl(int, VECTOR3&, VECTOR3&);
	virtual void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize)
	{ m_pGrid->BoundaryIntersection(intersectP, startP, endP, stepSize, oldStepSize); }
	virtual void GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort);
	virtual void GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap);
	virtual bool IsInRealBoundaries(PointInfo& p);
	virtual bool IsInRealBoundaries(PointInfo& p, float time);

	// feature computation - static
	MATRIX3 Jacobian(const VECTOR3& pos);
	MATRIX3 UnitJacobian(const VECTOR3& pos);

	void GenerateVortexMetrics(const VECTOR3& pos, float& lambda2, float& q, float& delta, float& gamma2);
	void GenerateVortexMetricsLine(VECTOR3* const fieldline, const int num, float* lambda2, float* q, float* delta, float* gamma2);
	void Curvature(VECTOR3* const fieldline, const int num, float* curvature);

protected:
	// reset
	virtual void Reset(void);
	// field functions
	virtual bool isTimeVarying(void);
	// curl
	virtual void curl(float, float, float, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&, VECTOR3&);
};

// file polynomials.c
extern double cube_root(double x);
extern int solve_cubic(float, float, float, float, float*, float*, float*);
extern int solve_quadratic(float, float, float, float*, float*);
extern int solve_linear(float, float, float*);

// file eigenvals.c
extern int compute_eigenvalues(float m[3][3], float eigenvalues[3]);
extern void compute_real_eigenvectors(float m[3][3], float vals[3], float vecs[3][3]);
extern void compute_complex_eigenvectors(float m[3][3], float vals[3], float vecs[3][3]);

#endif
