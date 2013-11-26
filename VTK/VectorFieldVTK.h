#ifndef VECTOR_FIELD_VTK
#define VECTOR_FIELD_VTK

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>
#include <Field.h>

class VectorFieldVTK : public CVectorField
{
protected:
	vtkSmartPointer<vtkDataSet> sDataset;
	vector<vtkAbstractInterpolatedVelocityField*> interpolatorAry;
	int getThreadID();
	void push_interpolatorAry(vtkDataSet *);
	float scaleFactor;
public:
	VectorFieldVTK(vtkDataSet *sDataset_) ;
	virtual ~VectorFieldVTK () ;
	virtual int lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos) ;
	virtual int at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData) ;
	virtual int at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData) ;
	virtual int at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue) ;
	// get vector
	virtual int at_phys(const VECTOR3 &pos, float t, VECTOR3& vecData) ;
	// get vector
	virtual int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData) ;
	virtual int at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue) ;
	// get cell volume
	virtual float volume_of_cell(int cellId) ;
	virtual void NormalizeField(bool bLocal) ;
	virtual void ScaleField(float scale) ;

	virtual bool IsNormalized(void) ;
	virtual void getDimension(int& xdim, int& ydim, int& zdim) ;
	virtual CellType GetCellType(void) ;
	virtual void GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t) ;
	virtual void GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t) ;
	virtual void GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t) ;
	virtual void GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType);
	virtual void GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType);
	virtual void GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType) ;
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) ;
	virtual void SetBoundary(VECTOR3 minB, VECTOR3 maxB) ;
	virtual void at_curl(int, VECTOR3&, VECTOR3&) ;
	virtual void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize);
	virtual void GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort) ;
	virtual void GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap) ;
	virtual bool IsInRealBoundaries(PointInfo& p);

	virtual bool IsInRealBoundaries(PointInfo& p, float time) ;
};

#endif //VECTOR_FIELD_VTK
