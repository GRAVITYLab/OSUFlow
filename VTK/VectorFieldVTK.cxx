#include <assert.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>


#include <Field.h>
#include "VectorFieldVTK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

VectorFieldVTK::VectorFieldVTK(vtkDataSet *sDataset_)
: sDataset(sDataset_)
{

	this->Reset();

#ifdef _OPENMP
	for (size_t i=0; i<omp_get_max_threads(); i++)
	{
		vtkInterpolatedVelocityField *interpolator = vtkInterpolatedVelocityField::New();
		interpolator->AddDataSet(sDataset);
		this->interpolatorAry.push_back(interpolator);
	}
#else
	vtkInterpolatedVelocityField *interpolator = vtkInterpolatedVelocityField::New();
	interpolator->AddDataSet(sDataset);
	this->interpolatorAry.push_back(interpolator);
#endif
}
 VectorFieldVTK::~VectorFieldVTK () {
	for (size_t i=0; i < this->interpolatorAry.size(); i++)
		this->interpolatorAry[i]->Delete();
}
 int VectorFieldVTK::lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos) {
	printf("lerp_phys_coord Not implemented\n");
	assert(false);
	return -1;
}
 int VectorFieldVTK::at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData) {
	printf("at_cell Not implemented\n");
	assert(false);
	return -1;
}
 int VectorFieldVTK::at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData) {
	printf("at_slice Not implemented\n");
	assert(false);
	return -1;
}
 int VectorFieldVTK::at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue) {
	printf("at_vert Not implemented\n");
	assert(false);
	return -1;
}
// get vector
 int VectorFieldVTK::at_phys(VECTOR3 pos, float t, VECTOR3& vecData) {
	double  coords[4];
	coords[0] = pos[0];
	coords[1] = pos[1];
	coords[2] = pos[2];
	coords[3] = t;
	double vel[3];
	if ( ! this->interpolatorAry[ this->getThreadID() ]->FunctionValues(coords, vel) )
		return -1;
	vecData[0] = vel[0];
	vecData[1] = vel[1];
	vecData[2] = vel[2];
	return 1;
}
// get vector
 int VectorFieldVTK::at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData) {
	bool success ;
	double vel[3];
	pInfo.Set(pos, pInfo.interpolant, fromCell, -1);
	this->interpolatorAry[ this->getThreadID() ]->SetLastCellId(fromCell);

	double  coords[4];
	coords[0] = pos[0];
	coords[1] = pos[1];
	coords[2] = pos[2];
	coords[3] = t;
	success = this->interpolatorAry[ this->getThreadID() ]->FunctionValues(coords, vel);
	if (!success) return -1;

	nodeData[0] = vel[0];
	nodeData[1] = vel[1];
	nodeData[2] = vel[2];
	pInfo.inCell = this->interpolatorAry[ this->getThreadID() ]->GetLastCellId();
	return 1;
}
 int VectorFieldVTK::at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue) {
	printf("Not implemented\n");
	assert(false);
	return -1;
}
// get cell volume
 float VectorFieldVTK::volume_of_cell(int cellId) {
	double bounds[6];
	sDataset->GetCellBounds(cellId, bounds);
	float v = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]);
	//printf("v=%f\n", v);
	return v;
}
 void VectorFieldVTK::NormalizeField(bool bLocal) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::ScaleField(float scale) {
	printf("Not implemented\n");
	assert(false);
}

 bool VectorFieldVTK::IsNormalized(void) {
	printf("Not implemented\n");
	assert(false);
	return false;
}
 void VectorFieldVTK::getDimension(int& xdim, int& ydim, int& zdim) {
	printf("Not implemented\n");
	assert(false);
}
 CellType VectorFieldVTK::GetCellType(void) {
	int type = sDataset->GetCellType(0);
	switch (type) {
	case VTK_TRIANGLE:
		return TRIANGLE;
	case VTK_VOXEL:
		return CUBE;
	case VTK_POLYGON:
		return POLYGON;
	case VTK_TETRA:
		return TETRAHEDRON;
	default:
		printf("OSUFlow: Unsupported cell type: %d\n", type);
		assert(false);
	};
	return CELLTYPE_UNKNOWN;
}
 void VectorFieldVTK::GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t) {
	assert(false);
	printf("Not implemented\n");
}
 void VectorFieldVTK::GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType)
{
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType)
{
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::Boundary(VECTOR3& minB, VECTOR3& maxB) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::SetBoundary(VECTOR3 minB, VECTOR3 maxB) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::at_curl(int, VECTOR3&, VECTOR3&) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize)
{
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort) {
	printf("Not implemented\n");
	assert(false);
}
 void VectorFieldVTK::GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap) {
	printf("Not implemented\n");
	assert(false);
}
 bool VectorFieldVTK::IsInRealBoundaries(PointInfo& p)
{
	double *bounds = sDataset->GetBounds();
	return (p.phyCoord[0] >= bounds[0] && p.phyCoord[0] <= bounds[1] &&
			p.phyCoord[1] >= bounds[2] && p.phyCoord[1] <= bounds[3] &&
			p.phyCoord[2] >= bounds[4] && p.phyCoord[2] <= bounds[5]);
}
 bool VectorFieldVTK::IsInRealBoundaries(PointInfo& p, float time) {
	printf("Not implemented\n");
	assert(false);
	return false;
}

int VectorFieldVTK::getThreadID()
{
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}
