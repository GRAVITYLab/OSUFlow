#include <assert.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkOverlappingAMR.h>
#include <vtkAMRInterpolatedVelocityField.h>
#include <vtkCellLocatorInterpolatedVelocityField.h>
#include <vtkImageData.h>
#include <vtkImageInterpolator.h>
#include <Field.h>
#include "VectorFieldVTK.h"

// OPENMP currently not used, because core OSUFlow has to be modified too
#ifdef WITH_OPENMP
#include <omp.h>
#endif

VectorFieldVTK::VectorFieldVTK(vtkDataSet *sDataset_)
: scaleFactor(1)
{
	this->sDataset = sDataset_ ;
	this->Reset();

#ifdef WITH_OPENMP
	for (size_t i=0; i<omp_get_max_threads(); i++)
		this->push_interpolatorAry(sDataset_);
#else
	this->push_interpolatorAry(sDataset_);
#endif
}
 VectorFieldVTK::~VectorFieldVTK () {
	for (size_t i=0; i < this->interpolatorAry.size(); i++)
		this->interpolatorAry[i]->Delete();
}
void VectorFieldVTK::push_interpolatorAry(vtkDataSet *data)
{
	vtkOverlappingAMR* amrData = vtkOverlappingAMR::SafeDownCast(data);
	vtkAbstractInterpolatedVelocityField *interpolator = NULL;

	if(amrData)
	{
		vtkAMRInterpolatedVelocityField *func = vtkAMRInterpolatedVelocityField::New();
		func->SetAMRData(amrData);
		interpolator = func;
	}
	else
	{
		vtkUnstructuredGrid *unstructured = vtkUnstructuredGrid::SafeDownCast(data);
		vtkImageData *image = vtkImageData::SafeDownCast(data);
		if (unstructured) {
			printf("Unstrcutred data\n");
			vtkCellLocatorInterpolatedVelocityField *func = vtkCellLocatorInterpolatedVelocityField::New();
			func->AddDataSet(data);
			interpolator = func;
		} else if (image) { // regular grid
			printf("Error: Image Data : Use CVectorField class\n");

		} else // structured grid
		{
			printf("Structured Data\n");
			vtkInterpolatedVelocityField *func = vtkInterpolatedVelocityField::New();
			func->AddDataSet(data);
			interpolator = func;
		}
	}

	this->interpolatorAry.push_back(interpolator);

}

//////////////////////////////////
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
 int VectorFieldVTK::at_phys(const VECTOR3 &pos, float t, VECTOR3& vecData) {
	double  coords[4];
	coords[0] = pos[0];
	coords[1] = pos[1];
	coords[2] = pos[2];
	coords[3] = t;
	double vel[3];
	bool success = this->interpolatorAry[ this->getThreadID() ]->FunctionValues(coords, vel) ;
	//printf("first querry success = %d, coord = %lf %lf %lf %lf,  vec = %lf %lf %lf\n", success, coords[0], coords[1], coords[2], coords[3], vel[0], vel[1], vel[2]);
	if (!success)
		return -1;
	vecData[0] = vel[0]*scaleFactor;
	vecData[1] = vel[1]*scaleFactor;
	vecData[2] = vel[2]*scaleFactor;
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
	//printf("querry success = %d, coord = %lf %lf %lf %lf,  vec = %lf %lf %lf\n", success, coords[0], coords[1], coords[2], coords[3], vel[0], vel[1], vel[2]);
	if (!success) return -1;

	nodeData[0] = vel[0]*scaleFactor;
	nodeData[1] = vel[1]*scaleFactor;
	nodeData[2] = vel[2]*scaleFactor;
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
	this->scaleFactor = scale;
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
	vtkStructuredGrid *structuredGrid;

	if ((structuredGrid = vtkStructuredGrid::SafeDownCast(sDataset)) != NULL || vtkImageData::SafeDownCast(sDataset)!=NULL) {
		int *bounds = structuredGrid->GetExtent();
		minB.Set(bounds[0], bounds[2], bounds[4]);
		maxB.Set(bounds[1], bounds[3], bounds[5]);
		//printf("Structured Extent: %d %d %d %d %d %d\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	}
	//else {
		double *bounds = sDataset->GetBounds();
		minB.Set(bounds[0], bounds[2], bounds[4]);
		maxB.Set(bounds[1], bounds[3], bounds[5]);
		//printf("Bound: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
	//	}
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
	VECTOR3 minB, maxB;
	this->Boundary(minB, maxB);
	return (p.phyCoord[0] >= minB[0] && p.phyCoord[0] <= maxB[0] &&
		 p.phyCoord[1] >= minB[1] && p.phyCoord[1] <= maxB[1] &&
		 p.phyCoord[2] >= minB[2] && p.phyCoord[2] <= maxB[2]);
}
 bool VectorFieldVTK::IsInRealBoundaries(PointInfo& p, float time) {
	printf("Not implemented\n");
	assert(false);
	return false;
}

int VectorFieldVTK::getThreadID()
{
#ifdef WITH_OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}
