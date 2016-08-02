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
#include <vtkMultiBlockDataSet.h>
#include <Field.h>
#include <vtkPointData.h>
#include <vtkNew.h>
#include <vtkMath.h>
#include "MultiBlockVectorFieldVTK.h"

// OPENMP currently not used, because core OSUFlow has to be modified too
#ifdef WITH_OPENMP
#include <omp.h>
#endif

MultiBlockVectorFieldVTK::MultiBlockVectorFieldVTK(vtkMultiBlockDataSet *sDataset_)
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
 MultiBlockVectorFieldVTK::~MultiBlockVectorFieldVTK () {
	for (size_t i=0; i < this->interpolatorAry.size(); i++)
		this->interpolatorAry[i]->Delete();
}
void MultiBlockVectorFieldVTK::push_interpolatorAry(vtkMultiBlockDataSet *mbData)
{
	vtkAbstractInterpolatedVelocityField *interpolator = NULL;

    int b;
    vtkInterpolatedVelocityField *func = vtkInterpolatedVelocityField::New();
    bool first = true;
    for (b=0; b < mbData->GetNumberOfBlocks(); b++)
    {
        vtkDataSet *dataset = vtkDataSet::SafeDownCast( mbData->GetBlock(b) );
        if (!dataset) {
            fprintf(stderr, "Hierarchical multiblock dataset not supported yet\n");
            return;
        }
        double *bounds = dataset->GetBounds();
        printf("bounds: %f %f %f %f %f %f\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
        if (first) {
            gbounds[0] = bounds[0];
            gbounds[1] = bounds[1];
            gbounds[2] = bounds[2];
            gbounds[3] = bounds[3];
            gbounds[4] = bounds[4];
            gbounds[5] = bounds[5];
        } else {
			gbounds[0] = min(bounds[0], gbounds[0]);
			gbounds[1] = max(bounds[1], gbounds[1]);
			gbounds[2] = min(bounds[2], gbounds[2]);
			gbounds[3] = max(bounds[3], gbounds[3]);
			gbounds[4] = min(bounds[4], gbounds[4]);
			gbounds[5] = max(bounds[5], gbounds[5]);
        }
        first = false;

        // vector
        func->AddDataSet(dataset);
    }
    interpolator = func;

	this->interpolatorAry.push_back(interpolator);

}

//////////////////////////////////
 int MultiBlockVectorFieldVTK::lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos) {
	printf("lerp_phys_coord Not implemented\n");
	assert(false);
	return -1;
}
 int MultiBlockVectorFieldVTK::at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData) {
	printf("at_cell Not implemented\n");
	assert(false);
	return -1;
}
 int MultiBlockVectorFieldVTK::at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData) {
	printf("at_slice Not implemented\n");
	assert(false);
	return -1;
}
 int MultiBlockVectorFieldVTK::at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue) {
	printf("at_vert Not implemented\n");
	assert(false);
	return -1;
}
 int MultiBlockVectorFieldVTK::phys_coord(const int i, const int j, const int k, VECTOR3 &pos) {
     printf("at_vert Not implemented\n");
     assert(false);
     return -1;
 }

// get vector
 int MultiBlockVectorFieldVTK::at_phys(const VECTOR3 &pos, float t, VECTOR3& vecData) {
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
 int MultiBlockVectorFieldVTK::at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData) {
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
// get cell volume
 float MultiBlockVectorFieldVTK::volume_of_cell(int cellId) {
     return 1;
    //double bounds[6];
    //sDataset->GetCellBounds(cellId, bounds);
    //float v = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]);
    ////printf("v=%f\n", v);
    //return v;
}
 void MultiBlockVectorFieldVTK::NormalizeField(bool bLocal) {
  // GetVectors() assumes the active vectors are set by the main program.
  // If not, add this command: vtkdata->GetPointData()->SetActiveVectors(array);
  int i;
  for (i=0; i<this->sDataset->GetNumberOfBlocks(); i++)
  {
    vtkDataArray *ary = vtkDataSet::SafeDownCast( this->sDataset->GetBlock(i) )->GetPointData()->GetVectors();
    assert(ary);
    vtkNew<vtkMath> math;
    int id;
    int n = ary->GetNumberOfTuples();
    for (id=0; id<n; id++)
    {
      double vec[10]; // in case for many components
      ary->GetTuple(id, vec);
      math->Normalize(vec); // always normalize first 3 components
      ary->SetTuple(id, vec);
    }
  }
}
 void MultiBlockVectorFieldVTK::ScaleField(float scale) {
	this->scaleFactor = scale;
}

 bool MultiBlockVectorFieldVTK::IsNormalized(void) {
	printf("Not implemented\n");
	assert(false);
	return false;
}
 void MultiBlockVectorFieldVTK::getDimension(int& xdim, int& ydim, int& zdim) {
	printf("Not implemented\n");
	assert(false);
}
 CellType MultiBlockVectorFieldVTK::GetCellType(void) {
    int type = vtkDataSet::SafeDownCast( sDataset->GetBlock(0))->GetCellType(0);
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
 void MultiBlockVectorFieldVTK::GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t) {
	assert(false);
	printf("Not implemented\n");
}
 void MultiBlockVectorFieldVTK::GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType)
{
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType)
{
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::Boundary(VECTOR3& minB, VECTOR3& maxB) {
    minB.Set(gbounds[0], gbounds[2], gbounds[4]);
    maxB.Set(gbounds[1], gbounds[3], gbounds[5]);
	//printf("Bound: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
}
 void MultiBlockVectorFieldVTK::SetBoundary(VECTOR3 minB, VECTOR3 maxB) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::at_curl(int, VECTOR3&, VECTOR3&) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize)
{
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort) {
	printf("Not implemented\n");
	assert(false);
}
 void MultiBlockVectorFieldVTK::GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap) {
	printf("Not implemented\n");
	assert(false);
}
 bool MultiBlockVectorFieldVTK::IsInRealBoundaries(PointInfo& p)
{
	VECTOR3 minB, maxB;
	this->Boundary(minB, maxB);
	return (p.phyCoord[0] >= minB[0] && p.phyCoord[0] <= maxB[0] &&
		 p.phyCoord[1] >= minB[1] && p.phyCoord[1] <= maxB[1] &&
		 p.phyCoord[2] >= minB[2] && p.phyCoord[2] <= maxB[2]);
}
 bool MultiBlockVectorFieldVTK::IsInRealBoundaries(PointInfo& p, float time) {
	printf("Not implemented\n");
	assert(false);
	return false;
}

int MultiBlockVectorFieldVTK::getThreadID()
{
#ifdef WITH_OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}
