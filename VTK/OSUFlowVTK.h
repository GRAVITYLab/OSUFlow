#ifndef OSUFLOW_VTK
#define OSUFLOW_VTK

#include <assert.h>

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>

#include <Field.h>
#include "OSUFlow.h"

class vtkOSUFlowVectorField : public CVectorField
{
protected:
	vtkSmartPointer<vtkDataSet> sDataset;
	vtkInterpolatedVelocityField* interpolator;

public:
	vtkOSUFlowVectorField(vtkSmartPointer<vtkDataSet> sDataset_) {

		this->Reset();

		sDataset = sDataset_;

		interpolator = vtkInterpolatedVelocityField::New();
		interpolator->AddDataSet(sDataset.GetPointer());
	}
	virtual ~vtkOSUFlowVectorField () {
		interpolator->Delete();
	}
	virtual int lerp_phys_coord(int cellId, CellTopoType eCellTopoType, float* coeff, VECTOR3& pos) {
		printf("lerp_phys_coord Not implemented\n");
		assert(false);
		return -1;
	}
	virtual int at_cell(int cellId, CellTopoType eCellTopoType, const float t, vector<VECTOR3>& vNodeData) {
		printf("at_cell Not implemented\n");
		assert(false);
		return -1;
	}
	virtual int at_slice(int slice, SliceType eSliceType, const float t, vector<VECTOR3>&vSliceData) {
		printf("at_slice Not implemented\n");
		assert(false);
		return -1;
	}
	virtual int at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue) {
		printf("at_vert Not implemented\n");
		assert(false);
		return -1;
	}
	// get vector
	virtual int at_phys(VECTOR3 pos, float t, VECTOR3& vecData) {
		double  coords[4];
		coords[0] = pos[0];
		coords[1] = pos[1];
		coords[2] = pos[2];
		coords[3] = t;
		double vel[3];
        if ( !interpolator->FunctionValues(coords, vel) )
        	return -1;
        vecData[0] = vel[0];
        vecData[1] = vel[1];
        vecData[2] = vel[2];
        return 1;
	}
	// get vector
	virtual int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, VECTOR3& nodeData) {

		pInfo.Set(pos, pInfo.interpolant, fromCell, -1);
		interpolator->SetLastCellId(fromCell);

		double  coords[4];
		coords[0] = pos[0];
		coords[1] = pos[1];
		coords[2] = pos[2];
		coords[3] = t;
		double vel[3];
		if ( !interpolator->FunctionValues(coords, vel) )
			return -1;
		nodeData[0] = vel[0];
		nodeData[1] = vel[1];
		nodeData[2] = vel[2];
		pInfo.inCell = interpolator->GetLastCellId();
		return 1;
	}
	virtual int at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue) {
		printf("Not implemented\n");
		assert(false);
		return -1;
	}
	// get cell volume
	virtual float volume_of_cell(int cellId) {
		double bounds[6];
		sDataset->GetCellBounds(cellId, bounds);
		float v = (bounds[1]-bounds[0])*(bounds[3]-bounds[2])*(bounds[5]-bounds[4]);
		printf("v=%f\n", v);
		return v;
	}
	virtual void NormalizeField(bool bLocal) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void ScaleField(float scale) {
		printf("Not implemented\n");
		assert(false);
	}

	virtual bool IsNormalized(void) {
		printf("Not implemented\n");
		assert(false);
		return false;
	}
	virtual void getDimension(int& xdim, int& ydim, int& zdim) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual CellType GetCellType(void) {
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
	virtual void GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t) {
		assert(false);
		printf("Not implemented\n");
	}
	virtual void GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GetInflowSlice(vector<VECTOR3>& inflowVerts, const float t, const int slice, const SliceType eSliceType)
	{
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GetOutflowSlice(vector<VECTOR3>& outflowVerts, const float t, const int slice, const SliceType eSliceType)
	{
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts, const float t, const int slice, const SliceType eSliceType) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void SetBoundary(VECTOR3 minB, VECTOR3 maxB) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void at_curl(int, VECTOR3&, VECTOR3&) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize)
	{
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GenerateVortField(int t, bool bToNormalize, VECTOR3* pVort) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual void GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap) {
		printf("Not implemented\n");
		assert(false);
	}
	virtual bool IsInRealBoundaries(PointInfo& p)
	{
		double *bounds = sDataset->GetBounds();
		return (p.phyCoord[0] >= bounds[0] && p.phyCoord[0] <= bounds[1] &&
				p.phyCoord[1] >= bounds[2] && p.phyCoord[1] <= bounds[3] &&
				p.phyCoord[2] >= bounds[4] && p.phyCoord[2] <= bounds[5]);
	}
	virtual bool IsInRealBoundaries(PointInfo& p, float time) {
		printf("Not implemented\n");
		assert(false);
		return false;
	}
};


#include <vtkXMLFileReadTester.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkGenericDataObjectReader.h>

class OSUFlowVTK: public OSUFlow
{
public:
	void setData(vtkSmartPointer<vtkDataSet> input)
	{
		double *bounds = input->GetBounds();

		gMin.Set(bounds[0], bounds[2], bounds[4]);
		gMax.Set(bounds[1], bounds[3], bounds[5]);
		lMin = gMin;
		lMax = gMax;
		MinT = 0;
		MaxT = 0;

		CVectorField* field = new vtkOSUFlowVectorField(input);
		flowField = field;
		has_data = true;
	}

#if 0
	vtkObject *loadVTKData(const char *fileName)
	{
	    vtkXMLFileReadTester *vtkXMLFormatFileChecker = vtkXMLFileReadTester::New();

	    vtkXMLFormatFileChecker->SetFileName(fileName);

		if (vtkXMLFormatFileChecker->TestReadFile() > 0) {
			//Logger.getLogger(getClass().getName()).log(Level.INFO, "[{0}] is of type [{1}]",
			//		new Object[] { fileName, vtkXMLFormatFileChecker.GetFileDataType() });

			vtkXMLGenericDataObjectReader *vtkXMLFileReader = vtkXMLGenericDataObjectReader::New();

			vtkXMLFileReader->SetFileName(fileName);
			vtkXMLFileReader->Update();

			if (vtkXMLFileReader->GetOutput()->IsA("vtkMultiBlockDataSet")) {
				return vtkXMLFileReader->GetOutput();
			} else {
				vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

				makeMultiblock->SetInputConnection(vtkXMLFileReader->GetOutputPort());

				return makeMultiblock;
			}
		} else // legacy format
		{
			vtkGenericDataObjectReader *legacyVTKFileReader = vtkGenericDataObjectReader::New();

			legacyVTKFileReader->SetFileName(fileName);
			legacyVTKFileReader->Update();

			vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

			makeMultiblock->SetInput(legacyVTKFileReader->GetOutput());

			return makeMultiblock;
		}
	}
#endif
};

#endif //OSUFLOW_VTK
