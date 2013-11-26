// code reference: vtkStreamLine.cxx

#include "vtkOSUFlow.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkInformationVector.h"
#include "vtkNew.h"
#include "vtkImageData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

#include "VectorFieldVTK.h"

// Defines static vtkOSUFlow::New() here
vtkStandardNewMacro(vtkOSUFlow);

vtkOSUFlow::vtkOSUFlow()
: IntegratorOrder((int)RK45)
, MaximumIntegrationStep(0.2)
, MinimumIntegrationStep(0.001)
, MaximumError(1e-6)
, MaximumNumberOfSteps(1000)
, scale(1.0)
{
	osuflow = new OSUFlow();
}

vtkOSUFlow::~vtkOSUFlow()
{
	delete osuflow;
}

// Description
// Tells VTK pipeline that both input ports are optional (Data can be assigned to OSUFlow in advance)
// The number of ports (2) are assigned in the super class vtkStreamer
int vtkOSUFlow::FillInputPortInformation(int port, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}


int vtkOSUFlow::RequestData(
	vtkInformation *,
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	int i;
	//vtkFloatArray *newVectors;
	//vtkFloatArray *newScalars=NULL;
	vtkIdType ptId;
	int j;
	//vtkPolyLine* lineNormalGenerator = NULL;
	//vtkFloatArray* normals = NULL;
	//vtkFloatArray* rotation = 0;

	vtkSmartPointer<vtkInformation> sourceInfo = inputVector[0]->GetInformationObject(0); // data
	vtkSmartPointer<vtkInformation> inInfo = inputVector[1]->GetInformationObject(0);  // seeds
	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);	// trace

	vtkSmartPointer<vtkDataSet> source;
	if (sourceInfo.GetPointer()) {
		source = vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkSmartPointer<vtkDataSet> input;
	if (inInfo.GetPointer()) {
		input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkSmartPointer<vtkPolyData> output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	//
	// OSUFlow
	//
	// set data
	if (source) {
		vtkImageData * image = vtkImageData::SafeDownCast(source);
		if (image!=NULL) {
			// use OSUFlow CVectorField
			printf("Image data\n");
			CVectorField* field;
			Solution* pSolution;
			RegularCartesianGrid* pRegularCGrid;
			VECTOR3 min_b, max_b;
			VECTOR4 realMin_b, realMax_b;

			VECTOR3 *pVector = (VECTOR3*)image->GetScalarPointer();
			VECTOR3 **ppVector = new VECTOR3*[1];
			ppVector[0] = pVector;

			int *dims  = image->GetDimensions();
			int totalNum = dims[0]*dims[1]*dims[2];
			pSolution = new Solution(ppVector, totalNum, 1);
			pRegularCGrid = new RegularCartesianGrid(dims[0], dims[1], dims[2]);

			double *bounds = image->GetBounds();
			min_b[0] = bounds[0]; max_b[0] = bounds[1];
			min_b[1] = bounds[2]; max_b[1] = bounds[3];
			min_b[2] = bounds[4]; max_b[2] = bounds[5];
			// Not exactly. TODO:
			realMin_b[0] = min_b[0]; realMin_b[1] = min_b[1];
			realMin_b[2] = min_b[2]; realMin_b[3] = 0;
			realMax_b[0] = max_b[0]; realMax_b[1] = max_b[1];
			realMax_b[2] = max_b[2]; realMax_b[3] = 0;

			pRegularCGrid->SetBoundary(min_b, max_b);
			pRegularCGrid->SetRealBoundary(realMin_b, realMax_b);

			assert(pSolution != NULL && pRegularCGrid != NULL);

			field = new CVectorField(pRegularCGrid, pSolution, 1);
			osuflow->SetFlowField( field );

		} else {
			CVectorField *field = new VectorFieldVTK( source );
			osuflow->SetFlowField( field );
		}
	} else if (! osuflow->HasData() ) {
		printf("vtkOSUFlow: no data\n");
		return 0;
	}

	// assign seeds
	VECTOR3 *pSeed;
	if (input) {
		int num_seeds = input->GetNumberOfPoints();
		printf("num_seed=%d\n", num_seeds);
		if (num_seeds) {
			pSeed = new VECTOR3[num_seeds];
			for (i=0; i < num_seeds; i++) {
				double pos[3];
				input->GetPoint(i,pos);
				pSeed[i] = VECTOR3(pos[0], pos[1], pos[2]);
			}
			osuflow->SetSeedPoints(pSeed, num_seeds);
		}
	} else {
		// use pre-loaded seeds if exist
		int num_seeds;
		osuflow->GetSeeds(num_seeds);
		if (num_seeds==0) {
			printf("vtkOSUFlow: No seeds\n");
			return 0;
		}
	}

	// integrate
	TRACE_DIR dir;
	switch (this->IntegrationDirection)
	{
	case VTK_INTEGRATE_FORWARD: dir=FORWARD_DIR; break;
	case VTK_INTEGRATE_BACKWARD: dir=BACKWARD_DIR; break;
	default: dir = BACKWARD_AND_FORWARD; break;
	};
	list<vtListSeedTrace*> list;
	osuflow->SetIntegrationOrder((INTEG_ORD)this->IntegratorOrder);
	osuflow->SetMaxError(this->MaximumError);
	osuflow->SetIntegrationParams(this->IntegrationStepLength, this->MinimumIntegrationStep, this->MaximumIntegrationStep);

	osuflow->GenStreamLines(list , dir, this->MaximumNumberOfSteps, 0); // default: RK45

	delete[] pSeed;

	//
	// convert traces from list-of-list to vtkPolyData
	//
	vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkFloatArray> fieldData = vtkSmartPointer<vtkFloatArray>::New();

#if 1
	fieldData->SetNumberOfComponents(1);
	fieldData->Allocate(1000);
	fieldData->SetName("Curvature");

	std::list<vtListSeedTrace*>::iterator pIter;
	pIter = list.begin();
	for (; pIter!=list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter;
		std::list<VECTOR3*>::iterator pnIter;

		for (pnIter = trace->begin(); pnIter!=trace->end(); pnIter++) {
			VECTOR3 &p = **pnIter;
			//printf(" %f %f %f ", p[0], p[1], p[2]);

			//vtk
			ptId = newPts->InsertNextPoint((float *)&p[0]);
			pts->InsertNextId(ptId);

			// field data
			float q;
			VECTOR3 vec(0,0,0);
			//osuflow->GetFlowField()->Curvature(&p, 1, &q);
			osuflow->GetFlowField()->at_phys(p, 0, vec);  q = vec.GetMag();
			fieldData->InsertTuple(ptId, &q);

			// clear up
			delete *pnIter;
		}
		newLines->InsertNextCell(pts);
		pts->Reset();
		//printf("\n");

		// clear up
		delete trace;
	}
#else
	// test
	for (int j=0; j<3; j++)
	{
		for (int i=0; i<100; i++)
		{
			float p[3];
			p[0]=p[1]=p[2]=i;
			p[0]+=j*2;
			ptId = newPts->InsertNextPoint(p);
			//line->GetPointIds()->InsertNextId(ptId);
			pts->InsertNextId(ptId);
		}
		newLines->InsertNextCell(pts);
		pts->Reset();
	}
	output->GetPointData()->SetNormals()
#endif

	//
	// assign lines to output
	//
	//if (newPts->GetNumberOfPoints() > 0)
	{
		output->SetPoints(newPts);
		printf("points=%lld\n", output->GetPoints()->GetNumberOfPoints());
		output->SetLines(newLines);
		int idx = output->GetPointData()->AddArray(fieldData);
		output->GetPointData()->SetActiveAttribute(idx,  vtkDataSetAttributes::SCALARS);
	}
	output->Squeeze();  // need it?
	printf("Done\n");

	return 1; // success
}





void vtkOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

