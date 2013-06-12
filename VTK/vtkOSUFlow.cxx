// code reference: vtkStreamLine.cxx

#include "vtkOSUFlow.h"
#include "OSUFlowVTK.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"


vtkStandardNewMacro(vtkOSUFlow);

vtkOSUFlow::vtkOSUFlow()
{
	osuflow = new OSUFlowVTK();
	this->SetMaximumPropagationTime(500);
}

vtkOSUFlow::~vtkOSUFlow()
{
	delete osuflow;
}

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
	vtkIdList *pts;
	//vtkPolyLine* lineNormalGenerator = NULL;
	//vtkFloatArray* normals = NULL;
	//vtkFloatArray* rotation = 0;
	vtkCellArray *newLines;
	vtkPoints *newPts;

	vtkInformation *sourceInfo = inputVector[0]->GetInformationObject(0); // data
	vtkInformation *inInfo = inputVector[1]->GetInformationObject(0);  // seeds
	vtkInformation *outInfo = outputVector->GetInformationObject(0);	// trace

	vtkDataSet *source = 0;
	if (sourceInfo) {
		source = vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkDataSet *input = 0;
	if (inInfo) {
		input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	// make compatible to vtkStreamer
	//this->SavePointInterval = this->IntegrationStepLength;
	// Here we store every point

	//
	// OSUFlow
	//
	// set data
	if (source)
		osuflow->setData(source);
	else if (! osuflow->getHasData() ) {
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
	osuflow->SetMaxError(this->MaximumError);
	osuflow->SetIntegrationParams(this->IntegrationStepLength, this->MinimumIntegrationStep, this->MaximumIntegrationStep);
	osuflow->GenStreamLines(list , dir, this->MaximumPropagationTime, 0); // default: RK45

	delete[] pSeed;

	//
	// convert traces from list-of-list to vtkPolyData
	//
	newLines = vtkCellArray::New();
	newPts = vtkPoints::New();
	pts = vtkIdList::New();

#if 1
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
#endif

	//
	// assign lines to output
	//
	if (newPts->GetNumberOfPoints() > 0)
	{
		output->SetPoints(newPts);
		newPts->Delete();
		printf("points=%lld\n", output->GetPoints()->GetNumberOfPoints());
		output->SetLines(newLines);
		newLines->Delete();
	}
	output->Squeeze();  // need it?
	printf("Done\n");

	return 1; // success
}



void vtkOSUFlow::SetIntegratorOrder(int order) {	this->integratorOrder = (INTEG_ORD)order; osuflow->SetIntegrationOrder((INTEG_ORD)order); }
int vtkOSUFlow::GetIntegratorOrder() { return (int)this->integratorOrder; }



void vtkOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

