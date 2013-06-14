// code reference: vtkStreamLine.cxx

#include "vtkPOSUFlow.h"
#include "OSUFlowVTK.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkMultiProcessController.h"
#include "vtkStreamingDemandDrivenPipeline.h"

// Defines static vtkPOSUFlow::New() here
vtkStandardNewMacro(vtkPOSUFlow);



#if 0
vtkPOSUFlow* vtkPOSUFlow::New()
{
	printf("!!New!!\n");
	VTK_STANDARD_NEW_BODY(vtkPOSUFlow)
}
vtkInstantiatorNewMacro(vtkPOSUFlow)
#endif


vtkPOSUFlow::vtkPOSUFlow()
{
}

vtkPOSUFlow::~vtkPOSUFlow()
{
}


int vtkPOSUFlow::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  int piece =      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces =  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int ghostLevel = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  printf("Piece=%d, numPieces=%d, ghostLevel=%d\n", piece, numPieces, ghostLevel);

  int numInputs = this->GetNumberOfInputConnections(0);
  for (int idx = 0; idx < numInputs; ++idx)
    {
    vtkInformation *info = inputVector[0]->GetInformationObject(idx);
    if (info)
      {
      info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),           piece);
      info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),       numPieces);
      info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevel);
      }
    }

  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  if (sourceInfo)
    {
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),           1);
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),       0);
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevel);
    }

  return 1;
}


int vtkPOSUFlow::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  return 1;
}

int vtkPOSUFlow::RequestData(
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

	// Parallel
	vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
	if (controller)
		printf("Rank=%d, numprocs=%d\n", controller->GetLocalProcessId(), controller->GetNumberOfProcesses());
	else
		printf("Controller=NULL\n");

	//
	// OSUFlow
	//
	// set data
	if (source)
		osuflow->setData(source);
	else if (! osuflow->getHasData() ) {
		printf("vtkPOSUFlow: no data\n");
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
			printf("vtkPOSUFlow: No seeds\n");
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
	output->GetPointData()->SetNormals()
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





void vtkPOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}

