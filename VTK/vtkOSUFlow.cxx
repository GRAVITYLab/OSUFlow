#include "vtkOSUFlow.h"
#include "OSUFlowVTK.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"


vtkStandardNewMacro(vtkOSUFlow);

vtkOSUFlow::vtkOSUFlow()
: StepLength(.1)
{
	osuflow = new OSUFlowVTK();
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

	vtkInformation *outInfo = outputVector->GetInformationObject(0);	// data
	vtkInformation *sourceInfo = inputVector[0]->GetInformationObject(0); // traces
	vtkInformation *inInfo = inputVector[1]->GetInformationObject(0);  // seeds

	vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *source = 0;
	if (sourceInfo)
	{
		source = vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	} else {
		printf("vtkOSUFlow: Output port not connected \n");
		return 0;
	}

	// make compatible to vtkStreamer
	this->SavePointInterval = this->StepLength;

	//
	// OSUFlow
	//
	// set data
	osuflow->setData(source);

	// assign seeds
	int num_seeds = input->GetNumberOfPoints();
	printf("num_seed=%d\n", num_seeds);
	VECTOR3 *pSeed = new VECTOR3[num_seeds];  //TODO: delete
	for (i=0; i < num_seeds; i++) {
		double pos[3];
		input->GetPoint(i,pos);
		pSeed[i] = VECTOR3(pos[0], pos[1], pos[2]);
	}
	osuflow->SetSeedPoints(pSeed, num_seeds);
	printf("!!");

	// integrate
	TRACE_DIR dir;
	switch (this->IntegrationDirection)
	{
	case VTK_INTEGRATE_FORWARD: dir=FORWARD_DIR; break;
	case VTK_INTEGRATE_BACKWARD: dir=BACKWARD_DIR; break;
	default: dir = BACKWARD_AND_FORWARD; break;
	};
	list<vtListSeedTrace*> list;
	osuflow->SetIntegrationParams(this->StepLength, this->StepLength);
	osuflow->GenStreamLines(list , dir, this->GetMaximumPropagationTime(), 0);

	delete[] pSeed;

	//
	// convert traces from list-of-list to vtkPolyData
	//
	newLines = vtkCellArray::New();
	newPts = vtkPoints::New();

	std::list<vtListSeedTrace*>::iterator pIter;
	pIter = list.begin();
	for (; pIter!=list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter;
		std::list<VECTOR3*>::iterator pnIter;

		vtkLine *line = vtkLine::New();

		for (pnIter = trace->begin(); pnIter!=trace->end(); pnIter++) {
			VECTOR3 &p = **pnIter;
			//printf(" %f %f %f ", p[0], p[1], p[2]);

			//vtk
			ptId = newPts->InsertNextPoint((float *)&p[0]);
			line->GetPointIds()->InsertNextId(ptId);

			// clear up
			delete *pnIter;
		}
		newLines->InsertNextCell(line);
		line->Delete();
		//printf("\n");

		// clear up
		delete trace;
	}

	//
	// assign lines to output
	//
	output->SetPoints(newPts);
	newPts->Delete();
	//printf("points=%d\n", output->GetPoints()->GetNumberOfPoints());
	//output->GetPointData()->SetVectors(newVectors);
	//newVectors->Delete();
	output->SetLines(newLines);
	newLines->Delete();
	//output->Squeeze();
	printf("Done\n");

	return 1;
}

// code reference: vtkStreamLine.cxx


void vtkOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Step Length: " << this->StepLength << "\n";

}
