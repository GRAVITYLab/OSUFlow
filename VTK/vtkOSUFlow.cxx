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

	//if (osuflow->seedPtr == NULL) // seeds not set
	//vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkInformation *sourceInfo = inputVector[0]->GetInformationObject(0);

	//vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkDataSet *source = 0;
	if (sourceInfo)
	{
		source = vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
	} else
		return 1;


	this->SavePointInterval = this->StepLength;

	//
	//  Convert streamer into lines. Lines may be dashed.
	//
	// set data
	vtkSmartPointer<vtkDataSet> sData = source;
	osuflow->setData(sData);

#if 0
	// set seeds
	int num_seeds = source->GetNumberOfPoints();
	VECTOR3 *pSeed = new VECTOR3[num_seeds];
	int i;
	for (i=0; i < num_seeds; i++) {
		double pos[3];
		source->GetPoint(i,pos);
		pSeed[i] = VECTOR3(pos[0], pos[1], pos[2]);
	}
	osuflow->SetSeedPoints(pSeed, num_seeds);
#endif

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


	// convert traces from list-of-list to vtkPolyData
	vtkOSUFlow::StreamPoint *sPrev, *sPtr;
	//vtkFloatArray *newVectors;
	//vtkFloatArray *newScalars=NULL;
	vtkIdType ptId;
	int j;
	vtkIdList *pts;
	//double tOffset, x[3], v[3], s, r;
	//double theta;
	//vtkPolyLine* lineNormalGenerator = NULL;
	//vtkFloatArray* normals = NULL;
	//vtkFloatArray* rotation = 0;
	vtkCellArray *newLines = vtkCellArray::New();
	vtkPoints *newPts = vtkPoints::New();

	std::list<vtListSeedTrace*>::iterator pIter;
	pIter = list.begin();
	for (; pIter!=list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter;
		std::list<VECTOR3*>::iterator pnIter;
		pnIter = trace->begin();

		vtkLine *line = vtkLine::New();

		for (i=0; pnIter!=trace->end(); pnIter++, i++) {
		  VECTOR3 p = **pnIter;
		  printf(" %f %f %f ", p[0], p[1], p[2]);
		  //vtk
		  ptId = newPts->InsertNextPoint((float *)&p[0]);
		  line->GetPointIds()->SetId(i, ptId);

		}
		newLines->InsertNextCell(line);
		line->Delete();
		printf("\n");
	}




	output->SetPoints(newPts);
	newPts->Delete();
	//output->GetPointData()->SetVectors(newVectors);
	//newVectors->Delete();
	output->SetLines(newLines);
	newLines->Delete();
	//output->Squeeze();

	return 0;
}

// code reference: vtkStreamLine.cxx


void vtkOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Step Length: " << this->StepLength << "\n";

}
