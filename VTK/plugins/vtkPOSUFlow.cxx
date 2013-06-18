// code reference: vtkStreamLine.cxx

#include <iostream>
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkMultiProcessController.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkExtentTranslator.h"
#include "vtkMPICommunicator.h"
#include "vtkMPI.h"
#include "vtkTableExtentTranslator.h"

#include "vtkPOSUFlow.h"
#include "OSUFlowVTK.h"

using namespace std;

////////////////////////////////////////////////////////
// Defines static vtkPOSUFlow::New() here
vtkStandardNewMacro(vtkPOSUFlow);

vtkPOSUFlow::vtkPOSUFlow()
{
	this->extentTable = vtkTableExtentTranslator::New();
	//this->pcontroller = new POSUFlowController;
}

vtkPOSUFlow::~vtkPOSUFlow()
{
	this->extentTable->Delete();
	//delete this->pcontroller;
}

// message sent downstream to Paraview
int vtkPOSUFlow::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

#if 0
  // test
  int wholeExtent[6];
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);
  printf("[vtkPOSUFlow::RequestInformation] whole extent=%d %d %d %d %d %d \n",
		  wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5]);

	// get local extent
	vtkExtentTranslator *translator = vtkStreamingDemandDrivenPipeline::GetExtentTranslator(inInfo);

	  translator->GetExtent(extent);
	  printf("* current extent: %d %d %d %d %d %d\n", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
#endif

  return 1;
}


// Message sent upstream to reader
int vtkPOSUFlow::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  int piece =      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numPieces =  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int ghostLevel = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  int pieceId =    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());

  printf("Piece=%d, numPieces=%d, ghostLevel=%d\n", piece, numPieces, ghostLevel);
  int wholeExtent[6], extent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);

  // data
  int numInputs = this->GetNumberOfInputConnections(0);
  printf("[vtkPOSUFlow::RequestUpdateExtent] whole extent=%d %d %d %d %d %d, numInputs=%d \n",
		  wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5], numInputs);

#if 1
  // determine all extents
  vtkExtentTranslator *translator = vtkStreamingDemandDrivenPipeline::GetExtentTranslator(outInfo);

  initExtentTable(translator);

  translator->Delete();

  this->extentTable->GetExtent(extent);
#endif

  for (int idx = 0; idx < numInputs; ++idx)
  {
    vtkInformation *info = inputVector[0]->GetInformationObject(idx);
    if (info)
    {
      // not needed
      //info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),           piece);
      //info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),       numPieces);
      //info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevel);

      info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);
      vtkStreamingDemandDrivenPipeline::SetExtentTranslator(info, this->extentTable);

      //info->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT_INITIALIZED(), 1); // not needed
    }
  }

  // seed
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  if (sourceInfo)
    {
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),           0);
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),       1);
    sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevel);
  }


  return 1;
}


int vtkPOSUFlow::RequestData(
	vtkInformation *request,
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{

	//
	// process inputs
	//
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

	vtkInformation *dataInfo = inputVector[0]->GetInformationObject(0); // data
	vtkInformation *inInfo = inputVector[1]->GetInformationObject(0);  // seeds
	vtkInformation *outInfo = outputVector->GetInformationObject(0);	// trace

	vtkDataSet *data = 0;
	if (dataInfo) {
		data = vtkDataSet::SafeDownCast(dataInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkDataSet *input = 0;
	if (inInfo) {
		input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	}
	vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));


	//
	// init
	//

	// get MPI controller, rank, nproc
	int rank, nproc, npart;
	vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
	if (controller && controller->GetNumberOfProcesses()>1) {
		rank = controller->GetLocalProcessId();
		nproc = controller->GetNumberOfProcesses();
	} else {
		// If only one process, use sequential codes
		int result = vtkOSUFlow::RequestData(request,inputVector,outputVector);
		return result;
	}
	npart = nproc;

	// Debug information
	//data->PrintSelf(cout, (vtkIndent)0);
	double *bounds = data->GetBounds();
	printf("[RequestData] Rank=%d, numprocs=%d\n", rank, nproc);
	printf("[RequestData] Data bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	//
	// init DIY
	//
	MPI_Comm comm;
	vtkMPICommunicator *vtkcomm = vtkMPICommunicator::SafeDownCast( controller->GetCommunicator() );
	if (vtkcomm)
		comm = *vtkcomm->GetMPIComm()->GetHandle();
	else {
		comm = MPI_COMM_WORLD;
		printf("communicator is not vtkMPICommunicator.  Use MPI_COMM_WORLD\n");
	}

	// determine extent
	int *wholeExtent, *extent;
	wholeExtent = extentTable->GetWholeExtent();
	printf("nPieces = %d, Piece = %d,  rank=%d\n", extentTable->GetNumberOfPieces(),  extentTable->GetPiece(), rank);


	int dims = 4; // 4D
	DIY_Init(dims, wholeExtent, nproc, comm);

	// set DIY decomposition
	// - currently one processor handles one block
	// - To process multiple blocks per processor, we will use vtkMultiBlock... in the pipeline
	{
		bb_t diyBound;
		diyBound.min[0] = bounds[0]; diyBound.max[0] = bounds[1];
		diyBound.min[1] = bounds[2]; diyBound.max[1] = bounds[3];
		diyBound.min[2] = bounds[4]; diyBound.max[2] = bounds[5];
		int glo_npart = nproc;
		vector<gb_t> neighborIdAry;
		getNeighborIds(neighborIdAry, this->extentTable, rank);
		gb_t *pNeighborId = &neighborIdAry[0];
		int num_neighborId = neighborIdAry.size();

		DIY_Decomposed(1, &rank, &diyBound,
				(ri_t**)NULL, (int*)NULL, (int **)NULL, (int *)NULL,
				&pNeighborId , &num_neighborId,
				0);

	}

	//
	// OSUFlow
	//
	// set data.  Since there is only one block per process, we create one osuflow
	OSUFlowVTK *osuflow = new OSUFlowVTK;
	osuflow->setData(data);
	{
		//Blocks *blocks = new Blocks(npart, (void *)osuflow, OSUFLOW, )
		//ParFlow *parflow = new ParFlow()
	}


    //request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );



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

#if 0
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
#elif 1
	// test
	newLines = vtkCellArray::New();
	newPts = vtkPoints::New();
	pts = vtkIdList::New();

	for (int i=0; i<input->GetNumberOfPoints(); i++)
	{
		ptId = newPts->InsertNextPoint(input->GetPoint(i));
		pts->InsertNextId(ptId);
	}
	newLines->InsertNextCell(pts);

	//output->GetPointData()->SetNormals()
#else
	// test
	newLines = vtkCellArray::New();
	newPts = vtkPoints::New();
	pts = vtkIdList::New();

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
	//output->GetPointData()->SetNormals()
#endif

	//
	// assign lines to output
	//
	if (newPts->GetNumberOfPoints() > 0)
	{
		output->SetPoints(newPts);
		newPts->Delete();
		printf("output points=%lld\n", output->GetPoints()->GetNumberOfPoints());
		output->SetLines(newLines);
		newLines->Delete();
	}
	output->Squeeze();  // need it?
	printf("Done\n");

	// clean up
	delete osuflow;

	// why deleting causes segfault?
	//newLines->Delete();
	//newPts->Delete();
	//pts->Delete();

	return 1; // success
}





void vtkPOSUFlow::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}


// protected:

void vtkPOSUFlow::initExtentTable(vtkExtentTranslator *translator)
{
    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    int pieces = translator->GetNumberOfPieces();
    int nproc = controller->GetNumberOfProcesses();
    int piece = translator->GetPiece();

	this->extentTable->SetNumberOfPieces(pieces);
	this->extentTable->SetNumberOfPiecesInTable(pieces);

    int *gatheredExtents = new int[nproc * 6];
    int extent[6], wholeExtent[6];
    translator->GetExtent(extent);
    translator->GetWholeExtent(wholeExtent);

    this->extentTable->SetWholeExtent(wholeExtent);

	int ghostlevel = translator->GetGhostLevel();
    //this->extentTable->SetMaximumGhostLevel() SetGhostLevel(ghostlevel);

    // determine extent with ghost cells
	if (ghostlevel<2) {
		extent[0] = max(extent[0]-1, wholeExtent[0]);
	    extent[1] = min(extent[1]+1, wholeExtent[1]);
	    extent[2] = max(extent[2]-1, wholeExtent[2]);
	    extent[3] = min(extent[3]+1, wholeExtent[3]);
	    extent[4] = max(extent[4]-1, wholeExtent[4]);
	    extent[5] = min(extent[5]+1, wholeExtent[5]);
	}

	// sync extent
	controller->AllGather(extent, gatheredExtents, 6);

	// assign extentTable
	for (int cc=0; cc < nproc; cc++)
    {
    	this->extentTable->SetExtentForPiece(cc, gatheredExtents + 6*cc);
    }
	this->extentTable->SetPiece(piece);
	this->extentTable->PieceToExtent();

    delete[] gatheredExtents;

#if 1
    if (controller->GetLocalProcessId() == 0)
    {
    	this->extentTable->Print(cout);
    }
#endif

}

void vtkPOSUFlow::getNeighborIds(vector<gb_t> &neighborIdAry, vtkExtentTranslator *translator, int rank)
{
	assert(translator->GetPiece()==rank);

	int *extent = translator->GetExtent();
	for (int i=0; i<translator->GetNumberOfPieces(); i++)
	{
		if (i==rank) continue;
		translator->SetPiece(i);
		int *tmp_extent = translator->GetExtent();
		int k;
		for (k=0; k<3; k++)
		{
			int minb=extent[k*2];
			int maxb=extent[k*2+1];
			int tmp_minb=tmp_extent[k*2];
			int tmp_maxb=tmp_extent[k*2+1];
			//  printf("* current extent: %d %d %d %d %d %d\n", tmp_extent[0], tmp_extent[1], tmp_extent[2], tmp_extent[3], tmp_extent[4], tmp_extent[5]);
			if (minb > tmp_maxb || maxb < tmp_minb)
				break;
		}
		if (k==3) // overlaps
		{
			gb_t gb;
			gb.gid = gb.proc = i;
			neighborIdAry.push_back(gb);
			printf("neighbor: %d->%d\n", rank, i);
		}
	}

	translator->SetPiece(rank);
}
