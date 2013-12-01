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
#include "VectorFieldVTK.h"
#include "Blocks.h"
#include "ParFlow.h"

using namespace std;

////////////////////////////////////////////////////////
// Defines static vtkPOSUFlow::New() here
vtkStandardNewMacro(vtkPOSUFlow);

vtkPOSUFlow::vtkPOSUFlow()
: UseDIYPartition(true)
, diy_initialized(false)
, WaitFactor(0.1)
, totTime(0)
, totInTime(0)
, totOutTime(0)
, totCompCommTime(0)
, MaxRounds(100)
{
	this->extentTable = vtkTableExtentTranslator::New();
	//this->pcontroller = new POSUFlowController;
}

vtkPOSUFlow::~vtkPOSUFlow()
{
	this->extentTable->Delete();
	//delete this->pcontroller;
	if (this->diy_initialized) {
		DIY_Finalize();
		this->diy_initialized = false;
	}
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

	printf("[RequestUpdateExtent] outInfo: Piece=%d, numPieces=%d, ghostLevel=%d\n", piece, numPieces, ghostLevel);
	int wholeExtent[6], extent[6];
	outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent);

	//
	// Deterine extents
	//
	vtkExtentTranslator *translator = vtkStreamingDemandDrivenPipeline::GetExtentTranslator(outInfo);
	if (this->UseDIYPartition) {
		// init DIY

		if (this->diy_initialized) {
			DIY_Finalize();
			this->diy_initialized = false;
		}

		vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
		int nproc = controller->GetNumberOfProcesses();
		MPI_Comm comm;
		vtkMPICommunicator *vtkcomm = vtkMPICommunicator::SafeDownCast( controller->GetCommunicator() );
		if (vtkcomm)
			comm = *vtkcomm->GetMPIComm()->GetHandle();
		else {
			comm = MPI_COMM_WORLD;
			printf("communicator is not vtkMPICommunicator!  Use MPI_COMM_WORLD\n");
		}

		int dims = 4; // 4D
		int data_size[3];
		data_size[0] = wholeExtent[1]-wholeExtent[0]+1;
		data_size[1] = wholeExtent[3]-wholeExtent[2]+1;
		data_size[2] = wholeExtent[5]-wholeExtent[4]+1;
		data_size[3] = 1;

		DIY_Init(dims, data_size, 1, comm);


		int npart = nproc;
		int loc_npart = 1;
		int given[4] = {0, 0, 0, 1}; // constraints in x, y, z, t
		int ghost_list[8] = {1, 1, 1, 1, 1, 1, 0, 0}; // -x, +x, -y, +y, -z, +z, -t, +t

		int diy_did = DIY_Decompose(ROUND_ROBIN_ORDER, npart, &loc_npart, 1, ghost_list, given);
		assert(diy_did==0);


		assert(loc_npart==1); // TODO

		initExtentTableByDIY(translator);

		this->diy_initialized = true;

	} else {

		// determine all extents
		initExtentTable(translator);

	}
	this->extentTable->GetExtent(extent);


	// data
	int numInputs = this->GetNumberOfInputConnections(0);
	printf("[vtkPOSUFlow::RequestUpdateExtent] whole extent=%d %d %d %d %d %d, numInputs=%d \n",
		  wholeExtent[0], wholeExtent[1], wholeExtent[2], wholeExtent[3], wholeExtent[4], wholeExtent[5], numInputs);

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
		  //vtkStreamingDemandDrivenPipeline::SetExtentTranslator(info, this->extentTable);

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
	// init.  Get comm
	//
	int rank=0, nproc=1, npart=1;
	MPI_Comm comm = MPI_COMM_WORLD;
	{
		// get MPI controller, rank, nproc
		vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
		if (controller && controller->GetNumberOfProcesses()>1) {
			rank = controller->GetLocalProcessId();
			nproc = controller->GetNumberOfProcesses();
		} else {
			// If only one process, use sequential codes
			printf("Use serial OSUFlow\n");
			int result = vtkOSUFlow::RequestData(request,inputVector,outputVector);
			return result;
		}
		npart = nproc; // for now

		// get comm
		vtkMPICommunicator *vtkcomm = vtkMPICommunicator::SafeDownCast( controller->GetCommunicator() );
		if (vtkcomm)
			comm = *vtkcomm->GetMPIComm()->GetHandle();
		else {
			printf("communicator is not vtkMPICommunicator!\n");
		}
	}

	// get extent
	int *wholeExtent, *extent;
	wholeExtent = extentTable->GetWholeExtent();
	extent = extentTable->GetExtent();

	// Debug information
	{
		//data->PrintSelf(cout, (vtkIndent)0);
		double *bounds = data->GetBounds();
		printf("[RequestData] Rank=%d, numprocs=%d\n", rank, nproc);
		printf("[RequestData] Data bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
		if (!(bounds[0] == extent[0] && bounds[1] == extent[1] && bounds[2] == extent[2]
		          && bounds[3] == extent[3] && bounds[4]==extent[4] && bounds[5] == extent[5]))
			printf("Warning: bounds are different from extent\n");
	}

	//
	// init DIY
	//
	if (this->UseDIYPartition == false ) {

		if (this->diy_initialized) {
			DIY_Finalize();
			this->diy_initialized = false;
		}

		int dims = 4; // 4D
		int data_size[3];
		data_size[0] = wholeExtent[1]-wholeExtent[0]+1;
		data_size[1] = wholeExtent[3]-wholeExtent[2]+1;
		data_size[2] = wholeExtent[5]-wholeExtent[4]+1;
		data_size[3] = 1;
		DIY_Init(dims, data_size, 1, comm);

		// set DIY decomposition
		// - currently one processor handles one block
		// - To process multiple blocks per processor, we will use vtkMultiBlock... in the pipeline
		{
			bb_t diyBound;
			diyBound.min[0] = extent[0]; diyBound.max[0] = extent[1];
			diyBound.min[1] = extent[2]; diyBound.max[1] = extent[3];
			diyBound.min[2] = extent[4]; diyBound.max[2] = extent[5];
			vector<gb_t> neighborIdAry;
			getNeighborIds(neighborIdAry, this->extentTable, rank);
			gb_t *pNeighborId = &neighborIdAry[0];
			int num_neighborId = neighborIdAry.size();

			int diy_did = DIY_Decomposed(1, &rank, &diyBound,
					(ri_t**)NULL, (int*)NULL, (int **)NULL, (int *)NULL,
					&pNeighborId , &num_neighborId,
					0);
			assert(diy_did==0);

		}
		this->diy_initialized = true;
	}
	else
		assert(this->diy_initialized);

	//
	// OSUFlow
	//
	// set data.  Since there is only one block per process, we create one osuflow
	int loc_npart = 1;
	OSUFlow **pposuflow = new OSUFlow*[loc_npart];
	for (i=0; i<loc_npart; i++)
	{
		pposuflow[i] = new OSUFlow;
	    CVectorField *field = new VectorFieldVTK( data );
		pposuflow[i]->SetFlowField(field); // TODO
	}
    //request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1 );


	// assign seeds
	int num_seeds = input->GetNumberOfPoints();
	printf("num_seed=%d\n", num_seeds);
	if (num_seeds==0)
		return 1;
	VECTOR3 *pSeed = new VECTOR3[num_seeds];
	for (i=0; i < num_seeds; i++) {
		double pos[3];
		input->GetPoint(i,pos);
		pSeed[i] = VECTOR3(pos[0], pos[1], pos[2]);  // convert to floats
	}

	// init parflow
	printf("init parflow\n");
	int tsize = 1;  // total time steps

	list<vtListTimeSeedTrace*> *sl_list = new list<vtListTimeSeedTrace*>[loc_npart];
	Blocks *blocks = new Blocks(loc_npart, (void *)pposuflow[0], OSUFLOW, 0, 0, (DataMode)0); // 0 : load OSUFlow manually
	VECTOR4 *pt = NULL; // resulting traced points
	int *npt = NULL; // resulting number of points per trace
	int tot_ntrace;
	ParFlow *parflow = new ParFlow(blocks, (OSUFlow **)pposuflow, sl_list, &pt, &npt, &tot_ntrace, loc_npart, 0);
	parflow->SetComm(comm);

	parflow->SetMaxError(this->MaximumError);
	parflow->SetInitialStepSize(this->IntegrationStepLength);
	parflow->SetMinStepSize(this->MinimumIntegrationStep);
	parflow->SetMaxStepSize(this->MaximumIntegrationStep);
	parflow->SetIntegrationOrder((INTEG_ORD)this->IntegratorOrder);
	parflow->SetUseAdaptiveStepSize(true);
	// integrate
	TRACE_DIR dir; // not used currently
	switch (this->IntegrationDirection)
	{
	case VTK_INTEGRATE_FORWARD: dir=FORWARD_DIR; break;
	case VTK_INTEGRATE_BACKWARD: dir=BACKWARD_DIR; break;
	default: dir = BACKWARD_AND_FORWARD; printf("[vktPOSUFlow] Warning: BOTH dir is not supported yet.\n"); break;
	};
	parflow->SetIntegrationDir(dir);


	vector< vector<Particle> > Seeds; // seeds in current round for all local blocks
	Seeds.resize(loc_npart);

	// set Seeds from pSeed array.  tf is ignored.
	parflow->InitTraces(Seeds, 0, loc_npart, tsize, 1, pSeed, num_seeds);

	//delete[] pSeed;



	MPI_Barrier(comm);

	//
	// RUN
	//
	printf("Run\n");
	this->totCompCommTime = 0;
	// for all time groups
	//for (g = 0; g < ntpart; g++)
	{
		double t0;

		// check if there is any work to be done for this time group
		//if(!isSeedInTimeGroupTotal(g))
		//  continue;  // go to next time group

		// synchronize before starting I/O
		//MPI_Barrier(comm);
		//t0 = MPI_Wtime();

		// delete blocks from previous time group
		//if (g > 0) {
		//  blocks->DeleteBlocks(g_io+1, tsize, ntpart, npart);
		//}

		// todo: change seeds to vector in repartition
		// #ifdef REPARTITION
		//     parflow->Repartition(g, &npart, &Seeds, 0, &osuflow,
		// 		       block, OSUFLOW, MPI_COMM_WORLD, wgts);
		// #endif

		// inform loaded blocks for this time group
		blocks->SetLoad(0); // lid

		//g_io = g;

		// synchronize after I/O
		//MPI_Barrier(comm);
		//TotInTime += (MPI_Wtime() - t0);
		t0 = MPI_Wtime();

		// scale blocks to improve visibility
		//if (this->scale!=1.0) {
		//	for (i = 0; i < npart; i++) {
		//		if (blocks->GetLoad(i))
		//			osuflow->ScaleField(scale);
		//	}
		//}

#ifdef REPARTITION
		assert(nblocks <= MAX_BLK);
		for (i = 0; i < nblocks; i++)
		  wgts[i] = 0;
#endif

		// for all rounds
		for (j = 0; j < this->MaxRounds; j++)
		{

			printf("Round %d, Rank %d: tracing seeds: %zu\n", j, rank, Seeds[0].size());
			// for all blocks
			for (i = 0; i < loc_npart; i++)
			{

#ifdef MPE
				MPE_Log_event(compute_begin, 0, NULL);
#endif

		// compute fieldlines
				if (tsize > 1) {
#ifdef REPARTITION
					assert(i < MAX_BLK);
					parflow->ComputePathlines(Seeds[i], i, this->MaximumNumberOfSteps, this->MaximumNumberOfSteps, &wgts[i]);
#else
					parflow->ComputePathlines(Seeds[i], i, this->MaximumNumberOfSteps, this->MaximumNumberOfSteps);
#endif
				}
				else {
#ifdef REPARTITION
					assert(i < MAX_BLK);
					parflow->ComputeStreamlines(Seeds[i], i, this->MaximumNumberOfSteps, this->MaximumNumberOfSteps, &wgts[i]);
#else
					parflow->ComputeStreamlines(Seeds[i], i, this->MaximumNumberOfSteps, this->MaximumNumberOfSteps);
#endif
				}

#ifdef MPE
				MPE_Log_event(compute_end, 0, NULL);
#endif

			} // for all blocks

			// exchange neighbors
			//printf("Start exchanging\n");
			parflow->ExchangeNeighbors(Seeds, WaitFactor);
			//printf("End exchanging\n");
			MPI_Barrier(comm);

		} // for all rounds

		// flush any remaining messages
		parflow->FlushNeighbors(Seeds);

#ifdef REPARTITION
		AdvanceWeights(g);
#endif

		// end time group synchronized to get accurate timing
		MPI_Barrier(comm);
		totCompCommTime += (MPI_Wtime() - t0);

    } //for all time groups


	// gather fieldlines for rendering
#if 1
	printf("Gather traces\n");
	{
		// synchronize prior to gathering
		MPI_Barrier(comm);
		this->totOutTime = MPI_Wtime();

		float size[3];
		size[0] = wholeExtent[1]-wholeExtent[0]+1;
		size[1] = wholeExtent[3]-wholeExtent[2]+1;
		size[2] = wholeExtent[5]-wholeExtent[4]+1;
		//xxxxparflow->GatherFieldlines(loc_npart, size, tsize);

		int *ntrace = NULL; // number of traces for each proc
		int n; // total number of my points

		// gather number of points in each trace at the root => ntrace
		int all_gather = 0;  // 0: only root collects the data
		n = parflow->GatherNumPts(ntrace, all_gather, loc_npart);

		// gather the actual points in each trace at the root
		parflow->GatherPts(ntrace, n, loc_npart);

		MPI_Barrier(comm);
		this->totOutTime = MPI_Wtime() - this->totOutTime;

		if (ntrace)
			delete[] ntrace;

		parflow->PrintPerf(this->totTime, this->totInTime, this->totOutTime, this->totCompCommTime, num_seeds, size);

	}
#endif

	//
	// convert traces from pt, npt, tot_ntrace to vtkPolyData
	//
	if (rank==0)
	{
		printf("convert traces\n");
		newLines = vtkCellArray::New();
		newPts = vtkPoints::New();
		pts = vtkIdList::New();
#if 1
		// normal output
		if (pt)
		{

			int count = 0;
			printf("ntrace=%d\n", tot_ntrace);
			for (i = 0; i < tot_ntrace; i++)
			{
				//printf("length[%d]=%d\n", i, npt[i]);
				for (j = 0; j < npt[i]; j++)
				{
					//printf("(%f %f %f) ", pt[count][0], pt[count][1], pt[count][2]);
					ptId = newPts->InsertNextPoint((float *)&pt[count++][0]);
					pts->InsertNextId(ptId);
				}
				newLines->InsertNextCell(pts);
				pts->Reset();
			}
		}

#elif 0
	// test

		for (int i=0; i<input->GetNumberOfPoints(); i++)
		{
			ptId = newPts->InsertNextPoint(input->GetPoint(i));
			pts->InsertNextId(ptId);
		}
		newLines->InsertNextCell(pts);

		//output->GetPointData()->SetNormals()
#else
		// test #2

		if (rank==0) {
			for (int j=0; j<30; j++)
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
		}
		//output->GetPointData()->SetNormals()
#endif

		//
		// assign lines to output
		//
		if (newPts->GetNumberOfPoints() > 0)
		{
			output->SetPoints(newPts);
			printf("output points=%lld\n", output->GetPoints()->GetNumberOfPoints());
			output->SetLines(newLines);
		}

		output->Squeeze();

		newLines->Delete();
		newPts->Delete();
		pts->Delete();

	} // rank==0

	//
	// Clean up
	//


// segfault happens sometimes after deleting these stuff
	printf("rank:%d, cleanup\n", rank);
	if (pt)
		delete[] pt;
	if (npt)
		delete[] npt;

	for(i=0; i<loc_npart; i++)
	{
		list<vtListTimeSeedTrace*>::iterator trace_iter;
		for(trace_iter=sl_list[i].begin();trace_iter!=sl_list[i].end();trace_iter++)
		{
			vtListTimeSeedTrace::iterator pt_iter;
			for(pt_iter = (*trace_iter)->begin(); pt_iter != (*trace_iter)->end();
					pt_iter++)
				delete *pt_iter;
			(*trace_iter)->clear();
			delete *trace_iter;
		}
		sl_list[i].clear();
	}
	delete [] sl_list;
	for (i = 0; i < Seeds.size(); i++)
		Seeds[i].clear();
	Seeds.clear();

	delete blocks;
	delete parflow;

	for (i = 0; i < loc_npart; i++)
		if (pposuflow[i] != NULL)
			delete pposuflow[i];
	delete[] pposuflow;

	MPI_Barrier(comm);

	// clean up
	DIY_Finalize();
	this->diy_initialized = false;

	printf("Done\n");

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

void vtkPOSUFlow::initExtentTableByDIY(vtkExtentTranslator *translator)
{
    vtkMultiProcessController* controller = vtkMultiProcessController::GetGlobalController();
    int nproc = controller->GetNumberOfProcesses();
    int rank = controller->GetLocalProcessId();
    int npart = nproc;
    int loc_npart = 1;

	this->extentTable->SetNumberOfPieces(npart);
	this->extentTable->SetNumberOfPiecesInTable(npart);

    int wholeExtent[6], extent[6];

    translator->GetWholeExtent(wholeExtent);
    this->extentTable->SetWholeExtent(wholeExtent);

	// assign extentTable
	bb_t bounds;
	//DIY_Block_bounds(0, rank, &bounds);  // did is always 0
	//extent[0] = bounds.min[0]; extent[1] = bounds.max[0];
	//extent[2] = bounds.min[1]; extent[3] = bounds.max[1];
	//extent[4] = bounds.min[2]; extent[5] = bounds.max[2];
	int starts[3], sizes[3];
	DIY_Block_starts_sizes(0, 0, starts, sizes);
	extent[0] = starts[0]; extent[1] = starts[0] + sizes[0] -1;
	extent[2] = starts[1]; extent[3] = starts[1] + sizes[1] -1;
	extent[4] = starts[2]; extent[5] = starts[2] + sizes[2] -1;
	printf("DIY: rank=%d, Extent=%d %d %d %d %d %d\n", rank, extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);

	// sync extent
    int *gatheredExtents = new int[nproc * 6];
	controller->AllGather(extent, gatheredExtents, 6);

	// assign extentTable
	for (int cc=0; cc < nproc; cc++)
	{
		this->extentTable->SetExtentForPiece(cc, gatheredExtents + 6*cc);
	}
	this->extentTable->SetPiece(rank);
	this->extentTable->PieceToExtent();
	delete[] gatheredExtents;

#if 1
    if (rank == 0)
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
