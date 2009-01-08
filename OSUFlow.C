

#include "OSUFlow.h"

#ifdef MPI
#include <mpi.h>
#endif

#pragma warning(disable : 4251 4100 4244 4101)

OSUFlow::OSUFlow()
{
	bUseRandomSeeds = false;
	flowName = NULL;
	flowField = NULL;
	bStaticFlow = true;
	seedPtr = NULL; 
	nSeeds = 0; 
	pStreakLine = NULL; 
}

OSUFlow::~OSUFlow()
{
	delete[] flowName;
	delete flowField;
	if (seedPtr!=NULL) delete[] seedPtr; 
	flowName = NULL;
	flowField = NULL;
}

void OSUFlow::Reset(void)
{
	bUseRandomSeeds = false;
}

void OSUFlow::LoadData(const char* fname, bool bStatic)
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;

	if(bStaticFlow) {
		numTimesteps = 1; 
		MinT = MaxT = 0; 
		InitStaticFlowField();
	}
	else
		InitTimeVaryingFlowField();
}


//sMin/sMax are local min and max range of the data that are held within 
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;

	if(bStaticFlow) {
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  InitStaticFlowField(sMin, sMax);
	}
	else
	  InitTimeVaryingFlowField(sMin, sMax); // to be implemented 
}

//sMin/sMax are local min and max range of the data that are held within 
// t_min/t_max are the time range (for time-varying field) 
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, int min_t, int max_t)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;

	if (max_t >= min_t) {
	  numTimesteps = max_t-min_t; 
	  MinT = min_t; MaxT = max_t; 
	}
	  
	if(bStaticFlow) {
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  InitStaticFlowField(sMin, sMax);
	}
	else
	  InitTimeVaryingFlowField(sMin, sMax, min_t, max_t); 
}

//---------------------------------------------------------------------------
//
/////////////////////////////////////////////////////////////////
//
//Read the whole datafile and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(void)
{
	FILE *fIn;
	int dimension[3], totalNum;
	float* pData = NULL;
	
	fIn = fopen(flowName, "rb");
	assert(fIn != NULL);
	fread(dimension, sizeof(int), 3, fIn);
	totalNum = dimension[0] * dimension[1] * dimension[2];
	pData = new float[totalNum * 3];
	fread(pData, sizeof(float), totalNum*3, fIn);
	fclose(fIn);

	// update the domain bounds 
	gMin.Set(0.0, 0.0, 0.0);
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	lMin.Set(0.0, 0.0, 0.0);
	lMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	
	// create field object
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	VECTOR3* pVector;
	VECTOR3** ppVector;
	pVector = new VECTOR3[totalNum];
	for(int iFor = 0; iFor < totalNum; iFor++)
		pVector[iFor].Set(pData[iFor*3], pData[iFor*3+1], pData[iFor*3+2]);
	delete[] pData;
	ppVector = new VECTOR3*[1];
	ppVector[0] = pVector;
	pSolution = new Solution(ppVector, totalNum, 1);
	pRegularCGrid = new RegularCartesianGrid(dimension[0], dimension[1], dimension[2]);
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, 1);
	if(!flowField->IsNormalized())
		flowField->NormalizeField(true);
}

/////////////////////////////////////////////////////////////////
//
//Read only a subdomain and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(VECTOR3 sMin, VECTOR3 sMax)
{
	FILE *fIn;
	int dimension[3], totalNum;
	float* pData = NULL;
	int lxdim, lydim, lzdim; 
	
	fIn = fopen(flowName, "rb");
	assert(fIn != NULL);
	fread(dimension, sizeof(int), 3, fIn);
	gMin.Set(0.0,0.0,0.0); 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	lMin = sMin; lMax = sMax; //local data min/max range

	lxdim = sMax[0]-sMin[0]+1; 
	lydim = sMax[1]-sMin[1]+1; 
	lzdim = sMax[2]-sMin[2]+1; 
	
	totalNum = lxdim*lydim*lzdim; 
	pData = new float[totalNum * 3];
	float *p = pData; 
	for (int z = sMin[2]; z<=sMax[2]; z++) {
	  for (int y = sMin[1]; y<=sMax[1]; y++) {
	    long offset = (z*dimension[0]*dimension[1]+y*dimension[0]+sMin[0])*3*4; 
	    fseek(fIn, offset, SEEK_SET); 
	    int size = (sMax[0]-sMin[0]+1)*3; 
	    fread(p, sizeof(float), size, fIn); 
	    p+=size; 
	  }
	}
	fclose(fIn);
        CreateStaticFlowField(pData, sMin, sMax); 
}


///////////////////////////////////////////////
//
//   Create a static flow field 
//   Input: (float) vector data, minB, maxB
//   note this flow field can be a subfield so 
//   minB[0/1/2] may not be zero, and maxB[0/1/2]
//   may not be the max of the entire field 
//   
//
void OSUFlow::CreateStaticFlowField(float* pData, VECTOR3 minB, 
				  VECTOR3 maxB)
{
	int dimension[3], totalNum;

	dimension[0] = maxB[0]-minB[0]+1; 
	dimension[1] = maxB[1]-minB[1]+1; 
	dimension[2] = maxB[2]-minB[2]+1; 

	totalNum = dimension[0] * dimension[1] * dimension[2];

	// create field object
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	VECTOR3* pVector;
	VECTOR3** ppVector;
	pVector = new VECTOR3[totalNum];

	for(int iFor = 0; iFor < totalNum; iFor++)
		pVector[iFor].Set(pData[iFor*3], pData[iFor*3+1], pData[iFor*3+2]);
	ppVector = new VECTOR3*[1];
	ppVector[0] = pVector;
	pSolution = new Solution(ppVector, totalNum, 1);
	pRegularCGrid = new RegularCartesianGrid(dimension[0], dimension[1], dimension[2]);
	lMin = minB; lMax = maxB; //local data min/max range
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, 1);
	if(!flowField->IsNormalized())
		flowField->NormalizeField(true);
}

// read the whole time sequence and create a time-varying vector field 
void OSUFlow:: InitTimeVaryingFlowField(void)
{
	FILE *fIn;
	FILE *fVecIn;
	int timesteps;
	char* filename = new char[100];
	int dimension[3], totalNum, tStep;
	//	float** ppData = NULL;
	float* pData = NULL;
	VECTOR3** ppVector;
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	
	// load in data
	fIn = fopen(flowName, "r");
	assert(fIn != NULL);
	fscanf(fIn, "%d", &timesteps);
	//	ppData = new float*[timesteps];
	ppVector = new VECTOR3*[timesteps];
	numTimesteps = timesteps; 
	int first = 1; 
	for(int iFor = 0; iFor < timesteps; iFor++)
	{
		VECTOR3* pVector;
		fscanf(fIn, "%s", filename);
		printf(" to read %s ...\n", filename); 
		fVecIn = fopen(filename, "rb");
		if (fVecIn==NULL) {
		  printf(" problem opening %s. skip\n", filename); 
		  continue; 
		}
		fread(dimension, sizeof(int), 3, fVecIn);
		//fread(&tStep, sizeof(int), 1, fVecIn);
		totalNum = dimension[0] * dimension[1] * dimension[2];
		if (first == 1) {
		  pData = new float[totalNum * 3];
		  first = 0; 
		}
		fread(pData, sizeof(float), totalNum*3, fVecIn);
		pVector = new VECTOR3[totalNum];
		for(int jFor = 0; jFor < totalNum; jFor++)
		  pVector[jFor].Set(pData[jFor*3], pData[jFor*3+1], pData[jFor*3+2]);
		ppVector[iFor] = pVector;
		fclose(fVecIn);
	}
	delete[] pData;
	fclose(fIn);

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
        lMin = gMin; lMax = gMax; 
	pSolution = new Solution(ppVector, totalNum, timesteps);
	//	pRegularCGrid = new RegularCartesianGrid(dimension[0], dimension[1], dimension[2]);
	pRegularCGrid = new RegularCartesianGrid(lMax[0]-lMin[0]+1, 
						 lMax[1]-lMin[1]+1, 
						 lMax[2]-lMin[2]+1); 
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, timesteps);
}

// read the whole time sequence within the subdomain 
void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB)
{
	FILE *fIn;
	FILE *fVecIn;
	int timesteps;
	char* filename = new char[100];
	int dimension[3], totalNum, tStep;
	float* pData = NULL;
	VECTOR3** ppVector;
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	int lxdim, lydim, lzdim; 
	
	// load in data
	fIn = fopen(flowName, "r");
	assert(fIn != NULL);
	fscanf(fIn, "%d", &timesteps);
	ppVector = new VECTOR3*[timesteps];
	numTimesteps = timesteps;

	lxdim = maxB[0]-minB[0]+1; 
	lydim = maxB[1]-minB[1]+1; 
	lzdim = maxB[2]-minB[2]+1; 
	int first = 1; 
	for(int iFor = 0; iFor < timesteps; iFor++)
	{
		VECTOR3* pVector;

		fscanf(fIn, "%s", filename);
		printf(" to read %s ...\n", filename); 
		fVecIn = fopen(filename, "rb");
		if (fVecIn==NULL) {
		  printf(" problem opening %s. skip\n", filename); 
		  continue; 
		}
		fread(dimension, sizeof(int), 3, fVecIn);
		totalNum = lxdim*lydim*lzdim; 
		if (first == 1) {
		  pData = new float[totalNum * 3];
		  first = 0; 
		}
		float *p = pData; 
		for (int z = minB[2]; z<=maxB[2]; z++) {
		  for (int y = minB[1]; y<=maxB[1]; y++) {
		    long offset = (z*dimension[0]*dimension[1]+y*dimension[0]+minB[0])*3*4;
		    fseek(fVecIn, offset, SEEK_SET);
		    int size = (maxB[0]-minB[0]+1)*3;
		    fread(p, sizeof(float), size, fVecIn);
		    p+=size;
		  }
		}
		fclose(fVecIn);
		pVector = new VECTOR3[totalNum];
		for(int jFor = 0; jFor < totalNum; jFor++)
			pVector[jFor].Set(pData[jFor*3], pData[jFor*3+1], pData[jFor*3+2]);
		ppVector[iFor] = pVector;
	}
	delete[] pData;
	fclose(fIn);

	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	lMin = minB; lMax = maxB; //local data min/max range
	
	pSolution = new Solution(ppVector, totalNum, timesteps);
	pRegularCGrid = new RegularCartesianGrid(lMax[0]-lMin[0]+1, 
						 lMax[1]-lMin[1]+1, 
						 lMax[2]-lMin[2]+1); 
	// set the boundary of physical grid
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, timesteps);
}


void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB, int min_t, int max_t)
{
	FILE *fIn;
	FILE *fVecIn;
	int timesteps;
	char* filename = new char[100];
	int dimension[3], totalNum, tStep;
	float* pData = NULL;
	VECTOR3** ppVector;
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	int lxdim, lydim, lzdim; 
	
	// load in data
	fIn = fopen(flowName, "r");
	assert(fIn != NULL);
	fscanf(fIn, "%d", &timesteps);

	numTimesteps = max_t-min_t+1; 

	ppVector = new VECTOR3*[numTimesteps];

	lxdim = maxB[0]-minB[0]+1; 
	lydim = maxB[1]-minB[1]+1; 
	lzdim = maxB[2]-minB[2]+1; 
	int first = 1; 
	for(int iFor = 0; iFor < timesteps; iFor++)
	{
		VECTOR3* pVector;

		fscanf(fIn, "%s", filename);
		if (iFor <min_t || iFor >max_t) 
		  continue; 
		printf(" to read %s ...\n", filename); 
		fVecIn = fopen(filename, "rb");
		if (fVecIn==NULL) {
		  printf(" problem opening %s. skip\n", filename); 
		  continue; 
		}
		fread(dimension, sizeof(int), 3, fVecIn);
		//fread(&tStep, sizeof(int), 1, fVecIn);
		totalNum = lxdim*lydim*lzdim; 
		if (first == 1) {
		  pData = new float[totalNum * 3];
		  first = 0; 
		}
		float *p = pData; 
		for (int z = minB[2]; z<=maxB[2]; z++) {
		  for (int y = minB[1]; y<=maxB[1]; y++) {
		    long offset = (z*dimension[0]*dimension[1]+y*dimension[0]+minB[0])*3*4;
		    fseek(fVecIn, offset, SEEK_SET);
		    int size = (maxB[0]-minB[0]+1)*3;
		    fread(p, sizeof(float), size, fVecIn);
		    p+=size;
		  }
		}
		fclose(fVecIn);
		pVector = new VECTOR3[totalNum];
		for(int jFor = 0; jFor < totalNum; jFor++)
			pVector[jFor].Set(pData[jFor*3], pData[jFor*3+1], pData[jFor*3+2]);

		ppVector[iFor-min_t] = pVector;
	}
	fclose(fIn);
	delete[] pData;

	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	lMin = minB; lMax = maxB; //local data min/max range
	MinT = min_t; MaxT = max_t; 

	pSolution = new Solution(ppVector, totalNum, numTimesteps, min_t, max_t);
	pRegularCGrid = new RegularCartesianGrid(lMax[0]-lMin[0]+1,   
						 lMax[1]-lMin[1]+1, 
						 lMax[2]-lMin[2]+1); 
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, numTimesteps);
}
//////////////////////////////////////////////////////////////////////////
// specify a set of seed points randomly generated over the specified
// spatial interval. Points can be in axis aligned dimension 0, 1, 2, 3
//////////////////////////////////////////////////////////////////////////
void OSUFlow::SetRandomSeedPoints(const float min[3], 
				  const float max[3], 
				  int num)
{
	for(int iFor = 0; iFor < 3; iFor++)
	{
		minRakeExt[iFor] = min[iFor];
		maxRakeExt[iFor] = max[iFor];
	}
	this->numSeeds[0] = num;
	this->numSeeds[1] = 1;
	this->numSeeds[2] = 1;

	bUseRandomSeeds = true;

	// generate seeds
	nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];

	if (seedPtr!=NULL) delete[]seedPtr; 
	seedPtr = new VECTOR3[nSeeds];

	size_t SeedSize[3]; 
	SeedSize[0] = numSeeds[0]; 
	SeedSize[1] = numSeeds[1]; 
	SeedSize[2] = numSeeds[2]; 

	SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
							  (const float*)maxRakeExt, 
							  (const size_t*)SeedSize); 
	pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	delete pSeedGenerator;

}

//////////////////////////////////////////////////////////////////////////
// specify a set of seed points regularly generated over the specified
// spatial interval. Points can be in axis aligned dimension 0, 1, 2, 3
//////////////////////////////////////////////////////////////////////////
void OSUFlow::SetRegularSeedPoints(const float min[3], 
								   const float max[3], 
								   const size_t numSeeds[3])
{
	for(int iFor = 0; iFor < 3; iFor++)
	{
		minRakeExt[iFor] = min[iFor];
		maxRakeExt[iFor] = max[iFor];
	}
	this->numSeeds[0] = (unsigned int)numSeeds[0];
	this->numSeeds[1] = (unsigned int)numSeeds[1];
	this->numSeeds[2] = (unsigned int)numSeeds[2];

	bUseRandomSeeds = false;

	// generate seeds
	nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
	seedPtr = new VECTOR3[nSeeds];
	if (seedPtr!=NULL) delete[]seedPtr; 

	SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
							  (const float*)maxRakeExt, 
							  (const size_t*)numSeeds);
	pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	delete pSeedGenerator;
}

void OSUFlow::SetIntegrationParams(float initStepSize, float maxStepSize)
{
	initialStepSize = initStepSize;
	this->maxStepSize = maxStepSize;
}

//////////////////////////////////////////////////////////////////////////
// generate streamlines
// input
// listSeedTraces: STL list to keep the advection result
// maxPoints: how many advection points each streamline
// randomSeed: seed for random number generator
//////////////////////////////////////////////////////////////////////////
bool OSUFlow::GenStreamLines(list<vtListSeedTrace*>& listSeedTraces, 
			     TRACE_DIR traceDir,
			     int maxPoints,
			     unsigned int randomSeed)
{
	// first generate seeds

        if (seedPtr==NULL)  {
	  nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
	  seedPtr = new VECTOR3[nSeeds];
	  SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
							    (const float*)maxRakeExt, 
							    (const size_t*)numSeeds);
	  pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	  delete pSeedGenerator;
	}

	listSeedTraces.clear();

	// execute streamline
	vtCStreamLine* pStreamLine;
	float currentT = 0.0;
	pStreamLine = new vtCStreamLine(flowField);
	switch(traceDir)
	{
	case BACKWARD_DIR:
		pStreamLine->setForwardTracing(false);
		break;
	case FORWARD_DIR:
		pStreamLine->setBackwardTracing(false);
		break;
	case BACKWARD_AND_FORWARD:
		break;
	}
	pStreamLine->SetLowerUpperAngle(3.0, 15.0);
	pStreamLine->setMaxPoints(maxPoints);
	pStreamLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(FOURTH);
	pStreamLine->execute((void *)&currentT, listSeedTraces);
	// release resource
	delete pStreamLine;
	return true;
}




/////////////////////////////////////////////////////////
// Generate streamlines from the given seeds
bool OSUFlow::GenStreamLines(VECTOR3* seeds, 
			     TRACE_DIR traceDir,
			     const int seedNum,
			     const int maxPoints, 
			     list<vtListSeedTrace*>& listSeedTraces)
{
	// execute streamline
	vtCStreamLine* pStreamLine;
	float currentT = 0.0;
	pStreamLine = new vtCStreamLine(flowField);
	switch(traceDir)
	{
	case BACKWARD_DIR:
		pStreamLine->setForwardTracing(false);
		break;
	case FORWARD_DIR:
		pStreamLine->setBackwardTracing(false);
		break;
	case BACKWARD_AND_FORWARD:
		break;
	}
	pStreamLine->SetLowerUpperAngle(3.0, 15.0);
	pStreamLine->setMaxPoints(maxPoints);
	pStreamLine->setSeedPoints(seeds, seedNum, currentT);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(FOURTH);
	pStreamLine->execute((void *)&currentT, listSeedTraces);
	// release resource
	delete pStreamLine;
	return true;
}

//////////////////////////////////////////////////////////////////////////////

bool OSUFlow::GenPathLines(list<vtListSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			     int maxPoints,
			     unsigned int randomSeed)
{
	// first generate seeds

        if (seedPtr==NULL)  {
	  nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
	  seedPtr = new VECTOR3[nSeeds];
	  SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
							    (const float*)maxRakeExt, 
							    (const size_t*)numSeeds);
	  pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
	  delete pSeedGenerator;
	}

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;
	float currentT = 0.0;
	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 

	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute((void *)&currentT, listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}


bool OSUFlow::GenStreakLines(vtStreakTraces& streakTraces, TIME_DIR dir, 
			     float current_time, bool is_existing)
{

  if (seedPtr==NULL)  {
    nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
    seedPtr = new VECTOR3[nSeeds];
    SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
						      (const float*)maxRakeExt, 
						      (const size_t*)numSeeds);
    pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
    delete pSeedGenerator;
  }
  if (is_existing == false) {
    if (pStreakLine!=NULL) delete pStreakLine; 
    pStreakLine = new vtCStreakLine(flowField); 
  }
  //otherwise one sterakline has already been created before 
  float currentT = current_time; 
  
  pStreakLine->SetTimeDir(dir); 
  pStreakLine->SetLowerUpperAngle(3.0, 15.0);
  pStreakLine->setSeedPoints(seedPtr, nSeeds, currentT);
  pStreakLine->SetInitStepSize(initialStepSize);
  pStreakLine->SetMaxStepSize(maxStepSize);
  pStreakLine->setIntegrationOrder(FOURTH);

  pStreakLine->execute((void*) &currentT, streakTraces); 

  printf(" back from streakline exec\n"); 

} 


//------------------------------------------------------------------------
//
// Error()
// mpi error handler
//
void Error(const char *fmt, ...){

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  sleep(5);
#ifdef MPI
  MPI_Abort(MPI_COMM_WORLD, 0);
#else
  exit(0);
#endif

}
//-----------------------------------------------------------------------

// MPI functions

#ifdef MPI

//-----------------------------------------------------------------------
//
// ReadData()
//
// collectively reads the dataset
//
// sMin, sMax are local subdomain min and max
// dim is the total size of the domain
// pad is the number of ghost cells (per side)
//
// Tom Peterka, 11/24/08
//
void OSUFlow::ReadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim) {

  flowName = new char[255];
  strcpy(flowName, fname);

  bStaticFlow = bStatic;

  if(bStaticFlow) {
    MinT = MaxT = 0;
    ReadStaticFlowField(sMin, sMax, dim);
  }

//   else
//     InitTimeVaryingFlowField(); // to be implemented 
}
//---------------------------------------------------------------------------
//
// ReadStaticFlowField
//
// Read data collectively and create a vectorfield object based on my subdomain
// used MPI-IO collectives to perform the I/O
// reads a single subdomain
//
// sMin, sMax: corners of subdomain
// dim: size of total domain
//
// Tom Peterka, 11/24/08
//
void OSUFlow::ReadStaticFlowField(VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim) {

  float* pData = NULL; // the data
  MPI_File fd;
  MPI_Datatype filetype;
  int size[3]; // sizes of entire dataset
  int subsize[3]; // sizes of my subdomain
  int start[3]; // starting indices of my subdomain
  MPI_Status status;
  int err; // error status
  int rank;
  int i;

  // init	
  gMin.Set(0.0,0.0,0.0); 
  gMax.Set((float)(dim[0] - 1), (float)(dim[1] - 1), (float)(dim[2] - 1));
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // open the file
  if (MPI_File_open(MPI_COMM_WORLD,flowName,MPI_MODE_RDONLY,MPI_INFO_NULL,&fd) != MPI_SUCCESS)
    Error("Error: ReadRawData() cannot open %s\n",flowName);

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // sMin and sMax are [x][y][z]
  for (i = 0; i < 3; i++) {
    size[2 - i] = dim[i];
    start[2 - i] = sMin[i];
    subsize[2 - i] = sMax[i] - sMin[i] + 1;
  }

  // allocate data space
  pData = new float[3 * subsize[0] * subsize[1] * subsize[2]];
  if (!pData)
    Error("Error: ReadStaticFlowField() unable to allocate data space\n");

  // debug
//   fprintf(stderr,"rank = %d dim = %.0f %.0f %.0f sMin = %.0f %.0f %.0f sMax = %.0f %.0f %.0f size = %d %d %d start = %d %d %d subsize = %d %d %d\n",rank,dim[0],dim[1],dim[2],sMin[0],sMin[1],sMin[2],sMax[0],sMax[1],sMax[2],size[2],size[1],size[0],start[2],start[1],start[0],subsize[2],subsize[1],subsize[0]);

  // do the actual collective io
  MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
			   MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);
  MPI_File_set_view(fd, 0, MPI_FLOAT, filetype, (char *)"native", 
		    MPI_INFO_NULL);
  err = MPI_File_read_all(fd, pData,
			  subsize[0] * subsize[1] * subsize[2], MPI_FLOAT, &status);
  if (err != MPI_SUCCESS)
    Error("Error: ReadStaticFlowField() rank %d error reading file\n", rank);

  // check the count
  if (status.count != sizeof(float) * subsize[0] * subsize[1] * subsize[2])
    Error("Error: ReadStaticFlowField() error rank %d read %d bytes instead of %d bytes from file\n",
	  rank, status.count, sizeof(float) * subsize[0] * subsize[1] * subsize[2]);

  //swap bytes
#ifdef BYTE_SWAP
  for (i = 0; i < subsize[0] * subsize[1] * subsize[2]; i++)
    swap4((char *)&(pData[i]));
#endif

  MPI_File_close(&fd);

  // create the field
  CreateStaticFlowField(pData, sMin, sMax); 

}
//---------------------------------------------------------------------------

#endif

// utility functions

//---------------------------------------------------------------------------
//
// swap4(n)
//
// Swaps 4 bytes from 1-2-3-4 to 4-3-2-1 order.
// cast the input as a char and use on any 4 byte variable
//
void swap4(char *n) {

  char *n1;
  char c;

  n1 = n + 3;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//----------------------------------------------------------------------------
