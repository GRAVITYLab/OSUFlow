#include "OSUFlow.h"

#pragma warning(disable : 4251 4100 4244 4101)

OSUFlow::OSUFlow()
{
	bUseRandomSeeds = false;
	flowName = NULL;
	flowField = NULL;
	bStaticFlow = true;
	seedPtr = NULL; 
	nSeeds = 0; 
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

	if(bStaticFlow)
		InitStaticFlowField();
	else
		InitTimeVaryingFlowField();
}

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
	// set the boundary of physical grid
	VECTOR3 minB, maxB;
	minB.Set(0.0, 0.0, 0.0);
	maxB.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	pRegularCGrid->SetBoundary(minB, maxB);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, 1);
	if(!flowField->IsNormalized())
		flowField->NormalizeField(true);
}

void OSUFlow:: InitTimeVaryingFlowField(void)
{
	FILE *fIn;
	FILE *fVecIn;
	int timesteps;
	char* filename = new char[100];
	int dimension[3], totalNum, tStep;
	float** ppData = NULL;
	float* pData = NULL;
	VECTOR3** ppVector;
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	
	// load in data
	fIn = fopen(flowName, "r");
	assert(fIn != NULL);
	fscanf(fIn, "%d", &timesteps);
	ppData = new float*[timesteps];
	ppVector = new VECTOR3*[timesteps];
	for(int iFor = 0; iFor < timesteps; iFor++)
	{
		VECTOR3* pVector;

		fscanf(fIn, "%s", filename);
		fVecIn = fopen(filename, "rb");
		fread(dimension, sizeof(int), 3, fVecIn);
		fread(&tStep, sizeof(int), 1, fVecIn);
		totalNum = dimension[0] * dimension[1] * dimension[2];
		pData = new float[totalNum * 3];
		fread(pData, sizeof(float), totalNum*3, fVecIn);
		fclose(fVecIn);

		pVector = new VECTOR3[totalNum];
		for(int jFor = 0; jFor < totalNum; jFor++)
			pVector[jFor].Set(pData[jFor*3], pData[jFor*3+1], pData[jFor*3+2]);
		delete[] pData;
		ppVector[iFor] = pVector;
	}
	fclose(fIn);
	
	pSolution = new Solution(ppVector, totalNum, timesteps);
	pRegularCGrid = new RegularCartesianGrid(dimension[0], dimension[1], dimension[2]);
	// set the boundary of physical grid
	VECTOR3 minB, maxB;
	minB.Set(0.0, 0.0, 0.0);
	maxB.Set((float)dimension[0], (float)dimension[1], (float)dimension[2]);
	pRegularCGrid->SetBoundary(minB, maxB);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, timesteps);
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
	printf(" done\n"); 


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

// only generate one streamline from the given seed
bool OSUFlow::GenStreamLines(VECTOR3* seed, 
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
	pStreamLine->setSeedPoints(seed, seedNum, currentT);
	pStreamLine->SetInitStepSize(initialStepSize);
	pStreamLine->SetMaxStepSize(maxStepSize);
	pStreamLine->setIntegrationOrder(FOURTH);
	pStreamLine->execute((void *)&currentT, listSeedTraces);

	// release resource
	delete pStreamLine;
	return true;
}
