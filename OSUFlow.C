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
	seedTimeArray = NULL; 
	nSeeds = 0; 
	pStreakLine = NULL; 
	has_data = false; 
	deferred_load_case = -1; 
}

OSUFlow::~OSUFlow()
{
  delete[] flowName;
  delete flowField;
  if (seedPtr!=NULL) delete[] seedPtr; 
  if (seedTimeArray!=NULL) delete [] seedTimeArray; 
  flowName = NULL;
  flowField = NULL;
}

void OSUFlow::DeleteData(void)
{
	delete flowField; 
	has_data = false; 
}

///////////////////////////////////////////////////////////
// 
// Load the whole static or time-varying data set
//
void OSUFlow::LoadData(const char* fname, bool bStatic, bool deferred)
{
	flowName = new char[255];
	strcpy(flowName, fname);
	bStaticFlow = bStatic;
	
	has_data = false; 

	if(bStaticFlow) {
		numTimesteps = 1; 
		MinT = MaxT = 0; 
		if (deferred == true) {
		  deferred_load_case = 0; 
		  return; 
		}
		InitStaticFlowField();
	}
	else {
		if (deferred == true) {
		  deferred_load_case = 0; 
		  return; 
		}
		InitTimeVaryingFlowField();
	}
	has_data = true; 
}

/////////////////////////////////////////////////////////////
//
//   Read a partial data set 
//   sMin/sMax are local min and max range of the data that are held 
//
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, bool deferred)  
{
	flowName = new char[255];
	strcpy(flowName, fname);
	bStaticFlow = bStatic;
	lMin = sMin; lMax = sMax; 
	has_data = false; 

	if(bStaticFlow) {
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  if (deferred == true) {
	    deferred_load_case = 1; 
	    return; 
	  }
	  InitStaticFlowField(sMin, sMax);
	}
	else {
	  if (deferred == true) {
	    deferred_load_case = 1; 
	    return; 
	  }
	  InitTimeVaryingFlowField(sMin, sMax); 
	}
	has_data = true; 
}

/////////////////////////////////////////////////////////////
//
//  Load a partial time-varying data set 
//  sMin/sMax are local min and max range of the data that are held
//  t_min/t_max are the time range (for time-varying field) 
//
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, int min_t, int max_t, 
		       bool deferred)  
{
	flowName = new char[255];
	strcpy(flowName, fname);

	bStaticFlow = bStatic;
	lMin = sMin; lMax = sMax; 
	has_data = false; 

	if (max_t >= min_t) {
	  numTimesteps = max_t-min_t+1; 
	  MinT = min_t; MaxT = max_t; 
	}
	else {   //exception. goes back to default 
	  numTimesteps = 1; 
	  MinT = MaxT = min_t; 
	}
	  
	if(bStaticFlow) {  // ignore the time range 
	  numTimesteps = 1; 
	  MinT = MaxT = 0; 
	  if (deferred == true) {
	    deferred_load_case = 2; 
	    return; 
	  }
	  InitStaticFlowField(sMin, sMax);
	}
	else {
	  if (deferred == true) {
	    deferred_load_case = 2; 
	    return; 
	  }
	  InitTimeVaryingFlowField(sMin, sMax, min_t, max_t); 
	}
	has_data = true; 
}


bool OSUFlow::DeferredLoadData() 
{
  if (deferred_load_case == -1) return(false); 

  switch(deferred_load_case) {
  case 0: 
    if(bStaticFlow) {
      InitStaticFlowField();
    }
    else
      InitTimeVaryingFlowField();
    has_data = true; 
    break; 
  case 1: 
    if(bStaticFlow) {
      InitStaticFlowField(lMin, lMax);
    }
    else
      InitTimeVaryingFlowField(lMin, lMax); 
    has_data = true; 
    break; 
  case 2: 
    if(bStaticFlow) {  // ignore the time range 
      InitStaticFlowField(lMin, lMax);
    }
    else
      InitTimeVaryingFlowField(lMin, lMax, MinT, MaxT); 
    has_data = true; 
    break; 
  }
  ScaleField(10.0); 
  return(true); 
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
	delete [] pData; 
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
	delete []pVector; 
	delete[] ppVector; 
}

///////////////////////////////////////////////////////////////
//
//  This function is used when the data dimensions are independent of 
//  the physical range (minB, maxB) 
//
void OSUFlow::CreateStaticFlowField(float* pData, int xsize, int ysize, int zsize, 
				    VECTOR3 minB, VECTOR3 maxB)
{
	int dimension[3], totalNum;

	dimension[0] = xsize; 
	dimension[1] = ysize; 
	dimension[2] = zsize; 

	totalNum = dimension[0] * dimension[1] * dimension[2];

	// create field object
	Solution* pSolution;
	RegularCartesianGrid* pRegularCGrid;
	VECTOR3* pVector;
	VECTOR3** ppVector;
	pVector = new VECTOR3[totalNum];

	for(int iFor = 0; iFor < totalNum; iFor++)
		pVector[iFor].Set(pData[iFor*3], pData[iFor*3+1], pData[iFor*3+2]);
	delete [] pData; 
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
	delete []pVector; 
	delete[] ppVector; 
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
	pRegularCGrid = new RegularCartesianGrid(lMax[0]-lMin[0]+1, 
						 lMax[1]-lMin[1]+1, 
						 lMax[2]-lMin[2]+1); 
	pRegularCGrid->SetBoundary(lMin, lMax);
	assert(pSolution != NULL && pRegularCGrid != NULL);
	flowField = new CVectorField(pRegularCGrid, pSolution, timesteps);
	for (int i=0; i<timesteps; i++) {
	  delete [] ppVector[i]; 
	}
	delete[] ppVector; 
}

/////////////////////////////////////////////////////////////////
//
// read the whole time sequence within the subdomain 
//
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
	for (int i=0; i<timesteps; i++) {
	  delete [] ppVector[i]; 
	}
	delete[] ppVector; 
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
	flowField = new CVectorField(pRegularCGrid, pSolution, numTimesteps, MinT);
	for (int i=0; i<numTimesteps; i++) {
	  delete [] ppVector[i]; 
	}
	delete[] ppVector; 
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
  
  if (has_data == false) DeferredLoadData(); 

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
	// or, seeds have already been generated previously, do nothing. 

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

  if (has_data == false) DeferredLoadData(); 

//   if (seedPtr!=NULL) delete [] seedPtr; 
        nSeeds = seedNum; 
	seedPtr = seeds; 

	listSeedTraces.clear();

	// execute streamline
	vtCStreamLine* pStreamLine;
	float currentT = 0.0;  // always starts from the default time 0 
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

//////////////////////////////////////////////////////////////////////////////
//   all particles start at the same time: currentT 
bool OSUFlow::GenPathLines(list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int maxPoints,
			   float currentT)
{

  if (has_data == false) DeferredLoadData(); 

	// first generate seeds if not exist before 
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

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 

	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
// Take an input list of seeds, which all start at the same time currentT 
bool OSUFlow::GenPathLines(VECTOR3* seeds, list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints,
			   float currentT)
{

  if (has_data == false) DeferredLoadData(); 
  if (seedPtr!=NULL) delete [] seedPtr; 

        seedPtr = seeds; 
	nSeeds = num_seeds; 

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 

	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, currentT);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
//    Take an input list of seeds, which 
//    can start from different times (tarray) 
//
bool OSUFlow::GenPathLines(VECTOR3* seeds, list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints,
			   float* tarray)
{

  if (has_data == false) DeferredLoadData(); 

  if (seedPtr!=NULL) delete[] seedPtr; 

        nSeeds = num_seeds; 
	seedPtr = seeds; 
	seedTimeArray = tarray; 

	listSeedTraces.clear();

	// execute streamline
	vtCPathLine* pPathLine;

	pPathLine = new vtCPathLine(flowField);

	pPathLine->SetTimeDir(dir); 
	pPathLine->SetLowerUpperAngle(3.0, 15.0);
	pPathLine->setMaxPoints(maxPoints);
	pPathLine->setSeedPoints(seedPtr, nSeeds, tarray);
	pPathLine->SetInitStepSize(initialStepSize);
	pPathLine->SetMaxStepSize(maxStepSize);
	pPathLine->setIntegrationOrder(FOURTH);
	pPathLine->execute(listSeedTraces);
	// release resource
	delete pPathLine;
	return true;
}

////////////////////////////////////////////////////////////////////////////
//
//    Take an input list of seeds (VECTOR4) consisting of 
//    positions and times 
//
bool OSUFlow::GenPathLines(VECTOR4* seeds, 
			   list<vtListTimeSeedTrace*>& listSeedTraces, 
			   TIME_DIR  dir, 
			   int num_seeds, 
			   int maxPoints)
{

  if (has_data == false)
    DeferredLoadData(); 

  nSeeds = num_seeds; 
  if (seedPtr!=NULL) delete [] seedPtr; 
  seedPtr = new VECTOR3[nSeeds]; 
  if (seedTimeArray!=NULL) delete [] seedTimeArray; 
  seedTimeArray = new float[nSeeds]; 
  for (int i=0; i<nSeeds; i++) {
    seedPtr[i][0] = seeds[i][0]; 
    seedPtr[i][1] = seeds[i][1]; 
    seedPtr[i][2] = seeds[i][2]; 
    seedTimeArray[i] = seeds[i][3]; 
  }

  listSeedTraces.clear();

  // execute streamline
  vtCPathLine* pPathLine;

  pPathLine = new vtCPathLine(flowField);
  pPathLine->SetTimeDir(dir); 
  pPathLine->SetLowerUpperAngle(3.0, 15.0);
  pPathLine->setMaxPoints(maxPoints);
  pPathLine->setSeedPoints(seeds, nSeeds);
  pPathLine->SetInitStepSize(initialStepSize);
  pPathLine->SetMaxStepSize(maxStepSize);
  pPathLine->setIntegrationOrder(FOURTH);
  pPathLine->execute(listSeedTraces);
  // release resource
  delete pPathLine;
  return true;
}

///////////////////////////////////////////////////////////////
//
// Use preset streakline seeds, all starting from current_time 
//
bool OSUFlow::GenStreakLines(vtStreakTraces& streakTraces, TIME_DIR dir, 
			     float current_time)
{

  if (has_data == false) DeferredLoadData(); 

  if (seedPtr==NULL)  {
    nSeeds = numSeeds[0]*numSeeds[1]*numSeeds[2];
    seedPtr = new VECTOR3[nSeeds];
    SeedGenerator* pSeedGenerator = new SeedGenerator((const float*)minRakeExt, 
						      (const float*)maxRakeExt, 
						      (const size_t*)numSeeds);
    pSeedGenerator->GetSeeds(seedPtr, bUseRandomSeeds);
    delete pSeedGenerator;
  }

  pStreakLine = new vtCStreakLine(flowField); 

  float currentT = current_time; 
  
  pStreakLine->SetTimeDir(dir); 
  pStreakLine->SetLowerUpperAngle(3.0, 15.0);
  pStreakLine->setSeedPoints(seedPtr, nSeeds, currentT);
  pStreakLine->SetInitStepSize(initialStepSize);
  pStreakLine->SetMaxStepSize(maxStepSize);
  pStreakLine->setIntegrationOrder(FOURTH);
  pStreakLine->execute((void*) &currentT, streakTraces); 
  delete pStreakLine; 
  return true; 
} 


bool OSUFlow::GenStreakLines(VECTOR3* seeds, vtStreakTraces& streakTraces, TIME_DIR dir,
			     int num_seeds, float current_time)
{

  if (has_data == false) DeferredLoadData(); 

  nSeeds = num_seeds; 
  seedPtr = seeds; 

  pStreakLine = new vtCStreakLine(flowField); 

  //otherwise one sterakline has already been created before 
  float currentT = current_time; 
  
  pStreakLine->SetTimeDir(dir); 
  pStreakLine->SetLowerUpperAngle(3.0, 15.0);
  pStreakLine->setSeedPoints(seedPtr, nSeeds, currentT);
  pStreakLine->SetInitStepSize(initialStepSize);
  pStreakLine->SetMaxStepSize(maxStepSize);
  pStreakLine->setIntegrationOrder(FOURTH);

  pStreakLine->execute((void*) &currentT, streakTraces); 
  delete pStreakLine; 
  return true; 
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
// fname is dataset file name
// sMin, sMax are local subdomain min and max
// dim is the total size of the domain
// t_min, t_max are min and max time steps (ignored if bStatic == true)
//
// Tom Peterka, 11/24/08
//
void OSUFlow::ReadData(const char* fname, bool bStatic, 
      VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim, int t_min, int t_max, 
      MPI_Comm comm) {

  flowName = new char[255];
  strcpy(flowName, fname);

  bStaticFlow = bStatic;

  if(bStaticFlow)
    ReadStaticFlowField(sMin, sMax, dim, comm);
  else
    ReadTimeVaryingFlowField(sMin, sMax, dim, t_min, t_max, comm);

  has_data = true;

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
void OSUFlow::ReadStaticFlowField(VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim, MPI_Comm comm) {

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
  MPI_Comm_rank(comm, &rank);

  // open the file
  if (MPI_File_open(comm, flowName, MPI_MODE_RDONLY, MPI_INFO_NULL,
       &fd) != MPI_SUCCESS)
    Error("Error: ReadRawData() cannot open %s\n",flowName);

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // sMin and sMax are [x][y][z]
  for (i = 0; i < 3; i++) {
    size[2 - i] = dim[i];
    start[2 - i] = sMin[i];
    subsize[2 - i] = sMax[i] - sMin[i] + 1;
  }
  size[2] *= 3;
  start[2] *= 3;
  subsize[2] *= 3;

  // allocate data space
  pData = new float[subsize[0] * subsize[1] * subsize[2]];
  if (!pData)
    Error("Error: ReadStaticFlowField() unable to allocate data space\n");

  // debug
//   fprintf(stderr,"dim = %.0f %.0f %.0f sMin = %.0f %.0f %.0f sMax = %.0f %.0f %.0f size = %d %d %d start = %d %d %d subsize = %d %d %d\n",dim[0],dim[1],dim[2],sMin[0],sMin[1],sMin[2],sMax[0],sMax[1],sMax[2],size[2],size[1],size[0],start[2],start[1],start[0],subsize[2],subsize[1],subsize[0]);

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
    Error("Error: ReadStaticFlowField() error rank %d read %d bytes instead of %d floats from file\n",
	  rank, status.count, sizeof(float) * subsize[0] * subsize[1] * subsize[2]);

  //swap bytes
#ifdef BYTE_SWAP
  for (i = 0; i < subsize[0] * subsize[1] * subsize[2]; i++)
    swap4((char *)&(pData[i]));
#endif

  MPI_File_close(&fd);

  // create the field
  CreateStaticFlowField(pData, sMin, sMax); 

  // debug
//   if (rank == 0) {
//     for (i = 0; i < subsize[0] / 3.0f * subsize[1] * subsize[2]; i += 100)
//       fprintf(stderr,"%.3f %.3f %.3f\n",pData[3 * i + 0],pData[3 * i + 1],pData[3 * i + 2]);
//   }

}
//---------------------------------------------------------------------------
//
// ReadTimeVaryingFlowField
//
// Read data collectively and create a vectorfield object based on my subdomain
// used MPI-IO collectives to perform the I/O
// reads a single subdomain
//
// sMin, sMax: corners of subdomain
// dim: size of total domain
// t_min, t_max: time range of subdomain
//
// Tom Peterka, 11/24/08
//
void OSUFlow::ReadTimeVaryingFlowField(VECTOR3 sMin, VECTOR3 sMax, 
	VECTOR3 dim, int t_min, int t_max, MPI_Comm comm) {

  float *Data; // the data
  MPI_File fd;
  MPI_Datatype filetype;
  int size[3]; // sizes of entire dataset
  int subsize[3]; // sizes of my subdomain
  int start[3]; // starting indices of my subdomain
  MPI_Status status;
  int err; // error status
  int rank;
  FILE *meta; // metadata file with file names of time step files
  int timesteps; // number of timesteps in the entire dataset
  char filename[256]; // individual timestep file name
  VECTOR3 **AllPt; // points in my time-space domain
  int npt; // number of points in my spatial domain
  int nt; // number of time steps in my time-space domain
  Solution* Soln;
  RegularCartesianGrid* RegCGrid;
  int i, t;

  // read metadata
  meta = fopen(flowName, "r");
  assert(meta != NULL);
  fscanf(meta, "%d", &timesteps);

  // init	
  MPI_Comm_rank(comm, &rank);
  AllPt = new VECTOR3 *[t_max - t_min + 1];
  assert(AllPt != NULL);
  gMin.Set(0.0,0.0,0.0); 
  gMax.Set((float)(dim[0] - 1), (float)(dim[1] - 1), (float)(dim[2] - 1));
  lMin = sMin;
  lMax = sMax;
  MinT = t_min;
  MaxT = t_max; 
  nt = t_max - t_min + 1;

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // sMin and sMax are [x][y][z]
  for (i = 0; i < 3; i++) {
    size[2 - i] = dim[i];
    start[2 - i] = sMin[i];
    subsize[2 - i] = sMax[i] - sMin[i] + 1;
  }
  npt = subsize[0] * subsize[1] * subsize[2];
  size[2] *= 3;
  start[2] *= 3;
  subsize[2] *= 3;

  // allocate data space
  Data = new float[3 * npt];
  assert(Data != NULL);

  // for all time steps
  for (t = 0; t < timesteps; t++) {

    fscanf(meta, "%s", filename);
    if (t < t_min || t > t_max)
      continue;

    // open the file
//     if (rank == 0)
//       fprintf(stderr,"%s\n",filename);
    err = MPI_File_open(comm, filename, MPI_MODE_RDONLY,
			MPI_INFO_NULL, &fd);
    assert(err == MPI_SUCCESS);

    // do the actual collective io
    MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
			     MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(fd, 0, MPI_FLOAT, filetype, (char *)"native", 
		      MPI_INFO_NULL);
    err = MPI_File_read_all(fd, Data,
	 subsize[0] * subsize[1] * subsize[2], MPI_FLOAT, &status);
    assert(err == MPI_SUCCESS);
    assert(status.count == sizeof(float) * 3 * npt);

    //swap bytes
#ifdef BYTE_SWAP
    for (i = 0; i < 3 * npt; i++)
      swap4((char *)&(Data[i]));
#endif

    MPI_File_close(&fd);

    // copy data from floats to VEC3s
    AllPt[t - t_min] = new VECTOR3[npt];
    for (i = 0; i < npt; i++)
      AllPt[t - t_min][i].Set(Data[i * 3 + 0], Data[i * 3 + 1], 
           Data[i * 3 + 2]);

  } // for all time steps

  // create the field
  Soln = new Solution(AllPt, npt, nt, t_min, t_max);
  RegCGrid = new RegularCartesianGrid(lMax[0] - lMin[0] + 1,   
	 lMax[1] - lMin[1] + 1, lMax[2] - lMin[2] + 1); 
  assert(Soln != NULL && RegCGrid != NULL);
  RegCGrid->SetBoundary(lMin, lMax);
  flowField = new CVectorField(RegCGrid, Soln, nt, MinT);

  // clean up
  fclose(meta);
  MPI_Type_free(&filetype);
  for (i=0; i < nt; i++)
    delete[] AllPt[i]; 
  delete[] AllPt; 
  delete[] Data;

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
