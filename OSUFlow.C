#include "OSUFlow.h"

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
  if (has_data)
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

////////////////////////////////////////////////////////////////////////
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

/////////////////////////////////////////////////////////////////
//
//Read the whole datafile and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(void)
{
  int dimension[3]; 
  float* pData = NULL;
	
  pData = ReadStaticDataRaw(flowName, dimension); 
	
  // update the domain bounds 
  gMin.Set(0.0, 0.0, 0.0);
  gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
  
  lMin.Set(0.0, 0.0, 0.0);
  lMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));

  float minB[3], maxB[3]; 
	
  minB[0] = minB[1] = minB[2] = 0; 
  maxB[0] = dimension[0]-1; 
  maxB[1] = dimension[1]-1; 
  maxB[2] = dimension[2]-1; 

  flowField = CreateStaticFlowField(pData, dimension[0], dimension[1], 
				    dimension[2],  minB, maxB); 
  
  if(!flowField->IsNormalized())
    flowField->NormalizeField(true);
}

/////////////////////////////////////////////////////////////////
//
// Read only a subdomain and create a vectorfield object based on that 
//
void OSUFlow::InitStaticFlowField(VECTOR3 sMin, VECTOR3 sMax)
{
	int dimension[3], totalNum;
	float* pData = NULL;
	int lxdim, lydim, lzdim; 
	float minB[3], maxB[3]; 

	minB[0] = sMin[0]; minB[1] = sMin[1]; minB[2] = sMin[2];
	maxB[0] = sMax[0]; maxB[1] = sMax[1]; maxB[2] = sMax[2];

	pData = ReadStaticDataRaw(flowName, dimension, minB, maxB); 

	gMin.Set(0.0,0.0,0.0); 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
	lMin = sMin; lMax = sMax; //local data min/max range

	dimension[0] = maxB[0]-minB[0]+1; 
	dimension[1] = maxB[1]-minB[1]+1; 
	dimension[2] = maxB[2]-minB[2]+1; 
	
	flowField = CreateStaticFlowField(pData, dimension[0], dimension[1], dimension[2], 
			      minB,maxB); 

}


///////////////////////////////////////////////////////////////////////
//
// read the whole time sequence and create a time-varying vector field 
void OSUFlow:: InitTimeVaryingFlowField(void)
{
	int n_timesteps;
	int dimension[3]; 
	float** ppData = NULL;
	float minB[3], maxB[3]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
        lMin = gMin; lMax = gMax; 

	minB[0]=minB[1]=minB[2] = 0; 

	maxB[0]= dimension[0]-1; 
	maxB[1]= dimension[1]-1; 
	maxB[2]= dimension[2]-1; 
       
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], minB, maxB, 
					       0, n_timesteps-1);  
}

/////////////////////////////////////////////////////////////////
//
// read the whole time sequence within the subdomain 
//
void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB)
{
  	int n_timesteps;
	int dimension[3]; 
	float** ppData = NULL;
	float min_B[3], max_B[3]; 

	min_B[0]= minB[0]; 
	min_B[1]= minB[1]; 
	min_B[2]= minB[2]; 

	max_B[0]= maxB[0]; 
	max_B[1]= maxB[1]; 
	max_B[2]= maxB[2]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension, 
				     min_B, max_B, -1, -1); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
        lMin = minB; lMax = maxB; 

	dimension[0] = maxB[0]-minB[0]+1; 
	dimension[1] = maxB[1]-minB[1]+1; 
	dimension[2] = maxB[2]-minB[2]+1; 
	
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], min_B, max_B, 
					       0, n_timesteps-1);  
}


void OSUFlow:: InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB, int min_t, int max_t)
{
  	int n_timesteps;
	int dimension[3]; 
	float** ppData = NULL;
	float min_B[3], max_B[3]; 

	min_B[0]= minB[0]; 
	min_B[1]= minB[1]; 
	min_B[2]= minB[2]; 

	max_B[0]= maxB[0]; 
	max_B[1]= maxB[1]; 
	max_B[2]= maxB[2]; 

	ppData = ReadTimeVaryingDataRaw(flowName, n_timesteps, dimension, 
				     min_B, max_B, min_t, max_t); 

	// update the global bounds of the field. Assuming the same for all time steps 
	gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dimensions 
	gMax.Set((float)(dimension[0]-1), (float)(dimension[1]-1), (float)(dimension[2]-1));
        lMin = minB; lMax = maxB; 

	dimension[0] = maxB[0]-minB[0]+1; 
	dimension[1] = maxB[1]-minB[1]+1; 
	dimension[2] = maxB[2]-minB[2]+1; 
	
	flowField = CreateTimeVaryingFlowField(ppData, dimension[0], dimension[1], 
					       dimension[2], min_B, max_B, 
					       min_t, max_t);  

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

///////////////////////////////////////////////////////////////////

CVectorField* OSUFlow::CreateStaticFlowField(float *pData, 
				    int xdim, int ydim, int zdim, 
				    float* minB, float* maxB) 
{
  CVectorField* field; 
  Solution* pSolution;
  RegularCartesianGrid* pRegularCGrid;
  VECTOR3* pVector;
  VECTOR3** ppVector;
  VECTOR3 min_b, max_b; 

  int totalNum = xdim*ydim*zdim; 
  pVector = new VECTOR3[totalNum]; 
  ppVector = new VECTOR3*[1]; 

  for(int i=0; i<totalNum; i++) {
    pVector[i].Set(pData[i*3], pData[i*3+1], pData[i*3+2]); 
  }
  delete [] pData; 
  ppVector[0] = pVector; 

  pSolution = new Solution(ppVector, totalNum, 1);
  pRegularCGrid = new RegularCartesianGrid(xdim, ydim, zdim); 

  min_b[0] = minB[0]; min_b[1] = minB[1]; min_b[2] = minB[2]; 
  max_b[0] = maxB[0]; max_b[1] = maxB[1]; max_b[2] = maxB[2]; 

  pRegularCGrid->SetBoundary(min_b, max_b);

  assert(pSolution != NULL && pRegularCGrid != NULL);
  
  field = new CVectorField(pRegularCGrid, pSolution, 1);

  delete [] ppVector; 

  flowField = field; 
  return(field); 
}


//////////////////////////////////////////////////////////////////
//
//   ppData are assumed to contain the vector data of 
//   max_t-min_t+1 time steps 
//
//   xdim, ydim, zdim are the resolutions of the data block
//   minB and maxB specify the actual physical bounds of the block 
//   min_t and max_t are the time interval that this block represents 
//
CVectorField* OSUFlow::CreateTimeVaryingFlowField(float** ppData, 
						  int xdim, int ydim, int zdim, 
						  float* minB, float* maxB, 
						  int min_t, int max_t) 
{
  CVectorField* field; 
  Solution* pSolution;
  RegularCartesianGrid* pRegularCGrid;
  VECTOR3* pVector;
  VECTOR3** ppVector;
  VECTOR3 min_b, max_b; 
  
  int totalNum = xdim*ydim*zdim; 

  numTimesteps = max_t-min_t;   //??????? potential problem 
  ppVector = new VECTOR3 *[numTimesteps]; 

  for (int i=0; i<numTimesteps; i++) {
    pVector = new VECTOR3[totalNum]; 
    for (int j = 0; j<totalNum; j++)
      pVector[j].Set(ppData[i][j*3], ppData[i][j*3+1], ppData[i][j*3+2]); 
    delete[] ppData[i]; 
    ppVector[i] = pVector; 
  }
  delete[] ppData;

  min_b[0] = minB[0]; min_b[1] = minB[1]; min_b[2] = minB[2]; 
  max_b[0] = maxB[0]; max_b[1] = maxB[1]; max_b[2] = maxB[2]; 

  // create the flow field now 
  pSolution = new Solution(ppVector, totalNum, numTimesteps, min_t, max_t);
  pRegularCGrid = new RegularCartesianGrid(xdim, ydim, zdim); 			

  pRegularCGrid->SetBoundary(min_b, max_b);
  assert(pSolution != NULL && pRegularCGrid != NULL);
  field = new CVectorField(pRegularCGrid, pSolution, numTimesteps, min_t);

  for (int i=0; i<numTimesteps; i++) {
    delete [] ppVector[i]; 
  }
  delete[] ppVector; 

  flowField = field; 
  return(field); 
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
  vfprintf(stderr, fmt, argp);
  sleep(5);
  exit(0);

}
//-----------------------------------------------------------------------
//
// LoadData()
//
// loads the dataset
//
// fname is dataset file name
// sMin, sMax are local subdomain min and max
// dim is the total size of the domain
// bt_max is the max number of time steps in any block
// t_min, t_max are min and max time steps (ignored if bStatic == true)
//
void OSUFlow::LoadData(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim, int bt_max, 
		       int min_t, int max_t) {

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
	  
  if (bStaticFlow) {  // ignore the time range 
    numTimesteps = 1; 
    MinT = MaxT = 0; 
    InitStaticFlowField(sMin, sMax, dim);
  }

  else
    InitTimeVaryingFlowField(sMin, sMax, dim, bt_max, min_t, max_t); 

  has_data = true; 

}
//-----------------------------------------------------------------------
//
// InitStaticFlowField
//
// initialize a vectorfield object based on my subdomain
//
// sMin, sMax: corners of subdomain
// dim: size of total domain
//
void OSUFlow::InitStaticFlowField(VECTOR3 sMin, VECTOR3 sMax, VECTOR3 dim) {

  int totalNum;
  float* pData = NULL;
  int lxdim, lydim, lzdim; 
  float minB[3], maxB[3]; 
  int dimension[3]; // domain size
  int sub[3]; // subdomain size

  dimension[0] = dim[0];
  dimension[1] = dim[1];
  dimension[2] = dim[2];

  minB[0] = sMin[0]; minB[1] = sMin[1]; minB[2] = sMin[2];
  maxB[0] = sMax[0]; maxB[1] = sMax[1]; maxB[2] = sMax[2];

  pData = ReadStaticDataRaw(flowName, minB, maxB, dimension); 

  gMin.Set(0.0,0.0,0.0); 
  gMax.Set(dimension[0] - 1.0f, dimension[1] - 1.0f, dimension[2] - 1.0f);
  lMin = sMin; lMax = sMax; //local data min/max range

  sub[0] = maxB[0] - minB[0]+1; 
  sub[1] = maxB[1] - minB[1]+1; 
  sub[2] = maxB[2] - minB[2]+1; 
	
  flowField = CreateStaticFlowField(pData, sub[0], sub[1], sub[2], minB, maxB); 

}
//--------------------------------------------------------------------------
//
// InitTimeVaryingFlowField
//
// initialize a time-varying vectorfield object based on my subdomain
//
// sMin, sMax: corners of subdomain
// dim: size of total domain
// bt_max: max number of time steps in any block
// t_min, t_max: time range of subdomain
//
void OSUFlow::InitTimeVaryingFlowField(VECTOR3 sMin, VECTOR3 sMax, 
				       VECTOR3 dim, int bt_max, int t_min,
				       int t_max) {

  int n_timesteps;
  float** ppData = NULL;
  float minB[3], maxB[3]; 
  int dimension[3]; // domain size
  int sub[3]; // subdomain size

  dimension[0] = dim[0];
  dimension[1] = dim[1];
  dimension[2] = dim[2];

  minB[0] = sMin[0]; 
  minB[1] = sMin[1]; 
  minB[2] = sMin[2]; 

  maxB[0] = sMax[0]; 
  maxB[1] = sMax[1]; 
  maxB[2] = sMax[2]; 

  ppData = ReadTimeVaryingDataRaw(flowName, minB, maxB, dimension, 
				  bt_max, t_min, t_max); 

  // update the global bounds of the field
  // Assuming the same for all time steps 
  gMin.Set(0.0,0.0,0.0);   // assume all time steps have the same dims
  gMax.Set(dimension[0] - 1.0f, dimension[1] - 1.0f, dimension[2] - 1.0f);
  lMin = sMin; lMax = sMax; 
  MinT = t_min;
  MaxT = t_max; 

  sub[0] = maxB[0] - minB[0]+1; 
  sub[1] = maxB[1] - minB[1]+1; 
  sub[2] = maxB[2] - minB[2]+1; 
	
  flowField = CreateTimeVaryingFlowField(ppData, sub[0], sub[1], sub[2], 
					 minB, maxB, t_min, t_max);  

}
//--------------------------------------------------------------------------

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
