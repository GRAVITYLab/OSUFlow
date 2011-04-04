//	File:		OSUFlow.h
//
//	Author:		Liya Li and Han-Wei Shen 
//
//	Date:		Sept 2005
//
//	Description:	Definition of OSUFlow class. It contains the interface
//					between gui and underlying flow functions.
//
#ifndef _OSU_FLOW_H_
#define _OSU_FLOW_H_

#include "header.h"
#include "VectorMatrix.h"
#include "Field.h"
#include "Rake.h"
#include "FieldLine.h"

#include "FileReader.h"

//
// a few utilities that are not part of the OSUFlow class
//
void Error(const char *fmt, ...);
void swap4(char *n);

class OSUFlow
{
 public:
  // constructor and destructor
  OSUFlow();
  ~OSUFlow();

  ////////////////
  // data reader helper functions
  void LoadData(const char* fname, bool bStatic, bool deferred = false);
  void LoadData(const char* fname, bool bStatic, VECTOR3 pMin, 
		VECTOR3 pMax, bool deferred = false);
  void LoadData(const char* fname, bool bStatic, VECTOR3 pMin, 
		VECTOR3 pMax, int min_t, int max_t, bool deferred = false);
  void LoadData(const char* fname, bool bStatic, 
		float *from, float *to, float *size, int bt_max, 
		int t_min, int t_max);
  void LoadData(const char* fname, float *sMin, float *sMax, 
		float *dim, int min_t, int max_t, DataMode mode, float **data = NULL);
  void LoadData(char **dataset_files, int num_dataset_files, 
		float *sMin, float *sMax, float *dim, int min_t, int max_t, 
		DataMode mode, float **data = NULL);

  // create a flow field from input data array 
  CVectorField* CreateStaticFlowField(float*, int xsize, int ysize, int zsize, 
				      float* minb, float* maxB); 
  CVectorField* CreateTimeVaryingFlowField(float**, int xres, int yres, int zres, float* minB, 
					   float* maxB, int min_t, int max_t); 

  //////////////
  CVectorField* GetFlowField(void) { return flowField; }
  void  SetFlowField(CVectorField* field) { flowField = field; }

  VECTOR3 *GetSeeds(int& num) {num = numSeeds[0]*numSeeds[1]*numSeeds[2]; 
    return seedPtr;}

  // --- streamline methods 
  bool GenStreamLines(list<vtListSeedTrace*>&, TRACE_DIR, int, unsigned int);
  bool GenStreamLines(VECTOR3*, TRACE_DIR,const int num_seeds, const int maxPoints, 
		      list<vtListSeedTrace*>&, int64_t *seedIds = NULL, list<int64_t> *listSeedIds = NULL);
  // ---  pathline methods 
  // use preset seedPtr, all seeds start at currentT 
  bool GenPathLines(VECTOR4* seeds, list<vtListTimeSeedTrace*>& listSeedTraces, TIME_DIR, 
		    int num_seeds, int maxPoints, 
				int64_t *seedIds = NULL, list<int64_t> *listSeedIds = NULL); 
  bool GenPathLines(list<vtListTimeSeedTrace*>& listSeedTraces, TIME_DIR, 
		    int maxPoints, float currentT = 0.0);
  // use the input seed list, all seeds start at currentT 
  bool GenPathLines(VECTOR3 * seeds, list<vtListTimeSeedTrace*>& listSeedTraces, TIME_DIR, 
		    int num_seeds, int maxPoints, float currentT = 0.0); 
  // use the input seed list, all seeds start at different time in tarray 
  bool GenPathLines(VECTOR3 * seeds, list<vtListTimeSeedTrace*>& listSeedTraces, TIME_DIR, 
		    int num_seeds, int maxPoints, float* tarray); 

  // ---  streakline methods 
  // use preset seedPtr, all seeds start at current_time 
  bool GenStreakLines(vtStreakTraces& StreakTraces, TIME_DIR, 
		      float current_time); 
  // use the input seed list, all start at current_time 
  bool GenStreakLines(VECTOR3* seeds, vtStreakTraces& StreakTraces, TIME_DIR, 
		      int num_seeds, float current_time); 
  //	bool GenStreakLines(vtStreakTraces& StreakTraces, TIME_DIR, 
  //			    float* tarray, bool is_existing); 
  //	bool GenStreakLines(VECTOR3* seeds, vtStreakTraces& StreakTraces, TIME_DIR, 
  //			    int num_seeds, float * tarray, bool is_existing);  

  void Boundary(VECTOR3& minB, VECTOR3& maxB) { flowField->Boundary(minB, maxB); };
  void SetBoundary(VECTOR3 minB, VECTOR3 maxB) {flowField->SetBoundary(minB, maxB);}; 

  void SetRandomSeedPoints(const float min[3], const float max[3], int numSeeds);
  void SetRegularSeedPoints(const float min[3], const float max[3], const size_t numSeeds[3]);
  void SetSeedPoints(VECTOR3* seeds, int num_seeds);
  void SetIntegrationParams(float initStepSize, float maxStepSize);

  void NormalizeField(bool bLocal) {flowField->NormalizeField(bLocal);}
  void ScaleField(float scaleF) {flowField->ScaleField(scaleF); }
  int  NumTimeSteps() {return numTimesteps; } 
  void GetMinMaxTime(int& min_t, int& max_t) {min_t = MinT; max_t = MaxT; }
  void GetGlobalBounds(VECTOR3 &minB, VECTOR3 &maxB) {minB = gMin; maxB = gMax;}

  void DeleteData(void);

	//used by curvilinear code  
  void UserPickedSeedPoints(VECTOR3* ptr, int num);//added by lijie
  void LoadDataCurvilinear(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax)  ;  //added by lijie
  void LoadDataIrregular(const char* fname, bool bStatic, 
		       VECTOR3 sMin, VECTOR3 sMax) ;  //added by lijie


 protected:

 private:
	//used by curvilinear code  
  void InitStaticCurvilinearFlowField(VECTOR3 sMin, VECTOR3 sMax);//added by lijie, for curvilinear grid
  void InitTimeVaryingCurvilinearFlowField(void); //added by lijie, for curvilinear grid
  void InitStaticFlowField(Solution* pSolution, Grid* pGrid, VECTOR3 minB, 
				  VECTOR3 maxB);//added by lijie
  void InitStaticIrregularFlowField(VECTOR3 sMin, VECTOR3 sMax);//added by lijie

  void InitStaticFlowField(void);
  void InitStaticFlowField(VECTOR3 minb, VECTOR3 maxB); 
  void CreateStaticFlowField(float*, VECTOR3 minb, VECTOR3 maxB); 

  void InitTimeVaryingFlowField(void); 
  void InitTimeVaryingFlowField(int min_t, int max_t);
  void InitTimeVaryingFlowField(VECTOR3 minB, VECTOR3 maxB); 
  void InitTimeVaryingFlowField(VECTOR3 minb, VECTOR3 maxB, int min_t, int max_t);
  void InitTimeVaryingFlowField(VECTOR3 sMin, VECTOR3 sMax, 
				VECTOR3 dim, int bt_max, int t_min,
				int t_max);
  void InitFlowField(float *sMin, float *sMax, 
		     float *dim, int t_min, int t_max,
		     DataMode dm, float **data = NULL);
  void InitStaticFlowField(VECTOR3 sMin, VECTOR3 sMax, 
			   VECTOR3 dim);


  bool DeferredLoadData(); 

  float minRakeExt[3]; // minimal rake range 
  float maxRakeExt[3]; // maximal rake range
  unsigned int numSeeds[3]; // number of seeds
  int nSeeds; 
  VECTOR3 *seedPtr; 
  float *seedTimeArray; // the time associated with each seed
  bool bUseRandomSeeds;	// whether use randomly or regularly generated seeds
  CVectorField* flowField; // flow field data
  char* flowName; // name of file including information about field
  char **dataset_files; // individual file names of timesteps
  int num_dataset_files; // number of timestep files
  float initialStepSize; // for integration
  float maxStepSize;
  VECTOR3 gMin, gMax; // global min/max range 
  VECTOR3 lMin, lMax; // local min/max range
  int MinT, MaxT; //local time range 
  bool bStaticFlow; // static flow
  bool has_data; 
  int deferred_load_case; 
  int numTimesteps; 
  vtCStreakLine *pStreakLine; 

};

#endif
