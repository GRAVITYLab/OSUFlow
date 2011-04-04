/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Vector Field Data
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _VECTOR_FIELD_LINE_H_
#define _VECTOR_FIELD_LINE_H_

#include "Field.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

//////////////////////////////////////////////////////////////////////////
// definition
//////////////////////////////////////////////////////////////////////////
#define MAX_LENGTH 1000000
const float STREAM_ACCURACY = EPS;
#define INFINITE_LINE -1				// advect fieldlines as far as possible

enum INTEG_ORD{ SECOND = 2, FOURTH = 4};		// integration order
enum TIME_DIR{ BACKWARD = -1, FORWARD = 1};		// advection direction
enum TIME_DEP{ STEADY=0,UNSTEADY=1 };	
enum TRACE_DIR{OFF=0, BACKWARD_DIR=1, FORWARD_DIR=2, BACKWARD_AND_FORWARD=3};
enum ADVECT_STATUS{NONE = -2, OUT_OF_BOUND = -1, CRITICAL_POINT = 0, OKAY = 1};

//////////////////////////////////////////////////////////////////////////
// information about particles
//////////////////////////////////////////////////////////////////////////
class vtParticleInfo
{
public:
	PointInfo m_pointInfo;		// basic information about this particle
	float m_fStartTime;			// start time
	int itsValidFlag;			// whether this particle is valid or not
	int itsNumStepsAlive;		// number of steps alive
	int ptId;					// particle ID

public:
	vtParticleInfo(void)
	{
	}

	vtParticleInfo(vtParticleInfo* source)
	{
		m_pointInfo.Set(source->m_pointInfo);
		m_fStartTime = source->m_fStartTime;
		itsValidFlag = source->itsValidFlag;
		itsNumStepsAlive = 0;
		ptId = -1;
	}

	vtParticleInfo(vtParticleInfo& source)
	{
		m_pointInfo.Set(source.m_pointInfo);
		m_fStartTime = source.m_fStartTime;
		itsValidFlag = source.itsValidFlag;
		ptId = source.ptId;
		itsNumStepsAlive = 0;
	}
};

typedef list<vtParticleInfo*> vtListParticle;
typedef list<vtParticleInfo*>::iterator vtListParticleIter;
typedef list<VECTOR3*> vtListSeedTrace;      // positions
typedef list<VECTOR4*> vtListTimeSeedTrace;  // positions and times 

//////////////////////////////////////////////////////////////////////////
// base class for all field lines, and the structure is like:
//                                    vtCFieldLine
//	                                  /  \
//			      vtCStreamLine   vtCTimeVaryingFieldLine
//               	                 	|     /   \
//                                    vtCPathline  vtCTimeLine   vtCStreakLine
//////////////////////////////////////////////////////////////////////////
class vtCFieldLine
{
protected:
	int m_nNumSeeds;		// number of seeds
	INTEG_ORD m_integrationOrder;	// integration order
	TIME_DIR m_timeDir;	        // advection direction
	TIME_DEP m_itsTimeDep;
	float m_fInitTime;
	float m_fStepTime;
	float m_fInitStepSize;	 // initial advection step size of particle
	float m_fDurationTime;
	float m_fLowerAngleAccuracy;	// for adaptive stepsize 
	float m_fUpperAngleAccuracy;
	float m_fMaxStepSize;	        // maximal advection stepsize
	int m_nMaxsize;		// maximal number of particles this line advects
	vtListParticle m_lSeeds;	// list of seeds
	list<int64_t> m_lSeedIds;	// list of seed ids
	CVectorField* m_pField;	        // vector field
	float m_fStationaryCutoff;	// cutoff value for critical points

public:
	vtCFieldLine(CVectorField* pField);
	virtual ~vtCFieldLine(void);

	void setSeedPoints(VECTOR3* points, int numPoints, float t, int64_t *seedIds = NULL); 
	void setSeedPoints(VECTOR3* points, int numPoints, float* tarray);
	void setSeedPoints(VECTOR4* points, int numPoints, int64_t *seedIds = NULL);
	void setMaxPoints(int val) { m_nMaxsize = val; }
	void setIntegrationOrder(INTEG_ORD ord) { m_integrationOrder = ord; }
	int  getMaxPoints(void){ return m_nMaxsize; }
	INTEG_ORD getIntegrationOrder(void){ return m_integrationOrder; }
	void SetMaxStepSize(float stepsize) {m_fMaxStepSize = stepsize;}
	float GetMaxStepSize(void) {return m_fMaxStepSize;}
	void SetInitStepSize(float initStep) { m_fInitStepSize = initStep; }
	float GetInitStepSize(void) { return m_fInitStepSize; }
	void SetLowerUpperAngle(float lowerAngle, float upperAngle) {m_fLowerAngleAccuracy = lowerAngle; m_fUpperAngleAccuracy = upperAngle;}
	void SetStationaryCutoff(float cutoff) {m_fStationaryCutoff = cutoff;}

protected:
	void releaseSeedMemory(void);
	int euler_cauchy(TIME_DIR, TIME_DEP,float*, float);
	int runge_kutta4(TIME_DIR, TIME_DEP, PointInfo&, float*, float);
	int runge_kutta2(TIME_DIR, TIME_DEP, PointInfo&, float*, float);
	int adapt_step(const VECTOR3& p2, const VECTOR3& p1, const VECTOR3& p0, float dt_estimate,float* dt);
};

//////////////////////////////////////////////////////////////////////////
// class declaratin of timevaryingfieldline
//////////////////////////////////////////////////////////////////////////
class vtCTimeVaryingFieldLine : public vtCFieldLine
{
protected:
	int m_itsTimeAdaptionFlag;
	int m_itsMaxParticleLife;   // how long the particles be alive
	int m_itsMapWithTimeFlag;
	int m_itsNeedResetFlag;
	int m_itsWrapTimeFlag;
	float m_itsTimeInc;
	list<vtParticleInfo*> m_itsParticles;
public:
	vtCTimeVaryingFieldLine(CVectorField* pField);
	~vtCTimeVaryingFieldLine(void);

	void SetTimeDir(TIME_DIR dir) { m_timeDir = dir; }

	void SetInjectionTime(float init, float step, float duration, TIME_DIR dir = FORWARD)
	{
		m_fInitTime = init;
		m_fStepTime = step;
		m_fDurationTime = duration;
		m_timeDir = dir;
	}

	void SetTimeAdaptionMode(int onoff)
	{
		m_itsTimeAdaptionFlag = onoff;
	}

	virtual void setParticleLife(int steps);
	virtual int getParticleLife(void);
	virtual void setTimeMapping(int enabled);
	virtual int getTimeMapping(void);
	virtual void killAllParticles(void) { m_itsNeedResetFlag = 1; }

protected:
	// code shared by all time varying field lines here
	void releaseParticleMemory(void);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    float initialTime,
			    vtParticleInfo& finalPoint,
			    float finalTime);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    float initialTime,
			    float finalTime,
			    vtListSeedTrace& seedTrace);
	int advectParticle( INTEG_ORD int_order, 
			    vtParticleInfo& initialPoint,
			    float initialTime,
			    float finalTime,
			    vtListTimeSeedTrace& seedTrace);
};

//////////////////////////////////////////////////////////////////////////
// class declaration of pathline
//////////////////////////////////////////////////////////////////////////
typedef struct vtPathlineParticle
{
	VECTOR3 pos;
	int ptId;
        float time; 
}vtPathlineParticle;

class vtCPathLine : public vtCTimeVaryingFieldLine
{
public:
	vtCPathLine(CVectorField* pField);
	~vtCPathLine(void);

	void execute(list<vtListTimeSeedTrace*>& listSeedTraces, list<int64_t> *listSeedIds = NULL);
	void execute(list<vtPathlineParticle*>& listSeedTraces);
		
protected:
	// code specific to pathline
	void computePathLine(list<vtListTimeSeedTrace*>&, list<int64_t> *listSeedIds = NULL);
	void computePathLine(list<vtPathlineParticle*>&);
};

//////////////////////////////////////////////////////////////////////////
// class declaration of streakline
// inject particles from a point
//////////////////////////////////////////////////////////////////////////
typedef struct vtStreakParticle
{
	PointInfo itsPoint;			// basic information
	float itsTime;				// start time
	int traceId;				// which tracing path this particle belongs to, used to track advection
}vtStreakParticle;

typedef vector<vtStreakParticle*> vtListStreakParticle;	    // one pathline trace from a particle
typedef vector<vtListStreakParticle*> vtStreakTraces;	    // all pathline traces from all particles released
typedef vector<vtStreakParticle*>::iterator vtStreakParticleIter;
typedef vector<vtListStreakParticle*>::iterator vtStreakTracesIter; 

class vtCStreakLine : public vtCTimeVaryingFieldLine
{
public:
	vtCStreakLine(CVectorField* pField);
	~vtCStreakLine(void);
	void execute(const void* userData, vtStreakTraces& listSeedTraces);

protected:
	// code specific to streakline
	void computeStreakLine(const void*, vtStreakTraces&);
	void advectOldParticles(vtListParticleIter start, 
				vtListParticleIter end, 
				vtStreakTraces& listSeedTraces,
				float initialTime, float finalTime,
				vector<vtListParticleIter>& deadList);
private:
	int nHowManyTraces;
};

//////////////////////////////////////////////////////////////////////////
// class declaration of timeline
// release a set of particles, injected periodically
//////////////////////////////////////////////////////////////////////////
class vtCTimeLine : public vtCTimeVaryingFieldLine
{
public:
	vtCTimeLine(CVectorField* pField);
	~vtCTimeLine(void);

	void execute(const void* userData, vtListStreakParticle& listSeedTraces);
	void setTimeDelay(int delay);
	int getTimeDelay(void);

protected:
	// code specific to timeline
	void computeTimeLine(const void*, vtListStreakParticle&);
	void advectOldParticles(vtListParticleIter start, 
				vtListParticleIter end, 
				vtListStreakParticle& listSeedTraces,
				float initialTime, float finalTime,
				vector<vtListParticleIter>& deadList);
	int m_itsTimeDelay;
	int numTillRelease;
};

//////////////////////////////////////////////////////////////////////////
// class declaration of streamline
//////////////////////////////////////////////////////////////////////////
class vtCStreamLine : public vtCFieldLine
{
public:
	vtCStreamLine(CVectorField* pField);
	~vtCStreamLine(void);

	void execute(const void* userData, list<vtListSeedTrace*>& listSeedTraces,
				list<int64_t> *listSeedIds = NULL);
	void setForwardTracing(int enabled);
	void setBackwardTracing(int enabled);
	int  getForwardTracing(void);
	int  getBackwardTracing(void);
	int executeInfiniteAdvection(TIME_DIR, TIME_DEP, vtListSeedTrace&, float&, vector<float>*);
	int AdvectOneStep(TIME_DIR, INTEG_ORD, TIME_DEP, PointInfo&, VECTOR3&);

protected:
	void computeStreamLine(const void* userData, list<vtListSeedTrace*>& listSeedTraces, list<int64_t> *listSeedIds = NULL);
	int computeFieldLine(TIME_DIR, INTEG_ORD, TIME_DEP, vtListSeedTrace&, PointInfo&);

	TRACE_DIR m_itsTraceDir;
	float m_fPsuedoTime;
	float m_fCurrentTime;
};

#endif

