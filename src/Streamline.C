/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

FILE* fDebugOut;
//////////////////////////////////////////////////////////////////////////
// definition of class FieldLine
//////////////////////////////////////////////////////////////////////////
vtCStreamLine::vtCStreamLine(CVectorField* pField):
vtCFieldLine(pField),
m_itsTraceDir(BACKWARD_AND_FORWARD),
m_fPsuedoTime(0.0),
m_fCurrentTime(0.0)
{
#ifdef DEBUG
	fDebugOut = fopen(".\\DebugInfo\\debug.txt", "w");
#endif
}

vtCStreamLine::~vtCStreamLine(void)
{
#ifdef DEBUG
	fclose(fDebugOut);
#endif
}

//////////////////////////////////////////////////////////////////////////
// Set the advection direction
//////////////////////////////////////////////////////////////////////////
void vtCStreamLine::setForwardTracing(int enabled)
{
	switch( m_itsTraceDir )
	{
	case BACKWARD_DIR:
		if(enabled)
			m_itsTraceDir = BACKWARD_AND_FORWARD;
		break;
	case FORWARD_DIR:
		if(!enabled)
			m_itsTraceDir = OFF;
		break;
	case BACKWARD_AND_FORWARD:
		if(!enabled)
			m_itsTraceDir = BACKWARD_DIR;
		break;
	case OFF:
		if (enabled)
			m_itsTraceDir = FORWARD_DIR;
		break;
	}
}

void vtCStreamLine::setBackwardTracing(int enabled)
{
	switch( m_itsTraceDir )
	{
	  case BACKWARD_DIR:
		  if(!enabled)
			  m_itsTraceDir = OFF;
		  break;
	  case FORWARD_DIR:
		  if(enabled)
			  m_itsTraceDir = BACKWARD_AND_FORWARD;
		  break;
	  case BACKWARD_AND_FORWARD:
		  if(!enabled)
			  m_itsTraceDir = FORWARD_DIR;
		  break;
	  case OFF:
		  if (enabled)
			  m_itsTraceDir = BACKWARD_DIR;
		  break;
	}
}

int vtCStreamLine::getForwardTracing(void)
{
	if( m_itsTraceDir==FORWARD_DIR || m_itsTraceDir==BACKWARD_AND_FORWARD )
		return 1;
	else
		return 0;
}

int vtCStreamLine::getBackwardTracing(void)
{
	if( m_itsTraceDir==BACKWARD_DIR || m_itsTraceDir==BACKWARD_AND_FORWARD )
		return 1;
	else
		return 0;
}

//////////////////////////////////////////////////////////////////////////
// Compute streamlines.
// output
//	listSeedTraces: For each seed, return a list keeping the trace it
//					advects
//	listSeedIds: For each seed, a unique id may be associated with it. 
//					These unique ids will populate the list (defaults to NULL)
//////////////////////////////////////////////////////////////////////////
void vtCStreamLine::execute(const void* userData, 
			    list<vtListSeedTrace*>& listSeedTraces,
					list<int64_t> *listSeedIds)
{
	m_fCurrentTime = *(float *)userData;
	computeStreamLine(userData, listSeedTraces, listSeedIds);
}

void vtCStreamLine::computeStreamLine(const void* userData,
				      list<vtListSeedTrace*>& listSeedTraces,
							list<int64_t> *listSeedIds)
{

	vtListParticleIter sIter;
	list<int64_t>::iterator sIdIter;

	#if	0	// MOD-BY-LEETEN 02/07/2011-FROM:
	for(sIter = m_lSeeds.begin(), sIdIter = m_lSeedIds.begin(); 
			sIter != m_lSeeds.end(); ++sIter, ++sIdIter)
	#else	// MOD-BY-LEETEN 02/07/2011-TO:
	if( !m_lSeedIds.empty() ) 
			sIdIter = m_lSeedIds.begin();
	for(sIter = m_lSeeds.begin(); sIter != m_lSeeds.end(); ++sIter)
	#endif	// MOD-BY-LEETEN 02/07/2011-END
	{

		vtParticleInfo* thisSeed = *sIter;

		if(thisSeed->itsValidFlag == 1)			// valid seed
		{

			if(m_itsTraceDir & BACKWARD_DIR)
			{
				vtListSeedTrace* backTrace;
				backTrace = new vtListSeedTrace;
				computeFieldLine(BACKWARD,m_integrationOrder, STEADY, *backTrace, thisSeed->m_pointInfo);
				listSeedTraces.push_back(backTrace);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
			if(m_itsTraceDir & FORWARD_DIR)
			{
				vtListSeedTrace* forwardTrace;
				forwardTrace = new vtListSeedTrace;
				computeFieldLine(FORWARD,m_integrationOrder, STEADY, *forwardTrace, thisSeed->m_pointInfo);
				listSeedTraces.push_back(forwardTrace);
				if (listSeedIds != NULL)
					(*listSeedIds).push_back(*sIdIter);
			}
		}
		// ADD-BY-LEETEN 02/07/2011-BEGIN
		if( !m_lSeedIds.empty() ) 
			sIdIter++;
		// ADD-BY-LEETEN 02/07/2011-END
	}
}

int vtCStreamLine::computeFieldLine( TIME_DIR time_dir,
									 INTEG_ORD integ_ord,
									 TIME_DEP time_dep, 
									 vtListSeedTrace& seedTrace,
									 PointInfo& seedInfo)
{
	int count = 0, istat;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	float dt, dt_estimate, mag, curTime;
	VECTOR3 vel;
	float cell_volume;

	// the first particle
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if(istat == OUT_OF_BOUND)  {
		return OUT_OF_BOUND;
	}
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff)) {
		return CRITICAL_POINT;
	}

	thisParticle = seedInfo;

	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	curTime = m_fCurrentTime;
	count++;

	// get the initial stepsize
	// this is a bug I think ...
	cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
	mag = vel.GetMag();
// 		printf(" **** volume = %f mag = %f  ....", cell_volume, mag); 
	if(fabs(mag) < 1.0e-6f)
		dt_estimate = 1.0e-5f;
	else
		dt_estimate = pow(cell_volume, (float)0.3333333f) / mag;
		dt = m_fInitStepSize * dt_estimate;
		//	printf(" **** dt = %f  ....", dt); 

#ifdef DEBUG
	fprintf(fDebugOut, "****************new particle*****************\n");
	fprintf(fDebugOut, "seed: %f, %f, %f with step size %f\n", seedInfo.phyCoord[0], seedInfo.phyCoord[1], seedInfo.phyCoord[2], dt);
#endif

	// start to advect
	while(count < m_nMaxsize)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
		else
			istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND)			// out of boundary
		  {
		    //		    printf(" count = %d here out of bound %f %f %f \n", count, 
		    //			   thisParticle.phyCoord[0], thisParticle.phyCoord[1],
		    //	  		   thisParticle.phyCoord[2]); 
		    seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			return OUT_OF_BOUND;
		  }
		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, m_fCurrentTime, vel);
		if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		  {
		    //printf(" here critical point.\n");
			return CRITICAL_POINT;
		  }
		else
		{
		  //#ifdef DEBUG
		  //			fprintf(fDebugOut, "point: %f, %f, %f with step size %f\n", thisParticle.phyCoord[0], thisParticle.phyCoord[1], thisParticle.phyCoord[2], dt);
		  //		  printf("point: %f, %f, %f with step size %f\n", thisParticle.phyCoord[0], thisParticle.phyCoord[1], thisParticle.phyCoord[2], dt);
			//#endif
			seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			count++;
		}

		if(count > 2)
		{
			float tempt;
			tempt = dt;
			adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
		}
	}

	return OKAY;
}

// streamline advects as far as possible till the boundary or terminates at critical points
// only two cases will happen: OUT_OF_BOUNDARY & CRITICAL_POINT
int vtCStreamLine::executeInfiniteAdvection(TIME_DIR time_dir,
											TIME_DEP time_dep,
											vtListSeedTrace& seedTrace,
											float& totalStepsize,
											vector<float>* vStepsize)
{
	int istat;
	vtParticleInfo* thisSeed;
	PointInfo thisParticle, prevParticle, second_prevParticle, seedInfo;
	vtListParticleIter sIter;
	float dt, dt_estimate, cell_volume, mag, curTime;
	VECTOR3 vel;
	INTEG_ORD integ_ord;

	// initialize
	integ_ord = m_integrationOrder;
	sIter = m_lSeeds.begin();
	thisSeed = *sIter;
	seedInfo = thisSeed->m_pointInfo;
	seedTrace.clear();
	totalStepsize = 0.0;

	// the first particle
	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;

	thisParticle = seedInfo;
	curTime = m_fCurrentTime;

	// get the initial stepsize
	cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
	mag = vel.GetMag();
	dt_estimate = pow(cell_volume, (float)0.3333333f) / mag;
	dt = m_fInitStepSize * dt_estimate;

#ifdef DEBUG
	fprintf(fDebugOut, "****************new particle*****************\n");
	fprintf(fDebugOut, "seed: %f, %f, %f with step size %f\n", seedInfo.phyCoord[0], seedInfo.phyCoord[1], seedInfo.phyCoord[2], dt);
#endif

	// start to advect
	while(true)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
		else
			istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND)			// out of boundary
		{
			// find the boundary intersection point
			VECTOR3 intersectP, startP, endP;
			float oldStepsize, stepSize; 
			oldStepsize = stepSize = dt;
			startP = prevParticle.phyCoord;
			endP = thisParticle.phyCoord;
			m_pField->BoundaryIntersection(intersectP, startP, endP, &stepSize, oldStepsize);
			totalStepsize += stepSize;
			seedTrace.push_back(new VECTOR3(intersectP));
			if(vStepsize != NULL)
				vStepsize->push_back(stepSize);
			return OUT_OF_BOUND;
		}

		// find the loop
		list<VECTOR3*>::iterator searchIter;
		searchIter = seedTrace.end();
		searchIter--;
		for(; searchIter != seedTrace.begin(); searchIter--)
		{
			if(thisParticle.phyCoord == **searchIter)			// loop
				return CRITICAL_POINT;
		}

		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, m_fCurrentTime, vel);
		seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
		if(vStepsize != NULL)
			vStepsize->push_back(dt);
		totalStepsize += dt;

#ifdef DEBUG
		fprintf(fDebugOut, "****************advected particle*****************\n");
		fprintf(fDebugOut, "pos = (%f, %f, %f), vel = (%f, %f, %f) with step size %f, total step size %f\n", thisParticle.phyCoord[0], thisParticle.phyCoord[1], thisParticle.phyCoord[2], vel[0], vel[1], vel[2], dt, totalStepsize);
#endif

		if(totalStepsize > 4000.0)
		{
			printf("curl\n");
			vel.Set(0.0, 0.0, 0.0);
		}

		if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
			return CRITICAL_POINT;
		if((int)seedTrace.size() > 2)
			adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
	}

	return OUT_OF_BOUND;
}

int vtCStreamLine::AdvectOneStep(TIME_DIR time_dir,
								 INTEG_ORD integ_ord,
								 TIME_DEP time_dep,
								 PointInfo& seedInfo,
								 VECTOR3& finalP)
{
	int istat;
	PointInfo thisParticle;
	VECTOR3 vel;
	float curTime = m_fCurrentTime;
	float dt = m_fInitStepSize;

	finalP.Set(-1, -1, -1);
	thisParticle = seedInfo;
	
	// the first particle
	istat = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, m_fCurrentTime, vel);
	if(istat == OUT_OF_BOUND)
		return OUT_OF_BOUND;
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;
	if(integ_ord == SECOND)
		istat = runge_kutta2(time_dir, time_dep, thisParticle, &curTime, dt);
	else
		istat = runge_kutta4(time_dir, time_dep, thisParticle, &curTime, dt);

	if(istat == OUT_OF_BOUND)			// out of boundary
			return OUT_OF_BOUND;
	m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, m_fCurrentTime, vel);
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
		return CRITICAL_POINT;
	else
		finalP = thisParticle.phyCoord;
	return OKAY;
}
