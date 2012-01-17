//*************************************************************************
//			  OSU Flow Vector Field						
//                        Created:  Han-Wei Shen, Liya Li
//			  The Ohio State University 
//                        Date:06/2005									   
//                        TimeVaryingFieldLines							
//*************************************************************************

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// methods common to all time varying field lines
//////////////////////////////////////////////////////////////////////////
vtCTimeVaryingFieldLine::vtCTimeVaryingFieldLine(CVectorField* pField) : 
vtCFieldLine(pField),
m_itsTimeAdaptionFlag(1),
m_itsMaxParticleLife(200),
m_itsMapWithTimeFlag(0),
m_itsWrapTimeFlag(1),
m_itsTimeInc(1.0)
{
	m_timeDir = FORWARD;
}

vtCTimeVaryingFieldLine::~vtCTimeVaryingFieldLine(void)
{
	releaseParticleMemory();
}

///////////////////////////////////////////////

void vtCTimeVaryingFieldLine::setParticleLife(int steps)
{
	m_itsMaxParticleLife = steps;
}

int vtCTimeVaryingFieldLine::getParticleLife(void)
{
	return m_itsMaxParticleLife;
}

void vtCTimeVaryingFieldLine::setTimeMapping(int enabled)
{
	m_itsMapWithTimeFlag = enabled;
}

int vtCTimeVaryingFieldLine::getTimeMapping(void)
{
	return m_itsMapWithTimeFlag;
}

void vtCTimeVaryingFieldLine::releaseParticleMemory(void)
{
	vtListParticleIter pIter = m_itsParticles.begin();
	for( ; pIter != m_itsParticles.end(); ++pIter )
	{
		vtParticleInfo* thisPart = *pIter;
		delete thisPart;
	}
	m_itsParticles.erase(m_itsParticles.begin(), m_itsParticles.end());
}

//////////////////////////////////////////////////////////////////////////
// advect the particle from initialTime to finalTime
//////////////////////////////////////////////////////////////////////////
int vtCTimeVaryingFieldLine::advectParticle(INTEG_ORD int_order, 
					    vtParticleInfo& initialPoint,
					    float initialTime,
					    vtParticleInfo& finalPoint,
					    float finalTime)
{
	int istat;
	float curTime, dt;
	PointInfo pt;

	pt = initialPoint.m_pointInfo;

	if(m_itsTimeAdaptionFlag == 1)
		dt = (finalTime - initialTime) * 0.5;
	else
		dt = finalTime - initialTime;

	curTime = initialTime;

	while(curTime < finalTime)
	{
		float error;
		if(int_order == SECOND)
			istat = runge_kutta2(m_timeDir, UNSTEADY, pt, &curTime, dt);
		else if(int_order == FOURTH)
			istat = runge_kutta4(m_timeDir, UNSTEADY, pt, &curTime, dt);
		else if(int_order == RK45)
			istat = runge_kutta45(m_timeDir, UNSTEADY, pt, &curTime,dt,&error);
		else
			return OUT_OF_BOUND;

		if(istat != OKAY)
			return istat;
	}

	finalPoint.m_pointInfo = pt;
	return istat;
}

//////////////////////////////////////////////////////////////////////////
// advect the particle until it goes out of boundary or vel = 0
// return back the track of advection
//////////////////////////////////////////////////////////////////////////
int vtCTimeVaryingFieldLine::advectParticle(INTEG_ORD int_order, 
                                            vtParticleInfo& initialPoint,
                                            float initialTime,
                                            float finalTime,
                                            vtListSeedTrace& seedTrace)
{  
	int count = 0, istat, res;
	PointInfo seedInfo;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	float dt, dt_estimate, cell_volume, mag, curTime;
	VECTOR3 vel;

	// the first particle
	seedInfo = initialPoint.m_pointInfo;
	res = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, 
	                        initialTime, vel);
	if(res == OUT_OF_BOUND){
		return OUT_OF_BOUND;
	}
	if((fabs(vel[0]) < m_fStationaryCutoff) && 
	   (fabs(vel[1]) < m_fStationaryCutoff) && 
	   (fabs(vel[2]) < m_fStationaryCutoff)) {
	  return CRITICAL_POINT;
	}

	thisParticle = seedInfo;
	seedTrace.push_back(new VECTOR3(seedInfo.phyCoord));
	curTime = initialTime;
	count++;

	// get the initial stepsize
	switch(m_pField->GetCellType())
	{
	case CUBE:
		dt = dt_estimate = m_fInitialStepSize;
		break;

	case TETRAHEDRON:
		cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
		mag = vel.GetMag();
		if(fabs(mag) < 1.0e-6f)
			dt_estimate = 1.0e-5f;
		else
			dt_estimate = pow(cell_volume, (float)0.3333333f) / mag;
		dt = dt_estimate;
		break;

	default:
		break;
	}

	// start to advect
	while(count < m_nMaxsize)
	{
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;

		if(int_order == SECOND || int_order == FOURTH)
			istat = oneStepGeometric(int_order, m_timeDir, UNSTEADY,
			                         thisParticle, prevParticle, 
			                         second_prevParticle, &curTime, &dt, count);
		else if(int_order == RK45)
			istat = oneStepEmbedded(int_order, m_timeDir, UNSTEADY,
			                        thisParticle, &curTime, &dt);
		else
			return OUT_OF_BOUND;

		if(istat == OUT_OF_BOUND)  {			// out of boundary
			return OUT_OF_BOUND;
		}

		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, 
		                  thisParticle, curTime, vel);
		if((fabs(vel[0]) < m_fStationaryCutoff) && 
		   (fabs(vel[1]) < m_fStationaryCutoff) && 
		   (fabs(vel[2]) < m_fStationaryCutoff)) {
			return CRITICAL_POINT;			// arrives at a critical point
		}
		else
		{
			seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			count++;
		}

		if((curTime < m_pField->GetMinTimeStep()) || 
		   (curTime > m_pField->GetMaxTimeStep())) {
			return OUT_OF_BOUND;
		}
	}

	return OKAY;
}




int vtCTimeVaryingFieldLine::advectParticle(INTEG_ORD int_order, 
					    vtParticleInfo& initialPoint,
					    float initialTime,
					    float finalTime,
					    vtListTimeSeedTrace& seedTrace)
{  
	int count = 0, istat, res;
	PointInfo seedInfo;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	PointInfo third_prevParticle;
	float dt, dt_attempt, dt_estimate, cell_volume, mag, curTime, prevCurTime;
	VECTOR3 vel;

	// the first particle
	seedInfo = initialPoint.m_pointInfo;
	res = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo,
	                        initialTime, vel);
	if(res == OUT_OF_BOUND)  {
		return OUT_OF_BOUND;
	}
	if((fabs(vel[0]) < m_fStationaryCutoff) && 
	   (fabs(vel[1]) < m_fStationaryCutoff) &&
	   (fabs(vel[2]) < m_fStationaryCutoff)) {
		return CRITICAL_POINT;
	}

	thisParticle = seedInfo;
	VECTOR4 *p = new VECTOR4; 
	(*p)[0] = seedInfo.phyCoord[0]; 
	(*p)[1] = seedInfo.phyCoord[1]; 
	(*p)[2] = seedInfo.phyCoord[2]; 
	(*p)[3] = initialTime; 
	seedTrace.push_back(p); 
	curTime = initialTime;
	count++;

	// get the initial stepsize
	dt = m_fInitialStepSize;
	if(m_adaptStepSize)
	{
		switch(m_pField->GetCellType())
		{
		case CUBE:
			dt = dt_estimate = m_fInitialStepSize;
			break;

		case TETRAHEDRON:
			cell_volume = m_pField->volume_of_cell(seedInfo.inCell);
			mag = vel.GetMag();
			if(fabs(mag) < 1.0e-6f)
				dt_estimate = 1.0e-5f;
			else
				dt_estimate = pow(cell_volume, (float)0.3333333f) / mag;
			dt = dt_estimate;
			break;

		default:
			break;
		}
	}

	// start to advect
	while(count < m_nMaxsize)
	{
		third_prevParticle = second_prevParticle;
		second_prevParticle = prevParticle;
		prevParticle = thisParticle;
		prevCurTime = curTime;
		dt_attempt = dt;

		if(int_order == SECOND || int_order == FOURTH)
			istat = oneStepGeometric(int_order, m_timeDir, UNSTEADY,
									 thisParticle, prevParticle,
			                         second_prevParticle, &curTime, &dt, count);
		else if(int_order == RK45)
			istat = oneStepEmbedded(int_order, m_timeDir, UNSTEADY,
			                        thisParticle, &curTime, &dt);
		
		// check if the step failed
		if(istat == FAIL)
		{
			if(!m_adaptStepSize)
			{
				// can't change the step size, so advection just ends
				return OKAY;
			}
			if(dt_attempt == m_fMinStepSize)
			{
				// tried to take a step with the min step size, but failed,
				// can't go any further
				return OKAY;
			}
			else
			{
				// try to retake the step with a smaller step size
				dt = dt_attempt * 0.1;
				if(dt < m_fMinStepSize)
					dt = m_fMinStepSize;
				thisParticle = prevParticle;
				prevParticle = second_prevParticle;
				second_prevParticle = third_prevParticle;
				curTime = prevCurTime;
				continue;
			}
		}

		VECTOR4 *p = new VECTOR4; 
		(*p)[0] = thisParticle.phyCoord[0]; 
		(*p)[1] = thisParticle.phyCoord[1]; 
		(*p)[2] = thisParticle.phyCoord[2]; 
		(*p)[3] = curTime; 
		seedTrace.push_back(p); 
		count++;

		if(!m_adaptStepSize)
		{
			// change step size to prevously used, since the oneStep methods
			// will change the value of dt
			dt = dt_attempt;
		}

		// check if point is outside real bounds 
		// (bounds not counting ghost cells)
		if(!m_pField->IsInRealBoundaries(thisParticle, curTime))
		{
			return OUT_OF_BOUND;
		}

		// check if point is at critical point
		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, 
		                  thisParticle, curTime, vel);
		if((fabs(vel[0]) < m_fStationaryCutoff) && 
		   (fabs(vel[1]) < m_fStationaryCutoff) && 
		   (fabs(vel[2]) < m_fStationaryCutoff))
			return CRITICAL_POINT;
	}

	return OKAY;
}
