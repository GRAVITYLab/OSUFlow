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
		if(int_order == SECOND)
			istat = runge_kutta2(m_timeDir, UNSTEADY, pt, &curTime, dt);
		else
			istat = runge_kutta4(m_timeDir, UNSTEADY, pt, &curTime, dt);

		if(istat != 1)
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
  printf(" hello!\n"); 
	int count = 0, istat, res;
	PointInfo seedInfo;
	PointInfo thisParticle, prevParticle, second_prevParticle;
	float dt, dt_estimate, cell_volume, mag, curTime;
	VECTOR3 vel;

	// the first particle
	seedInfo = initialPoint.m_pointInfo;
	res = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, initialTime, vel);
	if(res == OUT_OF_BOUND){
		return OUT_OF_BOUND;
	}
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff)) {
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
		dt = dt_estimate = m_fInitStepSize;
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

		if(int_order == SECOND)
			istat = runge_kutta2(m_timeDir, UNSTEADY, thisParticle, &curTime, dt);
		else
			istat = runge_kutta4(m_timeDir, UNSTEADY, thisParticle, &curTime, dt);

		if(istat == OUT_OF_BOUND)  {			// out of boundary
		  printf(" foul!  "); 
			return OUT_OF_BOUND;
		}

		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, curTime, vel);
		if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff)) {
		  printf("trap!! "); 
			return CRITICAL_POINT;			// arrives at a critical point
		}
			
		else
		{
			seedTrace.push_back(new VECTOR3(thisParticle.phyCoord));
			count++;
			printf(" count %d  ", count); 
		}

		if((curTime < m_pField->GetMinTimeStep()) || (curTime > (float)(m_pField->GetMaxTimeStep()))) {
		  printf("out!  "); 
			return -1;
		}

		//---> temporary turn off function of adaptive stepsize
		//if(count > 2)
		//	adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
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
	float dt, dt_estimate, cell_volume, mag, curTime;
	VECTOR3 vel;

	// the first particle
	seedInfo = initialPoint.m_pointInfo;
	res = m_pField->at_phys(seedInfo.fromCell, seedInfo.phyCoord, seedInfo, initialTime, vel);
	if(res == OUT_OF_BOUND)  {
		return OUT_OF_BOUND;
	}
	if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff)) {
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
	switch(m_pField->GetCellType())
	{
	case CUBE:
		dt = dt_estimate = m_fInitStepSize;
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

		if(int_order == SECOND)
			istat = runge_kutta2(m_timeDir, UNSTEADY, thisParticle, &curTime, dt);
		else
			istat = runge_kutta4(m_timeDir, UNSTEADY, thisParticle, &curTime, dt);
		
		{
		  VECTOR4 *p = new VECTOR4; 
		  (*p)[0] = thisParticle.phyCoord[0]; 
		  (*p)[1] = thisParticle.phyCoord[1]; 
		  (*p)[2] = thisParticle.phyCoord[2]; 
		  (*p)[3] = curTime; 
		  seedTrace.push_back(p); 
		  count++;
		}
		

		if(istat == OUT_OF_BOUND)			// out of boundary
			return OUT_OF_BOUND;

		m_pField->at_phys(thisParticle.fromCell, thisParticle.phyCoord, thisParticle, curTime, vel);
		if((fabs(vel[0]) < m_fStationaryCutoff) && (fabs(vel[1]) < m_fStationaryCutoff) && (fabs(vel[2]) < m_fStationaryCutoff))
			return CRITICAL_POINT;			// arrives at a critical point
		/*
		else
		{
		  VECTOR4 *p = new VECTOR4; 
		  (*p)[0] = thisParticle.phyCoord[0]; 
		  (*p)[1] = thisParticle.phyCoord[1]; 
		  (*p)[2] = thisParticle.phyCoord[2]; 
		  (*p)[3] = curTime; 
		  seedTrace.push_back(p); 
		  count++;
		}
		*/

		if((curTime < m_pField->GetMinTimeStep()) || (curTime > (float)(m_pField->GetMaxTimeStep()))) {
		  printf("**** time out  "); 
			return -1;
		}

		//---> temporary turn off function of adaptive stepsize
		//if(count > 2)
		//	adapt_step(second_prevParticle.phyCoord, prevParticle.phyCoord, thisParticle.phyCoord, dt_estimate, &dt);
	}

	return OKAY;
}


