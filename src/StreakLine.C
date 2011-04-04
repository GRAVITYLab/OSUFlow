//*************************************************************************
//
//	OSU Flow Vector Field 
//	Created:	Han-Wei Shen, Liya Li	
//	The Ohio State University *	
//      Date:		06/2005
//	StreakLines						       
//
//**************************************************************************

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// class implementation of Numerical StreakLine
//////////////////////////////////////////////////////////////////////////

vtCStreakLine::vtCStreakLine(CVectorField* pField) : 
vtCTimeVaryingFieldLine(pField)
{
	nHowManyTraces = 0;
}

vtCStreakLine::~vtCStreakLine(void)
{
}

//////////////////////////////////////////////////////////////////////////
// this is an animated process
// input: 
// userData: current time
// listSeedTraces:	the output position of particles for the current time
//					step
//////////////////////////////////////////////////////////////////////////
void vtCStreakLine::execute(const void* userData,
			    vtStreakTraces& listSeedTraces)
{
	listSeedTraces.clear();

	// we advect streaklines step by step
	float currentT = *(float*)userData;
	//	computeStreakLine((void *)&currentT, listSeedTraces);

	while((currentT >= m_pField->GetMinTimeStep()) &&  
	      (currentT <= m_pField->GetMaxTimeStep()))
	{
		computeStreakLine((void *)&currentT, listSeedTraces);
		currentT += m_timeDir * m_itsTimeInc;
	}

}

void vtCStreakLine::computeStreakLine(const void* userData,
				      vtStreakTraces& listSeedTraces)
{
	float currentT = *(float*)userData;
	float finalT = currentT + m_timeDir * m_itsTimeInc;
	int numPreParticles;
	vector<vtListParticleIter> deadList;
	
	// how many particles alive from last time
	numPreParticles = (int)m_itsParticles.size();

	// whether the number of particles are more than the maximal allowed size
	if((numPreParticles + m_nNumSeeds) > m_nMaxsize)
	{
		vtListParticleIter pStart, pEnd;

		pStart = m_itsParticles.begin();
		pEnd = pStart;
		int killCount = numPreParticles + m_nNumSeeds - m_nMaxsize;

		for(int iFor = 0; iFor < killCount && pEnd != m_itsParticles.end(); iFor++, pEnd++)
		{
			vtParticleInfo* thisParticle = *pEnd;
			delete thisParticle;
		}

		m_itsParticles.erase(pStart, pEnd);
	}

	// advect the previous old particles
	deadList.clear();
	advectOldParticles(m_itsParticles.begin(), m_itsParticles.end(), listSeedTraces, currentT, finalT, deadList);

	// advect new particles from the seeds at this time step
	vtListParticleIter pIter = m_lSeeds.begin();
	for(; pIter != m_lSeeds.end(); pIter++)
	{
		vtParticleInfo *thisSeed = *pIter;
		vtParticleInfo nextP;

		if(thisSeed->itsValidFlag == 1)
		{
			int res = advectParticle(m_integrationOrder, *thisSeed, currentT, nextP, finalT);
			if(res == 1)
			{
				// every new particle starts a new trace
				vtListStreakParticle* pNewTrace = new vtListStreakParticle;
				vtStreakParticle *newParticle = new vtStreakParticle;
				vtStreakParticle *oldParticle = new vtStreakParticle;

				oldParticle->itsTime = currentT;
				oldParticle->itsPoint = thisSeed->m_pointInfo;
				oldParticle->traceId = nHowManyTraces;
				newParticle->itsTime = finalT;
				newParticle->itsPoint = nextP.m_pointInfo;
				newParticle->traceId = nHowManyTraces;

				pNewTrace->push_back(oldParticle);
				pNewTrace->push_back(newParticle);
				listSeedTraces.push_back(pNewTrace);
								
				// for next timestep's advection
				vtParticleInfo* activeParticle = new vtParticleInfo(nextP);
				activeParticle->ptId = nHowManyTraces;
				m_itsParticles.push_back(activeParticle);

				nHowManyTraces++;
			}
		}
	}
	// process those dead particles
	vector<vtListParticleIter>::iterator deadIter = deadList.begin();
	for(; deadIter != deadList.end(); deadIter++)
	{
		vtListParticleIter v = *deadIter;
		vtParticleInfo* pi = *v;
		delete pi;
		m_itsParticles.erase(v);
	}
}

//////////////////////////////////////////////////////////////////////////
// advect those particles from previous steps
//////////////////////////////////////////////////////////////////////////
void vtCStreakLine::advectOldParticles( vtListParticleIter start, 
					vtListParticleIter end, 
					vtStreakTraces& listSeedTraces,
					float initialTime,
					float finalTime,
					vector<vtListParticleIter>& deadList)
{
	// advect the old particles first
	vtListParticleIter pIter = start;
	int pathId;

	while(pIter != end)
	{
		vtParticleInfo* thisParticle = *pIter;
		pathId = thisParticle->ptId;
		
		// kill this particle, too old
		if(thisParticle->itsNumStepsAlive > m_itsMaxParticleLife)
		{
			deadList.push_back(pIter++);
			continue;
		}

		// advect more
		int res = advectParticle(m_integrationOrder, *thisParticle, initialTime, *thisParticle, finalTime);

		if(res == 1)
		{
			vtStreakParticle* newParticle = new vtStreakParticle;
			newParticle->itsTime = finalTime;
			newParticle->itsPoint = thisParticle->m_pointInfo;
			listSeedTraces[pathId]->push_back(newParticle);
			++thisParticle->itsNumStepsAlive;
			++pIter;
		}
		else
			deadList.push_back(pIter++);
	}
}
