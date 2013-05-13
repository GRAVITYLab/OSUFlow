/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 timelines
//
///////////////////////////////////////////////////////////////////////////////


#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// class definition of timeline
//////////////////////////////////////////////////////////////////////////

vtCTimeLine::vtCTimeLine(CVectorField* pField) : 
vtCTimeVaryingFieldLine(pField),
m_itsTimeDelay(5)
{
	numTillRelease = m_itsTimeDelay;
}

vtCTimeLine::~vtCTimeLine(void)
{
	releaseSeedMemory();
}

void vtCTimeLine::setTimeDelay(int delay)
{
	m_itsTimeDelay = delay;
	numTillRelease = m_itsTimeDelay;
}

int vtCTimeLine::getTimeDelay(void)
{
	return m_itsTimeDelay;
}

//////////////////////////////////////////////////////////////////////////
// new particles are injected at a specified time interval
//////////////////////////////////////////////////////////////////////////
void vtCTimeLine::execute(const void* userData, vtListStreakParticle& listSeedTraces)
{
	listSeedTraces.clear();
	computeTimeLine(userData, listSeedTraces);
}

void vtCTimeLine::computeTimeLine(const void* userData, vtListStreakParticle& listSeedTraces)
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

	--numTillRelease;
	if((numTillRelease > 0) && (!m_itsParticles.empty()))
		return;
	else
		numTillRelease = m_itsTimeDelay;

	// advect the new particles from this time step
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
				// for output
				vtStreakParticle* newParticle = new vtStreakParticle;
				newParticle->itsTime = finalT;
				newParticle->itsPoint = nextP.m_pointInfo;
				listSeedTraces.push_back(newParticle);

				// for next timestep's advection
				m_itsParticles.push_back(new vtParticleInfo(nextP));
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

void vtCTimeLine::advectOldParticles( vtListParticleIter start, 
									   vtListParticleIter end, 
									   vtListStreakParticle& listSeedTraces,
									   float initialTime,
									   float finalTime,
									   vector<vtListParticleIter>& deadList)
{
	// advect the old particles first
	vtListParticleIter pIter = start;
	while(pIter != end)
	{
		vtParticleInfo* thisParticle = *pIter;

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
			listSeedTraces.push_back(newParticle);
			++thisParticle->itsNumStepsAlive;
			++pIter;
		}
		else
			deadList.push_back(pIter++);
	}
}
