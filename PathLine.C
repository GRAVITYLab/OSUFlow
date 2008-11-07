/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 PathLines					  
//
///////////////////////////////////////////////////////////////////////////////

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

vtCPathLine::vtCPathLine(CVectorField* pField) : 
vtCTimeVaryingFieldLine(pField)
{
}

vtCPathLine::~vtCPathLine(void)
{
}

//////////////////////////////////////////////////////////////////////////
// Get the whole pathline, this can be adaptive stepsize. For each seed 
// from the initial time step, advects the seed until its trace is out
// of boundary (both spacial and time)
//
// output
//	listSeedTraces: for each seed, there is a list keeping all particles
//					it advects					
//////////////////////////////////////////////////////////////////////////
void vtCPathLine::execute(const void* userData,
						  list<vtListSeedTrace*>& listSeedTraces)
{
	listSeedTraces.clear();
	computePathLine(userData, listSeedTraces);
}

void vtCPathLine::computePathLine(const void* userData,
								  list<vtListSeedTrace*>& listSeedTraces)
{
	int res;
	float currentT = *(float *)userData;
	vtListParticleIter pIter = m_lSeeds.begin();
	for(; pIter != m_lSeeds.end(); pIter++)
	{
		vtParticleInfo* thisSeed = *pIter;
		m_itsParticles.push_back(new vtParticleInfo(thisSeed));
	}

	// advect
	pIter = m_itsParticles.begin();
	for(; pIter != m_itsParticles.end(); pIter++)
	{
		vtListSeedTrace *trace =  new vtListSeedTrace;
		vtParticleInfo* thisParticle = *pIter;
		
		if(thisParticle->itsValidFlag == 1)
		{
			res = advectParticle(m_integrationOrder, *thisParticle, currentT, FLT_MAX, *trace);
			if((int)trace->size() > 1)
				listSeedTraces.push_back(trace);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Get particle advecting one step, none-adaptive stepsize. For each seed,
// it only advects one step.
//////////////////////////////////////////////////////////////////////////
void vtCPathLine::execute(const void* userData, 
						  list<vtPathlineParticle*>& listSeedTraces)
{
	listSeedTraces.clear();

	if(m_itsParticles.size() == 0)
	{
		vtListParticleIter pIter = m_lSeeds.begin();
		for(; pIter != m_lSeeds.end(); pIter++)
		{
			vtParticleInfo* thisSeed = *pIter;
			m_itsParticles.push_back(new vtParticleInfo(thisSeed));
		}
	}
		
	computePathLine(userData, listSeedTraces);
}

void vtCPathLine::computePathLine(const void* userData,
								  list<vtPathlineParticle*>& listSeedTraces)
{
	float currentT = *(float*)userData;
	float finalT = currentT + m_timeDir * m_itsTimeInc;

	// advect the particles from this time step
	vtListParticleIter pIter = m_itsParticles.begin();
	for(; pIter != m_itsParticles.end(); pIter++)
	{
		vtParticleInfo *thisParticle = *pIter;
	
		if(thisParticle->itsNumStepsAlive > m_itsMaxParticleLife)
			continue;

		if(thisParticle->itsValidFlag == 1)
		{
			int res = advectParticle(m_integrationOrder, *thisParticle, currentT, *thisParticle, finalT);
			if(res == 1)
			{
				// for output
				vtPathlineParticle* newParticle = new vtPathlineParticle;
				newParticle->pos = thisParticle->m_pointInfo.phyCoord;
				newParticle->ptId = thisParticle->ptId;
				listSeedTraces.push_back(newParticle);
			}
		}
	}
}
