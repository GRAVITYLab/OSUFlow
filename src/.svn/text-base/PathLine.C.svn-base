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
//
//      vtListTimeSeedTrace will contain the time information for the particle traces
// 
//////////////////////////////////////////////////////////////////////////
void vtCPathLine::execute(list<vtListTimeSeedTrace*>& listSeedTraces,
		list<int64_t> *listSeedIds)
{
	if (listSeedIds != NULL)
	{
		(*listSeedIds).clear();
	}
	listSeedTraces.clear();
	computePathLine(listSeedTraces, listSeedIds);
}



void vtCPathLine::computePathLine(list<vtListTimeSeedTrace*>& listSeedTraces, list<int64_t> *listSeedIds)
{
	int res;
	float currentT; 
	vtListParticleIter pIter = m_lSeeds.begin();
	std::list<int64_t>::iterator idIter = m_lSeedIds.begin();
	for(; pIter != m_lSeeds.end(); pIter++)
	{
		vtParticleInfo* thisSeed = *pIter;
		m_itsParticles.push_back(new vtParticleInfo(thisSeed));
	}

	// advect
	pIter = m_itsParticles.begin();
	int which_point = 0;
	for(; pIter != m_itsParticles.end(); pIter++)
	{
		vtListTimeSeedTrace *trace =  new vtListTimeSeedTrace;
		vtParticleInfo* thisParticle = *pIter;
		currentT = thisParticle->m_fStartTime; 

		if(thisParticle->itsValidFlag == 1)
		{
			res = advectParticle(m_integrationOrder, *thisParticle, currentT, FLT_MAX, *trace);
			if((int)trace->size() > 1)
			{
				listSeedTraces.push_back(trace);
				if (listSeedIds != NULL)
				{
					(*listSeedIds).push_back(*idIter);
				}
			}
		}
		if (listSeedIds != NULL)
		{
			idIter++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Get particle advecting one step, none-adaptive stepsize. For each seed,
// it only advects one step.
//////////////////////////////////////////////////////////////////////////
void vtCPathLine::execute(list<vtPathlineParticle*>& listSeedTraces)
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
		
	computePathLine(listSeedTraces);
}

void vtCPathLine::computePathLine(list<vtPathlineParticle*>& listSeedTraces)
{
        float currentT; 
	float finalT; 

	// advect the particles from this time step
	vtListParticleIter pIter = m_itsParticles.begin();
	for(; pIter != m_itsParticles.end(); pIter++)
	{
		vtParticleInfo *thisParticle = *pIter;
	
		if(thisParticle->itsNumStepsAlive > m_itsMaxParticleLife)
			continue;

		if(thisParticle->itsValidFlag == 1)
		{
		  currentT = thisParticle->m_fStartTime; 
		  finalT = currentT + m_timeDir * m_itsTimeInc;
		  int res = advectParticle(m_integrationOrder, *thisParticle, currentT, *thisParticle, finalT);
		  if(res == 1)
		    {
		      // for output
		      vtPathlineParticle* newParticle = new vtPathlineParticle;
		      newParticle->pos = thisParticle->m_pointInfo.phyCoord;
		      newParticle->ptId = thisParticle->ptId;
		      newParticle->time = finalT; 
		      listSeedTraces.push_back(newParticle);
		    }
		}
	}
}
