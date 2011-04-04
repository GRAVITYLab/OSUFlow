/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 FieldLine
//
///////////////////////////////////////////////////////////////////////////////

#include "FieldLine.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// definition of class FieldLine
//////////////////////////////////////////////////////////////////////////
vtCFieldLine::vtCFieldLine(CVectorField* pField):
m_nNumSeeds(0),
m_integrationOrder(FOURTH),
m_timeDir(FORWARD),
m_fInitTime((float)0.0),
m_fStepTime((float)0.0),
m_fDurationTime((float)0.0),
m_fLowerAngleAccuracy((float)0.99),
m_fUpperAngleAccuracy((float)0.999),
m_fStationaryCutoff((float)0.00001), 
m_nMaxsize(MAX_LENGTH),
m_fInitStepSize(1.0),
m_pField(pField)
{
}

vtCFieldLine::~vtCFieldLine(void)
{
	releaseSeedMemory();
}

//////////////////////////////////////////////////////////////////////////
// release the memory allocated to seeds
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::releaseSeedMemory(void)
{
	vtListParticleIter pIter = m_lSeeds.begin();
	for( ; pIter != m_lSeeds.end(); ++pIter )
	{
		vtParticleInfo* thisPart = *pIter;
		delete thisPart;
	}
	m_lSeeds.erase(m_lSeeds.begin(), m_lSeeds.end() );
	m_lSeedIds.erase(m_lSeedIds.begin(), m_lSeedIds.end() );
	m_nNumSeeds = 0;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 2nd order Euler-Cauchy 
// predictor-corrector method. This routine is used for both steady and
// unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::euler_cauchy(TIME_DIR, TIME_DEP,float*, float)
{
	return 1;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 2th order Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta2(TIME_DIR time_dir, TIME_DEP time_dep, 
							   PointInfo& ci, 
							   float* t, float dt)
{
	int istat = 0;

	return istat;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 4th order Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta4(TIME_DIR time_dir, 
							   TIME_DEP time_dep, 
							   PointInfo& ci, 
							   float* t,			// initial time
							   float dt)			// stepsize
{
	int i, istat;
	VECTOR3 pt0;
	VECTOR3 vel;
	VECTOR3 k1, k2, k3;
	VECTOR3 pt;
	int fromCell;

	pt = ci.phyCoord;
	// 1st step of the Runge-Kutta scheme
	istat = m_pField->at_phys(ci.fromCell, pt, ci, *t, vel);
	if ( istat != 1 )
		return OUT_OF_BOUND;

	for( i=0; i<3; i++ )
	{
		pt0[i] = pt[i];
		k1[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k1[i]*(float)0.5;
	}

	// 2nd step of the Runge-Kutta scheme
	fromCell = ci.inCell;
	if ( time_dep  == UNSTEADY)
		*t += (float)0.5*time_dir*dt;
	
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel);
	if ( istat!= 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		k2[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k2[i]*(float)0.5;
	}

	// 3rd step of the Runge-Kutta scheme
	fromCell = ci.inCell;
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel);
	if ( istat != 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		k3[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k3[i];
	}

	//    4th step of the Runge-Kutta scheme
	if ( time_dep  == UNSTEADY)
		*t += (float)0.5*time_dir*dt;
	
	fromCell = ci.inCell;
	istat=m_pField->at_phys(fromCell, pt, ci, *t, vel);
	if ( istat != 1 )
	{
		ci.phyCoord = pt;
		return OUT_OF_BOUND;
	}
	
	for( i=0; i<3; i++ )
	{
		pt[i] = pt0[i]+(k1[i]+(float)2.0*(k2[i]+k3[i])+time_dir*dt*vel[i])/(float)6.0;
	}
	ci.phyCoord = pt;

	return( istat );
}

//////////////////////////////////////////////////////////////////////////
// aptive step size
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::adapt_step(const VECTOR3& p2,
							 const VECTOR3& p1,
							 const VECTOR3& p0,
							 float dt_estimate,
							 float* dt)
{
	float angle;
	VECTOR3 p2p1 = p2 - p1;
	VECTOR3 p1p0 = p1 - p0;
		
	angle = (float)acos(dot(p1p0, p2p1)/(p1p0.GetMag() * p2p1.GetMag()))*(float)RAD_TO_DEG;
	if(angle > m_fUpperAngleAccuracy)
		*dt = (*dt) * (float)0.5;
	else
		if(angle < m_fLowerAngleAccuracy)
		{
			*dt = (*dt) * (float)2.0;
			if(*dt >= m_fMaxStepSize)
				*dt = m_fMaxStepSize;
		}
	
	return true;
}

//////////////////////////////////////////////////////////////////////////
// initialize seeds
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR3* points, int numPoints, float t,
		int64_t *seedIds)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;
	
	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();

	if ( seedIds != NULL )
	{
		for (i = 0 ; i < numPoints ; i++ )
		{
			m_lSeedIds.push_back( seedIds[i] );
		}
	}

	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;

			newParticle->m_pointInfo.phyCoord = points[i];
			newParticle->m_fStartTime = t;
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information

			res = m_pField->at_phys(-1, points[i], newParticle->m_pointInfo, t, nodeData);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = points[i];
			thisSeed->m_fStartTime = t;
			res = m_pField->at_phys(-1, points[i], thisSeed->m_pointInfo, t, nodeData);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}

	m_nNumSeeds = numPoints;
}




//////////////////////////////////////////////////////////////////////////
// initialize seeds with possibly different start times in t array 
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR3* points, int numPoints, 
				 float* t)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;

	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();

	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;

			newParticle->m_pointInfo.phyCoord = points[i];
			newParticle->m_fStartTime = t[i];
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information

			res = m_pField->at_phys(-1, points[i], newParticle->m_pointInfo, t[i], nodeData);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = points[i];
			thisSeed->m_fStartTime = t[i];
			res = m_pField->at_phys(-1, points[i], thisSeed->m_pointInfo, t[i], nodeData);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}

	m_nNumSeeds = numPoints;
}



//////////////////////////////////////////////////////////////////////////
// initialize seeds with possibly different start times 
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::setSeedPoints(VECTOR4* points, int numPoints, 
		int64_t *seedIds)
{
	int i, res;
	VECTOR3 nodeData;

	if(points == NULL)
		return;

	// if the rake size has changed, forget the previous seed points
	if( m_nNumSeeds != numPoints )
		releaseSeedMemory();

	if ( seedIds != NULL )
	{
		for (i = 0 ; i < numPoints ; i++ )
		{
			m_lSeedIds.push_back( seedIds[i] );
		}
	}

	if( m_nNumSeeds == 0 )
	{
		for(i = 0 ; i < numPoints ; i++ )
		{
			vtParticleInfo* newParticle = new vtParticleInfo;

			VECTOR3 pos(points[i][0], points[i][1], points[i][2]); 
			newParticle->m_pointInfo.phyCoord = pos;
			newParticle->m_fStartTime = points[i][3]; 
			newParticle->ptId = i;

			// query the field in order to get the starting
			// cell interpolant for the seed point
			// which was passed to us without any prior information

			res = m_pField->at_phys(-1, pos, newParticle->m_pointInfo, points[i][3], nodeData);
// 			newParticle->itsValidFlag =  (res == 1) ? 1 : 0 ;
			newParticle->itsValidFlag =  1;
			m_lSeeds.push_back( newParticle );
		}
	} 
	else
	{
		vtListParticleIter sIter = m_lSeeds.begin();

		for(i = 0 ; sIter != m_lSeeds.end() ; ++sIter, ++i )
		{
			vtParticleInfo* thisSeed = *sIter;

			VECTOR3 pos(points[i][0], points[i][1], points[i][2]); 
			// set the new location for this seed point
			thisSeed->m_pointInfo.phyCoord = pos;
			thisSeed->m_fStartTime = points[i][3];
			res = m_pField->at_phys(-1, pos, thisSeed->m_pointInfo, points[i][3], nodeData);
			thisSeed->itsValidFlag =  (res == 1) ? 1 : 0 ;
		}    
	}

	m_nNumSeeds = numPoints;
}
