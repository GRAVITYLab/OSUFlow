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
m_integrationOrder(RK45),
m_timeDir(FORWARD),
m_fInitTime((float)0.0),
m_fStepTime((float)0.0),
m_fDurationTime((float)0.0),
m_fLowerAngleAccuracy((float)0.99),
m_fUpperAngleAccuracy((float)0.999),
m_fStationaryCutoff((float)0.00001), 
m_nMaxsize(MAX_LENGTH),
m_fInitialStepSize(1.0),
m_fMinStepSize(0.01),
m_fMaxStepSize(5.0),
m_fMaxError(0.01),
m_adaptStepSize(true),
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
// set the lower and upper angle (used for adaptive step size)
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::SetLowerUpperAngle(float lowerAngle, float upperAngle)
{
  m_fLowerAngleAccuracy = lowerAngle;
  m_fUpperAngleAccuracy = upperAngle;
}

//////////////////////////////////////////////////////////////////////////
// get the lower and upper angle (used for adaptive step size)
//////////////////////////////////////////////////////////////////////////
void vtCFieldLine::GetLowerUpperAngle(float* lowerAngle, float* upperAngle)
{
  *lowerAngle = m_fLowerAngleAccuracy;
  *upperAngle = m_fUpperAngleAccuracy;
}

//////////////////////////////////////////////////////////////////////////
// Perform one proper integration step, taking into account error metrics. This
// is used for integration methods that rely on a geometric error analysis
// based on the angle greated by the new step.
//
// return FAIL or OKAY
//
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::oneStepGeometric(INTEG_ORD integ_ord, TIME_DIR time_dir, 
                                   TIME_DEP time_dep, PointInfo& thisParticle,
                                   PointInfo prevParticle, 
                                   PointInfo second_prevParticle, 
                                   float* curTime, float* dt, int count)
{
	int istat;
	PointInfo tempParticle = thisParticle;
	float prevTime = *curTime;
	int pass = 0;
	while(!pass)
	{
		// take a step
		if(integ_ord == SECOND)
			istat = runge_kutta2(time_dir, time_dep, thisParticle, curTime,*dt);
		else if(integ_ord == FOURTH)
			istat = runge_kutta4(time_dir, time_dep, thisParticle, curTime,*dt);
		else
			return OUT_OF_BOUND;

		if(istat == FAIL)
		  return FAIL;

		if(count >= 2)
		{
			pass = adapt_step(second_prevParticle.phyCoord,
			                  prevParticle.phyCoord, thisParticle.phyCoord, dt);

			if(!pass)
			{
				// revert to before the step was taken
				*curTime = prevTime;
				thisParticle = tempParticle;
			}
		}
		else
		{
			pass = 1;
		}
	}

	return OKAY;
}

//////////////////////////////////////////////////////////////////////////
// Perform one proper integration step, taking into account error metrics.
// This is used for integration methods that have an embedded higher order
// method, and uses it to derive its error value. For an explanation of how the
// step size is increase or decreased, please see the book 
// "Numerical Recipes".
//
// return FAIL or OKAY
//
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::oneStepEmbedded(INTEG_ORD integ_ord, TIME_DIR time_dir, 
                                  TIME_DEP time_dep, PointInfo& thisParticle, 
                                  float* curTime, float* dt)
{
	// constants concerning how much the step size grows and shrinks
	const float pgrow = -0.2;
	const float pshrink = -0.25;
	const float safety = 0.9;

	// this value is equal to (5/safety) raised to the power (1/pgrow).
	// see use below.
	const float errcon = 1.89e-4; 

	int istat;
	PointInfo tempParticle = thisParticle;
	float prevTime = *curTime;
	float errmax;
	float h = *dt;

	// prevent input stepsize out of range
	h = max(h, m_fMinStepSize);
	h = min(h, m_fMaxStepSize);

	while(true)
	{
		if(integ_ord == RK45)
			istat = runge_kutta45(time_dir, time_dep, thisParticle, 
								  curTime, h, &errmax);
		else
			return OUT_OF_BOUND;

		if(istat == FAIL)
		  return FAIL;

		// convert errmax to a ratio to see how close it is to our maximum
		// allowed error
		errmax = errmax / m_fMaxError;

		if(errmax <= 1.0)
			// error is small enough, step succeeded
			break;

		if(h == m_fMinStepSize)
			// just took a step with the min step size, 
			// so automatically accept this step
			break;

		// error too large, make step size smaller.
		// don't decrease any more than 10 times smaller.
		float htemp = safety * h * pow(errmax, pshrink);
		if(h >= 0.0)
			h = max(htemp, 0.1f * h);
		else
			h = min(htemp, 0.1f * h);

		h = max(h, m_fMinStepSize);
        h = min(h, m_fMaxStepSize);

		thisParticle = tempParticle;
		*curTime = prevTime;
	}

	// step succeeded, increase step size. 
	// increase it no more than a factor of 5.
	float hnext;
	if(h == m_fMinStepSize && errmax > 1.0)
	{
		// step size should not increase
		hnext = h;
	}
	else
	{
		if(errmax > errcon)
			hnext = safety * h * pow(errmax, pgrow);
		else
			hnext = 5.0f * h;
	}

	// get variables ready for next step
	hnext = max(hnext, m_fMinStepSize);
	hnext = min(hnext, m_fMaxStepSize);
	*dt = hnext;

	return OKAY;
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
//
// time_dir: backwards or forwards integration (value is -1 or 1)
// time_dep: boolean indicating whether doing time dependent integration
// ci: the current point. will hold the new point at end of the function.
// t: the initial time. will be modified to be the correct time by the end of
//    the integration.
// dt: the step size
//
// returns FAIL or OKAY
//
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta2(TIME_DIR time_dir, 
                               TIME_DEP time_dep, 
                               PointInfo& ci, 
                               float* t,
                               float dt)
{
	int i, istat;
	VECTOR3 pt0;
	VECTOR3 vel;
	VECTOR3 k1, k2;
	VECTOR3 pt;
	int fromCell;

	pt = ci.phyCoord;

	// 1st step of the Runge-Kutta scheme
	istat = m_pField->at_phys(ci.fromCell, pt, ci, *t, vel);
	if(istat != OKAY)
		return FAIL;

	for(i=0; i<3; i++)
	{
		pt0[i] = pt[i];
		k1[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i]+k1[i];    // x0 + k1
	}

	// 2nd step of the Runge-Kutta scheme
	fromCell = ci.inCell;
	if(time_dep  == UNSTEADY)
		*t += time_dir*dt;    // t + h

	istat = m_pField->at_phys(fromCell, pt, ci, *t, vel);
	if(istat != OKAY)
	{
		ci.phyCoord = pt;
		return FAIL;
	}

	for( i=0; i<3; i++ )
	{
		k2[i] = time_dir*dt*vel[i];
		pt[i] = pt0[i] + (float)0.5*(k1[i] + k2[i]);  // new point
	}
	ci.phyCoord = pt;

	return OKAY;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using the 4th order Runge-Kutta method.
// This routine is used for both steady and unsteady vector fields. 
//
// time_dir: backwards or forwards integration (value is -1 or 1)
// time_dep: boolean indicating whether doing time dependent integration
// ci: the current point. will hold the new point at end of the function.
// t: the initial time. will be modified to be the correct time by the end of
//    the integration.
// dt: the step size
//
// returns FAIL or OKAY
//
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta4(TIME_DIR time_dir, 
							   TIME_DEP time_dep, 
							   PointInfo& ci, 
							   float* t,
							   float dt)
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
		return FAIL;

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
		return FAIL;
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
		return FAIL;
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
		return FAIL;
	}
	
	for( i=0; i<3; i++ )
	{
		pt[i] = pt0[i]+(k1[i]+(float)2.0*(k2[i]+k3[i])+time_dir*dt*vel[i])/(float)6.0;
	}
	ci.phyCoord = pt;

	return OKAY;
}

//////////////////////////////////////////////////////////////////////////
// Integrate along a field line using an embedded 4th and 5th order Runge-Kutta
// method. The actual method used is the Cash-Karp Runge-Kutta. An estimation
// of the error is also computed. This routine is used for both steady and
// unsteady vector fields.
//
// time_dir: backwards or forwards integration (value is -1 or 1)
// time_dep: boolean indicating whether doing time dependent integration
// ci: the current point. will hold the new point at end of the function.
// t: the initial time. will be modified to be the correct time by the end of
//    the integration.
// dt: the step size
// error: will hold the calculated error at the end of the function (output)
//
// returns FAIL or OKAY
//
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::runge_kutta45(TIME_DIR time_dir,
                                TIME_DEP time_dep,
                                PointInfo& ci,
                                float* t,			// initial time
                                float dt,			// step size
								float* error)       // error
{
	int i, istat;
	float t_tmp;
	VECTOR3 vel;
	VECTOR3 k1, k2, k3, k4, k5, k6;
	VECTOR3 pt0;	// the original point
	VECTOR3 pt;		// holds the last point at which the field was sampled
	int fromCell;

	// factors used for this runge-kutta
	// Note that c2, c5, and dc2 are zero
	float a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
	float b21 = 0.2; 
	float b31 = 3.0/40.0, b32 = 9.0/40.0;
	float b41 = 0.3, b42 = -0.9, b43 = 1.2;
	float b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0;
	float b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0;
	float b64 = 44275.0/110592.0, b65 = 253.0/4096.0;
	float c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0; 
	float c6 = 512.0/1771.0; 
	float dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0;
	float dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;

	pt = ci.phyCoord;
	pt0 = ci.phyCoord;
	t_tmp = *t;

	// 1st step
	istat = m_pField->at_phys(ci.fromCell, pt, ci, *t, vel);
	if ( istat != 1 ) {
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k1[i] = time_dir*dt*vel[i];
	}

	// 2nd step
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + k1[i]*b21;
	}
	if ( time_dep  == UNSTEADY ) {
		t_tmp = *t + a2*time_dir*dt;
	}
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, t_tmp, vel);
	if ( istat!= 1 ) {
		ci.phyCoord = pt;
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k2[i] = time_dir*dt*vel[i];
	}

	// 3rd step
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + k1[i]*b31 + k2[i]*b32;
	}
	if ( time_dep  == UNSTEADY ) {
		t_tmp = *t + a3*time_dir*dt;
	}
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, t_tmp, vel);
	if ( istat != 1 ) {
		ci.phyCoord = pt;
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k3[i] = time_dir*dt*vel[i];
	}

	// 4th step
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + k1[i]*b41 + k2[i]*b42 + k3[i]*b43;
	}
	if ( time_dep  == UNSTEADY ) {
		t_tmp = *t + a4*time_dir*dt;
	}
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, t_tmp, vel);
	if ( istat != 1 ) {
		ci.phyCoord = pt;
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k4[i] = time_dir*dt*vel[i];
	}

	// 5th step
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + k1[i]*b51 + k2[i]*b52 + k3[i]*b53 + k4[i]*b54;
	}
	if ( time_dep  == UNSTEADY ) {
		t_tmp = *t + a5*time_dir*dt;
	}
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, t_tmp, vel);
	if ( istat != 1 ) {
		ci.phyCoord = pt;
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k5[i] = time_dir*dt*vel[i];
	}

	// 6th step
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + k1[i]*b61 + k2[i]*b62 + k3[i]*b63 + k4[i]*b64 + 
	                                                         k5[i]*b65;
	}
	if ( time_dep  == UNSTEADY ) {
		t_tmp = *t + a6*time_dir*dt;
	}
	fromCell = ci.inCell;
	istat = m_pField->at_phys(fromCell, pt, ci, t_tmp, vel);
	if ( istat != 1 ) {
		ci.phyCoord = pt;
		return FAIL;
	}
	for( i=0; i<3; i++ ) {
		k6[i] = time_dir*dt*vel[i];
	}

	// the new point (resulting from the fifth order runge-kutta)
	for( i=0; i<3; i++ ) {
		pt[i] = pt0[i] + c1*k1[i] + c3*k3[i] + c4*k4[i] + c6*k6[i];
	}

	// error, which is the difference between the fifth order and the fourth
	// order answers. this is the vector from the 5th order and 4th order
	// answers.
	VECTOR3 err;
	for( i=0; i<3; i++ ) {
		err[i]= dc1*k1[i] + dc3*k3[i] + dc4*k4[i] + dc5*k5[i] + dc6*k6[i];
	}

	// total error is magnitude of err
	*error = sqrt(err[0]*err[0] + err[1]*err[1] + err[2]*err[2]);
	if ( time_dep  == UNSTEADY )
		*t += dt;

	ci.phyCoord = pt;
	return OKAY;
}

//////////////////////////////////////////////////////////////////////////
// adaptive step size
// returns 0 if the step size had to be reduced. otherwise returns 1, which
// means step size was kept the same or increased.
//
// Basically look at the angle formed by the current point and the two previous
// points. The smaller the angle formed, the more accurate it is. If the angle
// formed by them is too great, then the step size is decreased. If the angle
// is smaller than the min allowed angle, then the step size is increased. 
// This error method was taken from "Interactive Time-Dependent Particle
// Tracing Using Tetrahedral Decomposition" by Kenwright and Lane.
//////////////////////////////////////////////////////////////////////////
int vtCFieldLine::adapt_step(const VECTOR3& p2,
                             const VECTOR3& p1,
                             const VECTOR3& p0,
                             float* dt)
{
	float angle;
	VECTOR3 p2p1 = p2 - p1;
	VECTOR3 p1p0 = p1 - p0;
		
	angle = (float)acos(dot(p1p0, p2p1)/(p1p0.GetMag() * p2p1.GetMag()))*(float)RAD_TO_DEG;
	if(angle > m_fUpperAngleAccuracy)
	{
		if(*dt == m_fMinStepSize)
		{
			// if this is true, i know that the last step was taken using the
			// min step size, so returning fail will only repeat the same step,
			// so pass
			return 1;
		}

		*dt = (*dt) * 0.5f;

		if(*dt < m_fMinStepSize)
			*dt = m_fMinStepSize;
		if(*dt > m_fMaxStepSize)
			*dt = m_fMaxStepSize;
		return 0;
	}
	else if(angle < m_fLowerAngleAccuracy)
	{
		*dt = (*dt) * 2.0f;
	}
	
	if(*dt < m_fMinStepSize)
		*dt = m_fMinStepSize;
	if(*dt > m_fMaxStepSize)
		*dt = m_fMaxStepSize;

	return 1;
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

