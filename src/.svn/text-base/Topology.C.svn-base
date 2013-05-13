/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 
//
///////////////////////////////////////////////////////////////////////////////


#include "Topology.h"

#pragma warning(disable : 4251 4100 4244 4101)

CPTopology::CPTopology()
{
	m_pField = NULL;
}

CPTopology::CPTopology(CVectorField* pField)
{
	m_pField = pField;
}

CPTopology::~CPTopology()
{
}

//////////////////////////////////////////////////////////////////////////
// construct Jacobian matrix for 3D cartesian grid
//////////////////////////////////////////////////////////////////////////
void CPTopology::makeJacobian(tCriticalPoint& cp)
{
	int xdim, ydim, zdim;
	float xdim_1, ydim_1, zdim_1;
	m_pField->getDimension(xdim, ydim, zdim);
	xdim_1 = (float)(xdim-1);
	ydim_1 = (float)(ydim-1);
	zdim_1 = (float)(zdim-1);

	float stride = 1.0;
	VECTOR3 point = cp.physicalLocation;
	float min = 1000000.0f;
	float tmp = stride;
	if( point[0] - stride < 0.0f  )  tmp = point[0];				if(tmp < min)  min = tmp;
	if( point[0] + stride > xdim_1)  tmp = xdim_1  - point[0];		if(tmp < min)  min = tmp;
	if( point[1] - stride < 0.0f  )  tmp = point[1];				if(tmp < min)  min = tmp;
	if( point[1] + stride > ydim_1)  tmp = ydim_1  - point[1];		if(tmp < min)  min = tmp;
	if( point[2] - stride < 0.0f  )  tmp = point[2];				if(tmp < min)  min = tmp;
	if( point[2] + stride > zdim_1)  tmp = zdim_1  - point[2];		if(tmp < min)  min = tmp;

	if(tmp < stride)
		tmp = min;
	stride = tmp;

	VECTOR3  lftVec, rgtVec, botVec, topVec, fntVec, bckVec;
	m_pField->at_phys(VECTOR3(point[0]-stride, point[1],	    point[2]	   ), 0, lftVec);
	m_pField->at_phys(VECTOR3(point[0]+stride, point[1],		point[2]	   ), 0, rgtVec);
	m_pField->at_phys(VECTOR3(point[0],		   point[1]-stride, point[2]	   ), 0, botVec);
	m_pField->at_phys(VECTOR3(point[0],		   point[1]+stride, point[2]	   ), 0, topVec);
	m_pField->at_phys(VECTOR3(point[0],		   point[1],		point[2]-stride), 0, fntVec);
	m_pField->at_phys(VECTOR3(point[0],		   point[1],		point[2]+stride), 0, bckVec);

	stride = 0.5f / stride;
	cp.jacobian[0][0] = (rgtVec[0] - lftVec[0]) * stride;
	cp.jacobian[0][1] = (topVec[0] - botVec[0]) * stride;
	cp.jacobian[0][2] = (bckVec[0] - fntVec[0]) * stride;
	cp.jacobian[1][0] = (rgtVec[1] - lftVec[1]) * stride;
	cp.jacobian[1][1] = (topVec[1] - botVec[1]) * stride;
	cp.jacobian[1][2] = (bckVec[1] - fntVec[1]) * stride;
	cp.jacobian[2][0] = (rgtVec[2] - lftVec[2]) * stride;
	cp.jacobian[2][1] = (topVec[2] - botVec[2]) * stride;
	cp.jacobian[2][2] = (bckVec[2] - fntVec[2]) * stride;
}

int CPTopology::howManyPositiveEigenValues(tCriticalPoint cp)
{
	int ret = 0;
	int i;

	switch(cp.eigenType )
	{
	case ALL_REAL:
		for(i = 0; i <= 2; i++)
			if (cp.eigenValues[i] >= 0)		/* BUG: >= 0??? or > 0 */
				ret++;
		break;

	case TWO_COMPLEX:
		if (cp.eigenValues[0] >= 0)		/* BUG: >= 0??? or > 0 */
			ret++;
		if (cp.eigenValues[1] >= 0)		/* BUG: >= 0??? or > 0 */
			ret += 2;
		break;

	default:	
		assert( 0 );
	}

	return ret;
}

void CPTopology::putOddManOutEigenVectorFirst(tCriticalPoint& cp)
{
	float temp;
	int i;
	int numberPositiveEigenValues;

	numberPositiveEigenValues = howManyPositiveEigenValues(cp);
	if((cp.eigenType == ALL_REAL) && (numberPositiveEigenValues >= 1) && (numberPositiveEigenValues <= 2))
	{
		for(i = 0; i <= 2; i++)
			if((numberPositiveEigenValues == 1) ? (cp.eigenValues[i] >= 0) : (cp.eigenValues[i] < 0))
			{	/* swap */
				temp = cp.eigenValues[i];
				cp.eigenValues[i] = cp.eigenValues[0];
				cp.eigenValues[0] = temp;
				return;
			}
	}
}

t3dCPtype CPTopology::find3dType(tCriticalPoint cp)
{
	switch(cp.eigenType )
	{
	case ALL_REAL:
		switch(howManyPositiveEigenValues(cp))
		{
		case 0:		return AttractingNode;
		case 1:		return AttractingNodeSaddle;
		case 2:		return RepellingNodeSaddle;
		case 3:		return RepellingNode;
		default:	assert( 0 );
		}

	case TWO_COMPLEX:
		switch(howManyPositiveEigenValues(cp))
		{
		case 0:		return AttractingSpiral;			// no positive real part
		case 1:		return AttractingSpiralSaddle;		// the real root has real part
		case 2:		return RepellingSpiralSaddle;		// a pair of complex have real part
		case 3:		return RepellingSpiral;				// all have real part
		default:	assert( 0 );
		}

	default:
		return Degenerate3d;
	}
}

//////////////////////////////////////////////////////////////////////////
// classify 3D critical points
// positive real part signifies a repelling direction
// negative real part signifies an attracting direction
// imaginary part denotes circulation or spiral
//
// eigenvalues
// 1) all real roots
//		a. all positive		----	repelling node (source)
//		b. 2 positive, 1 negative	----	repelling saddle
//		c. 1 positive, 2 negative	----	attracting saddle
//		d. 0 positive	----	attracting node (sink)
//
// 2) 1 real root, a pair of complex conjugate roots
//		a. all positive		----	repelling spiral
//		b. 2 positive, 1 negative	----	spiral repelling saddle
//		c. 1 positive, 2 negative	----	spiral attracting saddle
//		d. 0 positive	----	spiral attracting
//////////////////////////////////////////////////////////////////////////
void CPTopology::CPClassify(tCriticalPoint& cp)
{
	int i;

	makeJacobian(cp);
	cp.eigenType = compute_eigenvalues(cp.jacobian, cp.eigenValues);
	if((cp.eigenType == ALL_REAL) || (cp.eigenType == TWO_COMPLEX))
	{
		putOddManOutEigenVectorFirst(cp);
		cp.cp3dType = find3dType(cp);
	}

	switch(cp.eigenType)
	{
	case ALL_REAL: 
		compute_real_eigenvectors(cp.jacobian, cp.eigenValues, cp.eigenVectors);
		break;

	case TWO_COMPLEX:
		compute_complex_eigenvectors(cp.jacobian, cp.eigenValues, cp.eigenVectors);
		break;

	default:
		cp.cp3dType = Degenerate3d;
		break;
	}
}

