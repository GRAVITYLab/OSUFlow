/*
THIS SOFTWARE BELONGS TO THE AMERICAN PEOPLE.  IT WAS DEVELOPED
AT THEIR EXPENSE AT NASA AMES RESEARCH CENTER.

NAME: Eigenvectors

REQUIREMENTS: Calculate the eigenvectors of a 3x3 matrix.

DISCUSSION:  This code used to have a hard time with certain
matrixes with all real eigenvalues.  It would find identical
eigenvectors for more than one of the eigenvalues.  This happened
because the algorithm would assign 1.0 to an element of the
eigenvector-to-be arbitrarily, so long as the resulting
determinant was non-zero.  As this process was not well controlled,
sometimes two eigenvectors would have 1.0 as the same element.
Al Globus added new code to be more robust about assigning the
1.0 to an element without replication.  Unfortunately, the
old method had to be kept around to find the one real eigenvector
when there is a complex conjugant pair of eigenvalues.

HISTORY
Modified for TOPOLOGY				Jan 06				Liya Li
Originated as part of eplane		????				Creon Levit
Modified for TOPOLOGY		 		late 90/early 91	Creon Levit
Modified for a more robust means
of finding all real eigenvectors	June 91				Al Globus
*/

#include "Topology.h"
#include "zero_test.h"

/*
* return the eigenvectors of a 3x3 matrix given the eigenvalues, assuming:
*
*  1) The egenvalues are all real
*  2) The matrix has three linearly independent eigenvectors
*
*/

/* used to calculate eigenvectors */
typedef struct
{
	double determinant[3];
	int useOne;		/* which memory of the eigenvector to set to 1 */
	double m[3][3];	/* matrix with eigenvalue subtracted from diagonal */
}tUsageStruct,*tUsage, SizeOfUsage;

#define UNASSIGNED -1

void AssignOnes( tUsageStruct usage[3] );
tUsage ThreeNewUsages( void );
void FreeThreeUsages( tUsage );
void ComputeDeterminants( tUsageStruct usage[3], float [3][3], float [3] );
void ComputeEigenVectors( tUsageStruct usage[3], float [3][3] );
int isUsed( tUsageStruct usage[3], int );
float determinant (float [3][3]);
void matrix_vector_product( float m[3][3], float x[3], float y[3] );
void normalize_vector( float v[3], float h );
void zero_small_components( float v[3] );


void compute_real_eigenvectors (float m[3][3], float vals[3], float vecs[3][3])
{
	tUsage usage = ThreeNewUsages();
	ComputeDeterminants( usage, m, vals );
	AssignOnes( usage );
	ComputeEigenVectors( usage, vecs );
	FreeThreeUsages( usage );
}

int isUsed( tUsageStruct usage[3], int which )
{
	int i;

	for( i = 0; i <= 2; i++ )
		if ( usage[i].useOne == which )
			return true;
	return false;
}

/*
assign useOnes: decide which element of each eigenvector is to be 1.0
Requirements:
Must not correspond with 0 determinant.
Should not correspond with another useOne.
Algorithm
1 If any eigenvalue has only one non-zero det, set useOne to that.
while everyone not assigned
2 If any remaining ones position has only one non-zero det, set
that useOne to that position.
3 Assign one with largest abs(det)
*/
void AssignOnes( tUsageStruct usage[3] )
{
	int u, d;				/* usage and determinant indexes */
	int areAssigned = 0;	/* how many 1.0s have been assigned */
	int whichu, whichd;		/* which usage and determinant is to have a 1 eigenvector member */
	int howMany;			/* number of cases that pass some test */
	double max;

	/* step 1 */
	for( u = 0; u <= 2; u++ )
	{
		howMany = 0; whichd = UNASSIGNED;
		for (d = 0; d <= 2; d++)
			if ( usage[u].determinant[d] != 0.0 )
			{ whichd = d; howMany++; }
			if ( howMany == 1 )
			{ usage[u].useOne = whichd; areAssigned++; assert(whichd != UNASSIGNED); }
	}

	while ( areAssigned < 3 )
	{
		/* step 2 */
		for( d = 0; d <= 2; d++ )
		{
			howMany = 0; whichu = UNASSIGNED;
			for (u = 0; u <= 2; u++)
				if ( (usage[u].useOne == UNASSIGNED) 
					&& (usage[u].determinant[d] != 0.0) 
					&& !isUsed(usage,d) )
				{ whichu = u; howMany++; }
				if ( howMany == 1 )
				{ usage[whichu].useOne = d; areAssigned++; assert(whichu != UNASSIGNED); }
		}

		/* step 3 */
		whichu = whichd = UNASSIGNED;
		max = -1;
		for( u = 0; u <= 2; u++ )
			for( d = 0; d <= 2; d++ )
				if ( (usage[u].useOne == UNASSIGNED) 
					&& (max <= fabs(usage[u].determinant[d])) 
					&& !isUsed(usage,d) )
				{ max = usage[u].determinant[d]; whichu = u; whichd = d; }
				if ( whichu != UNASSIGNED )
				{ usage[whichu].useOne = whichd; areAssigned++; }
	}
}

/* TS: this returns an ARRAY of three structs. Note the . dereference */
tUsage ThreeNewUsages( void )
{
	int i;
	tUsage self = (tUsage) malloc( 3 * sizeof(SizeOfUsage) );

	for( i = 0; i <= 2; i++ )
		self[i].useOne = UNASSIGNED;
	return self;
}

void FreeThreeUsages( tUsage self )
{
	free( self );
}

void ComputeDeterminants( tUsageStruct usage[3], float m[3][3], float vals[3] )
{
	int n, i,j;
	tUsage u;

	for (n = 0; n <= 2; n++)
	{
		u = &(usage[n]);

		for (i = 0; i <= 2; i++)
			for (j = 0; j <= 2; j++)
				u->m[i][j] = m[i][j];

		for (i = 0; i <= 2; i++)
		{
			u->m[i][i] -= vals[n];
			Zero_test2 (u->m[i][i], u->m[i][i], vals[n]);
		}

		for (i = 0; i <= 2; i++)
		{
			int o1 = (i + 1) % 3;
			int o2 = (i + 2) % 3;
			u->determinant[i] = u->m[o1][o1] * u->m[o2][o2] - u->m[o1][o2] * u->m[o2][o1];
			Zero_test4( u->determinant[i], u->m[o1][o1], u->m[o2][o2], u->m[o1][o2], u->m[o2][o1]);
		}
	}
}

void ComputeEigenVectors( tUsageStruct usage[3], float vecs[3][3] )
{
	int i;

	for( i = 0; i <= 2; i++ )
	{
		tUsage u = &(usage[i]);
		int one = u->useOne;
		int o1 = (u->useOne + 1) % 3;
		int o2 = (u->useOne + 2) % 3;

		assert(( 0 <= one) && (one <= 2 ));

		vecs[i][one] = 1.0;
		if ( u->determinant[one] != 0.0 )
		{
			vecs[i][o1] = ((u->m[o1][o2] * u->m[o2][one]) - (u->m[o2][o2] * u->m[o1][one])) / u->determinant[one];
			vecs[i][o2] = ((u->m[o2][o1] * u->m[o1][one]) - (u->m[o1][o1] * u->m[o2][one])) / u->determinant[one];
			zero_small_components (vecs[i]);
			normalize_vector (vecs[i], 1.0);
		}
		else
			vecs[i][o1] = vecs[i][o2] = 0.0;
	}

	if (determinant(vecs) == 0.0)
		printf ("Warning: real eigenvectors are not linearly independent!\n");
}

float L_infinity_norm (float m[3][3])
{
	float norm = 0.0, x;
	int i, j;
	for (i=0; i<2; i++)
		for (j=0; j<2; j++)
		{
			x = fabs(m[i][j]);
			if (x>norm)
				norm = x;
		}
	return norm;
}

float determinant (float m[3][3])
{
	double det = 
		-(m[0][2]*m[1][1]*m[2][0]) + m[0][1]*m[1][2]*m[2][0] +
		m[0][2]*m[1][0]*m[2][1]  - m[0][0]*m[1][2]*m[2][1] -
		m[0][1]*m[1][0]*m[2][2]  + m[0][0]*m[1][1]*m[2][2];
	float maxval = L_infinity_norm (m);
	Zero_test1 (det, maxval);
	return det;
}

void compute_real_eigenvector (float a[3][3], float val, float vec[3], int eigen_value_index)
{
	double det[3], m[3][3];
	int i,j;

	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			m[i][j] = a[i][j];

	for (i=0; i<3; i++)
	{
		m[i][i] = m[i][i] - val;
		Zero_test2 (m[i][i], m[i][i], val);
	}

	det[0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
	Zero_test4 (det[0], m[1][1],m[2][2],m[1][2],m[2][1]);

	det[1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
	Zero_test4 (det[1], m[0][0],m[2][2],m[0][2],m[2][0]);

	det[2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
	Zero_test4 (det[2], m[0][0],m[1][1],m[0][1],m[1][0]);

	switch (eigen_value_index)
	{
	case 0:
		if (det[0] != 0.0)
		{
			vec[0] = 1.0;
			vec[1] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) / det[0];
			vec[2] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det[0];
		}

		else if (det[1] != 0.0)
		{
			vec[0] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) / det[1];
			vec[1] = 1.0;
			vec[2] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) / det[1];
		}

		else if (det[2] != 0.0)
		{
			vec[0] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det[2];
			vec[1] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) / det[2];
			vec[2] = 1.0;
		}

		else
		{
			printf ("Warning: no nonzero determinant in compute_real_eigenvector\n");
			vec[0] = 1.0; vec[1]=0.0; vec[2]=0.0;
		}
		break;
	case 1:
		if (det[1] != 0.0)
		{
			vec[0] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) / det[1];
			vec[1] = 1.0;
			vec[2] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) / det[1];
		}

		else if (det[2] != 0.0)
		{
			vec[0] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det[2];
			vec[1] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) / det[2];
			vec[2] = 1.0;
		}
		else if (det[0] != 0.0)
		{
			vec[0] = 1.0;
			vec[1] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) / det[0];
			vec[2] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det[0];
		}
		else
		{
			printf ("Warning: no nonzero determinant in compute_real_eigenvector\n");
			vec[0] = 0.0; vec[1]=1.0; vec[2]=0.0;
		}
		break;
	case 2:
		if (det[2] != 0.0)
		{
			vec[0] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det[2];
			vec[1] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) / det[2];
			vec[2] = 1.0;
		}
		else if (det[0] != 0.0)
		{
			vec[0] = 1.0;
			vec[1] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) / det[0];
			vec[2] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det[0];
		}

		else if (det[1] != 0.0)
		{
			vec[0] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) / det[1];
			vec[1] = 1.0;
			vec[2] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) / det[1];
		}
		else
		{
			printf ("Warning: no nonzero determinant in compute_real_eigenvector\n");
			vec[0] = 0.0; vec[1]=0.0; vec[2]=1.0;
		}
		break;
	default:
		assert (0);
	}
	zero_small_components (vec);
	normalize_vector (vec, 1.0);
}

void compute_complex_eigenvectors( float m[3][3], float vals[3], float vecs[3][3] )
{
	/*
	* Find three vectors, the first, v[0], is a real eigenvector corresponding to the first (real)
	* eigenvalue.  The next two, v[1] and v[2],  span the two dimensional subspace of R^3 that is spanned by complex
	* linear cominations of the two complex conjugate eigenvectors corresponding to the two complex conjugate
	* eigenvalues.  If the complex eigenvectors are v and v*, then v[1] = Re{v} and v[2] = Im{v}.
	*/

	double det0, det1, det2;
	double a1, b1, a2, b2;
	double c0, c1, c2, d0, d1, d2;
	int i;

	compute_real_eigenvector (m, vals[0], vecs[0], 0);

	/* 
	* Find vec[1] and vec[2].  Mathematica hepled generate this code.  See file vecs.m
	*/

	c0 = -vals[2]*vals[2] + vals[1]*vals[1] - vals[1]*m[1][1] - m[1][2]*m[2][1] - vals[1]*m[2][2] + m[1][1]*m[2][2];
	Zero_test6 (c0, vals[2]*vals[2], vals[1]*vals[1], vals[1]*m[1][1], m[1][2]*m[2][1], vals[1]*m[2][2], m[1][1]*m[2][2]);

	d0 = 2*vals[2]*vals[1] - vals[2]*m[1][1] - vals[2]*m[2][2];
	Zero_test3 (d0, 2*vals[2]*vals[1], vals[2]*m[1][1], vals[2]*m[2][2]);

	det0 = c0*c0 + d0*d0;
	Zero_test2 (det0, c0*c0, d0*d0);

	c1 = -vals[2]*vals[2] + vals[1]*vals[1] - vals[1]*m[0][0] - m[0][2]*m[2][0] - vals[1]*m[2][2] + m[0][0]*m[2][2];
	Zero_test6 (c1, vals[2]*vals[2], vals[1]*vals[1], vals[1]*m[0][0], m[0][2]*m[2][0], vals[1]*m[2][2], m[0][0]*m[2][2]);

	d1 = 2*vals[2]*vals[1] - vals[2]*m[0][0] - vals[2]*m[2][2];
	Zero_test3 (d1, 2*vals[2]*vals[1], vals[2]*m[0][0], vals[2]*m[2][2]);

	det1 = c1*c1 + d1*d1;
	Zero_test2 (det1, c1*c1, d1*d1);

	c2 = -vals[2]*vals[2] + vals[1]*vals[1] - vals[1]*m[0][0] - m[0][1]*m[1][0] - vals[1]*m[1][1] + m[0][0]*m[1][1];
	Zero_test6 (c2, vals[2]*vals[2], vals[1]*vals[1], vals[1]*m[0][0], m[0][1]*m[1][0], vals[1]*m[1][1], m[0][0]*m[1][1]);

	d2 = 2*vals[2]*vals[1] - vals[2]*m[0][0] - vals[2]*m[1][1];
	Zero_test3 (d2, 2*vals[2]*vals[1], vals[2]*m[0][0], vals[2]*m[1][1]);

	det2 = c2*c2 + d2*d2;
	Zero_test2 (det2, c2*c2, d2*d2);

	if (det0 != 0.0) 
	{
		a1 = vals[1]*m[1][0] + m[1][2]*m[2][0] - m[1][0]*m[2][2];
		b1 = vals[2]*m[1][0];
		a2 = vals[1]*m[2][0] - m[1][1]*m[2][0] + m[1][0]*m[2][1];
		b2 = vals[2]*m[2][0];

		vecs[1][0] = 1;
		vecs[1][1] = (a1*c0 + b1*d0)/det0;
		vecs[1][2] = (d0*vals[2]*m[2][0] + c0*(vals[1]*m[2][0] - m[1][1]*m[2][0] + m[1][0]*m[2][1]))/det0;

		vecs[2][0] = 0;
		vecs[2][1] = (c0*vals[2]*m[1][0] - d0*(vals[1]*m[1][0] + m[1][2]*m[2][0] - m[1][0]*m[2][2]))/det0;
		vecs[2][2] = (b2*c0 - a2*d0)/det0;
	}

	else if (det1 != 0.0) 
	{
		a1 = vals[1]*m[0][1] + m[0][2]*m[2][1] - m[0][1]*m[2][2];
		b1 = vals[2]*m[0][1];
		a2 = vals[1]*m[2][1] - m[0][0]*m[2][1] + m[0][1]*m[2][0];
		b2 = vals[2]*m[2][1];

		vecs[1][1] = 1;
		vecs[1][0] = (a1*c1 + b1*d1)/det1;
		vecs[1][2] = (d1*vals[2]*m[2][1] + c1*(vals[1]*m[2][1] - m[0][0]*m[2][1] + m[0][1]*m[2][0]))/det1;

		vecs[2][1] = 0;
		vecs[2][0] = (c1*vals[2]*m[0][1] - d1*(vals[1]*m[0][1] + m[0][2]*m[2][1] - m[0][1]*m[2][2]))/det1;
		vecs[2][2] = (b2*c1 - a2*d1)/det1;
	}

	else if (det2 != 0.0) 
	{
		a1 = vals[1]*m[0][2] + m[0][1]*m[1][2] - m[0][2]*m[1][1];
		b1 = vals[2]*m[0][2];
		a2 = vals[1]*m[1][2] - m[0][0]*m[1][2] + m[0][2]*m[1][0];
		b2 = vals[2]*m[1][2];

		vecs[1][2] = 1;
		vecs[1][0] = (a1*c2 + b1*d2)/det2;
		vecs[1][1] = (d2*vals[2]*m[1][2] + c2*(vals[1]*m[1][2] - m[0][0]*m[1][2] + m[0][2]*m[1][0]))/det2;

		vecs[2][2] = 0;
		vecs[2][0] = (c2*vals[2]*m[0][2] - d2*(vals[1]*m[0][2] + m[0][1]*m[1][2] - m[0][2]*m[1][1]))/det2;
		vecs[2][1] = (b2*c2 - a2*d2)/det2;
	}

	else
	{
		printf ("Warning: no nonzero determinant in compute_complex_eigenvectors\n");
		vecs[1][0] = 0.0; vecs[1][1] = 1.0; vecs[1][2] = 0.0;
		vecs[2][0] = 0.0; vecs[2][1] = 0.0; vecs[2][2] = 1.0;
	}

	for (i=0; i<3; i++)
	{
		zero_small_components (vecs[i]);
		normalize_vector (vecs[i], 1.0);
	}

	if (determinant (vecs) == 0.0)
	{
		printf ("Warning: complex eigenvectors are not linearly independent!\n");
	}
}

void matrix_vector_product( float m[3][3], float x[3], float y[3] )
{
	y[0] = m[0][0]*x[0] + m[0][1]*x[1] + m[0][2]*x[2];
	y[1] = m[1][0]*x[0] + m[1][1]*x[1] + m[1][2]*x[2];
	y[2] = m[2][0]*x[0] + m[2][1]*x[1] + m[2][2]*x[2];
}


void normalize_vector( float v[3], float h )
{
	float mag = sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	float h_over_mag;
	if (mag==0.0)
		mag = 1e-9;
	h_over_mag = h/mag;
	v[0] = v[0]*h_over_mag;
	v[1] = v[1]*h_over_mag;
	v[2] = v[2]*h_over_mag;
}

void zero_small_components( float v[3] )
{
	int i;
	for (i=0; i<3; i++)
		Zero_test3 (v[i], v[0], v[1], v[2]);
}

