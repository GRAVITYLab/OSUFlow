#include "VectorMatrix.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// VECTOR3 class definitions
//////////////////////////////////////////////////////////////////////////

//added by lijie
MATRIX3 MATRIX3::transpose()
{
	MATRIX3 inv;
	for(int i=0; i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			inv[j][i]=mat[i][j];
		}
	}
	return inv;
}
//added by lijie
int MATRIX3::inverse(MATRIX3& m) 
{


    // 9 floating-point ops
    float d01d12md11d02 = mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2];
    float d01d22md21d02 = mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2];
    float d11d22md21d12 = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];

    // 5 floating-point ops
    float det =
      mat[0][0] * d11d22md21d12 -
      mat[1][0] * d01d22md21d02 +
      mat[2][0] * d01d12md11d02;

	if (0.0 == det) return 0;

    float det_inv = (float) 1.0 / det;

    // 19 floating-point ops
	MATRIX3 tm;
    tm[0]=VECTOR3(d11d22md21d12, -d01d22md21d02, d01d12md11d02);
	tm[1]=VECTOR3(mat[2][0] * mat[1][2] - mat[1][0] * mat[2][2],
			mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2],
			mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]);
	tm[2]=VECTOR3(mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1],
	     mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1],
	     mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]);

	m=tm*det_inv;

   
    return 1;


}
// normalize vector
void VECTOR3::Normalize()
{
	float norm = (float)sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	if(norm != 0.0)
	{
		for (int i = 0; i < Dimension(); i++)
			vec[i] = vec[i]/norm;
	}
}

// get the maximum value
float VECTOR3::GetMax() 
{
	float maxval = vec[0];
	if (vec[1] > maxval) maxval = vec[1];
	if (vec[2] > maxval) maxval = vec[2];
	return maxval;
}

// make sure all dimension <=1.0
void VECTOR3::Clamp() 
{
	for (int i = 0; i < Dimension(); i++)
		if (vec[i]>1.0) vec[i] = 1.0;
}

bool VECTOR3::IsSame(VECTOR3& a)
{
	for(int i = 0; i < Dimension(); i++)
		if(vec[i] != a(i))
			return false;
	return true;
}

void VECTOR3::scale(const float s)
{
	for(int i = 0; i < Dimension(); i++)
		vec[i] *= s;
}

void VECTOR3::minus(VECTOR3& v1, VECTOR3& v2)
{
	for(int i = 0; i < Dimension(); i++)
		vec[i] = v1(i) - v2(i);
}
//////////////////////////////////////////////////////////////////////////
// VECTOR4 class definitions
//////////////////////////////////////////////////////////////////////////
// normalize vector
void VECTOR4::Normalize()
{
	float norm = (float)sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	if(norm != 0.0)
	{
		for (int i = 0; i < Dimension(); i++)
			vec[i] = vec[i]/norm;
	}
}

//////////////////////////////////////////////////////////////////////////
// MATRIX3 class definitions
//////////////////////////////////////////////////////////////////////////

// set matrix to identity matrix
void MATRIX3::Identity()
{
	for (int i = 0; i < Dimension(); i++) {
		mat[i].Zero();
		mat[i][i] = 1.0;
	};
}

//////////////////////////////////////////////////////////////////////////
// MATRIX4 class definitions
//////////////////////////////////////////////////////////////////////////

// set matrix to identity matrix
void MATRIX4::Identity()
{
	for (int i = 0; i < Dimension(); i++) {
		mat[i].Zero();
		mat[i][i] = 1.0;
	};
}

//////////////////////////////////////////////////////////////////////////
// MATRIX3 operations
//////////////////////////////////////////////////////////////////////////

// return m0 + m1
MATRIX3 operator +(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	result[0] = m0(0) + m1(0);
	result[1] = m0(1) + m1(1);
	result[2] = m0(2) + m1(2);

	return(result);
}

// return m0 - m1
MATRIX3 operator -(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	result[0] = m0(0) - m1(0);
	result[1] = m0(1) - m1(1);
	result[2] = m0(2) - m1(2);

	return(result);
}

// return m0 * m1
MATRIX3 operator *(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	for (int i = 0; i < m0.Dimension(); i++)
		for (int j = 0; j < m0.Dimension(); j++) {
			result[i][j] = 0;
			for (int k = 0; k < m0.Dimension(); k++)
				result[i][j] += m0(i,k) * m1(k,j);
		};

	return(result);
}

// return x0 * m0
MATRIX3 operator *(const float x0, const MATRIX3 & m0)
{
	MATRIX3 result;

	result[0] = x0*m0(0);
	result[1] = x0*m0(1);
	result[2] = x0*m0(2);

	return(result);
}

// return m0 * x0
MATRIX3 operator *(const MATRIX3 & m0, const float x0)
{
	return(x0 * m0);
}

// return m0 * v0
VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0)
{
	VECTOR3 result;

	result[0] = dot(m0(0),v0);
	result[1] = dot(m0(1),v0);
	result[2] = dot(m0(2),v0);

	return(result);
}

// return v0 * m0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX3 & m0)
{
	VECTOR3 result;

	result[0] = v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0);
	result[1] = v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1);
	result[2] = v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2);

	return(result);
}

//////////////////////////////////////////////////////////////////////////
// MATRIX4 operations
//////////////////////////////////////////////////////////////////////////

// get the rotation matrix
// type: 1:x, 2:y, 3:z
MATRIX4 rotate_matrix(int type, float angle) {
	MATRIX4 rot_mat;
	switch(type) {
		case 1: // rotate by x
			{
				float anglex = (float)(angle/180.0)*(float)3.14159;
				float sinx = (float)sin((float)anglex);
				float cosx = (float)cos((float)anglex);
				rot_mat[1][1] = rot_mat[2][2] = cosx;
				rot_mat[1][2] = -sinx;
				rot_mat[2][1] = sinx;
				break;
			}
		case 2: // rotate by y
			{
				float angley = (float)(angle/180.0)*(float)3.14159;
				float siny = (float)sin((float)angley);
				float cosy = (float)cos((float)angley);
				rot_mat[0][0] = rot_mat[2][2] = cosy;
				rot_mat[0][2] = siny;
				rot_mat[2][0] = -siny;
				break;
			}
		case 3: // rotate by z
			{
				float anglez = (float)(angle/180.0)*(float)3.14159;
				float sinz = (float)sin((float)anglez);
				float cosz = (float)cos((float)anglez);
				rot_mat[0][0] = rot_mat[1][1] = cosz;
				rot_mat[0][1] = -sinz;
				rot_mat[1][0] = sinz;
				break;
			}
	}
	return rot_mat;
}

// get translation matrix
MATRIX4 translate_matrix(float tx, float ty, float tz){
	MATRIX4 translation_mat;
	translation_mat[0][3] = tx;
	translation_mat[1][3] = ty;
	translation_mat[2][3] = tz;
	return translation_mat;
}

// get scaling matrix
MATRIX4 scale_matrix(float sx, float sy, float sz){
	MATRIX4 scale_mat;
	scale_mat[0][0] = sx;
	scale_mat[1][1] = sy;
	scale_mat[2][2] = sz;
	return scale_mat;
}

// return m0 + m1
MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	result[0] = m0(0) + m1(0);
	result[1] = m0(1) + m1(1);
	result[2] = m0(2) + m1(2);
	result[3] = m0(3) + m1(3);

	return(result);
}

// return m0 - m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	result[0] = m0(0) - m1(0);
	result[1] = m0(1) - m1(1);
	result[2] = m0(2) - m1(2);
	result[3] = m0(3) - m1(3);

	return(result);
}

// return m0 * m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	for (int i = 0; i < m0.Dimension(); i++)
		for (int j = 0; j < m0.Dimension(); j++) {
			result[i][j] = 0;
			for (int k = 0; k < m0.Dimension(); k++)
				result[i][j] += m0(i,k) * m1(k,j);
		};

	return(result);
}

// return x0 * m0
MATRIX4 operator *(const float x0, const MATRIX4 & m0)
{
	MATRIX4 result;

	result[0] = x0*m0(0);
	result[1] = x0*m0(1);
	result[2] = x0*m0(2);
	result[3] = x0*m0(3);

	return(result);
}

// return m0 * x0
MATRIX4 operator *(const MATRIX4 & m0, const float x0)
{ return(x0*m0); };

// return m0 * v0
VECTOR4 operator *(const MATRIX4 & m0, const VECTOR4 & v0)
{
	VECTOR4 result;

	result[0] = dot(m0(0),v0);
	result[1] = dot(m0(1),v0);
	result[2] = dot(m0(2),v0);
	result[3] = dot(m0(3),v0);

	return(result);
}


// return m0 * v0
VECTOR3 operator *(const MATRIX4 & m0, const VECTOR3 & v0)
{
	VECTOR4 v(v0);
	VECTOR3 result;	

	float temp = dot(m0(3),v);
	result[0] = dot(m0(0),v)/temp;
	result[1] = dot(m0(1),v)/temp;
	result[2] = dot(m0(2),v)/temp;
	return(result);
}


VECTOR4 operator *(const VECTOR4 & v0, const MATRIX4 & m0)
// return v0 * m0
{
	VECTOR4 result;

	result[0] = v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0) + v0(3)*m0(3,0);
	result[1] = v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1) + v0(3)*m0(3,1);
	result[2] = v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2) + v0(3)*m0(3,2);
	result[3] = v0(0)*m0(0,3) + v0(1)*m0(1,3) + v0(2)*m0(2,3) + v0(3)*m0(3,3);

	return(result);
}

VECTOR3 operator *(const VECTOR3 & v0, const MATRIX4 & m0)
// return v0 * m0
{
	VECTOR3 result;
	float temp = v0(0)*m0(0,3) + v0(1)*m0(1,3) + v0(2)*m0(2,3) + m0(3,3);

	result[0] = (v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0) + m0(3,0))/temp;
	result[1] = (v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1) + m0(3,1))/temp;
	result[2] = (v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2) + m0(3,2))/temp;

	return(result);
}

//Code was taken from the original 'edge' library written by dave ebert.
MATRIX4 inverse(const MATRIX4 & m) {
	register int lp,i,j,k;
	static float wrk[4][8];
	static float a, b;
	MATRIX4 result;

	for( i=0; i<4; i++ )	/* Set up matrices */
	{
		for( j=0; j<4; j++ )
		{
			wrk[i][j]=(float)m(i,j);
			wrk[i][j+4]=0.0;
			result[i][j] = 0.0;
		}
		wrk[i][i+4]=1.0;
	}

	for( lp=0; lp<4; lp++ )	/* Loop over all rows */
	{
		a=0.0;
		j=(-1);
		for( i=lp; i<4; i++ )	/* Find largest non-zero element */
		{
			b=wrk[i][lp];
			if( b< 0.0 )
				b=(-b);
			if( b>a )
			{
				a=b;
				j=i;
			}
		}
		if( j!=lp )			/* If not on diagonal, put it there */
		{
			if( j<0 )		/* Singular if none found */
				return(result);
			else			/* Exchange rows from lp to end */
				for( k=lp; k<8; k++ )
				{	
					a=wrk[j][k];
					wrk[j][k]=wrk[lp][k];
					wrk[lp][k]=a;
				}
		}
		a=wrk[lp][lp];		/* Normalize working row */
		for( i=lp; i<8; i++ )
			wrk[lp][i]/=a;

		for( i=lp+1; i<8; i++ )  /* Adjust rest of work space */
		{
			b=wrk[lp][i];
			for( j=0; j<4; j++ )	/* One column at a time */
				if( j!=lp )
					wrk[j][i]-=wrk[j][lp]*b;
		}
	}

	for( i=0; i<4; i++ )	/* Return result matrix */
		for( j=0; j<4; j++ )
			result[i][j]=(float)wrk[i][j+4];
	return(result);
}


