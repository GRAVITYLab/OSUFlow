//////////////////////////////////////////////////////////////////////////
// Vector and Matrix Class
//
// Created by: Matt Camuto
//
// Modification:
//		Time            Programmer
//		07-09-99        R. Wenger
//		08-20-99        K. Boner
//		07-15-00        J. Gao
// 		05-25-05		H-W Shen, Liya Li
//////////////////////////////////////////////////////////////////////////

#ifndef _VECTOR_MATRIX_H_
#define _VECTOR_MATRIX_H_

#include "header.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

typedef float SCALAR;

//////////////////////////////////////////////////////////////////////////
// 2d vector used for Point
//////////////////////////////////////////////////////////////////////////
class VECTOR2 {

private :
	float vec[2];

public :
	VECTOR2()                                    // constructor
	{ Zero(); };
	VECTOR2(const float x0, const float x1) // constructor
	{ vec[0] = x0; vec[1] = x1;};
	int Dimension() const
	{ return 2; };
	float & operator [](const int i)           // index i'th element
	{ return(vec[i]); };
	float operator ()(const int i) const        // return i'th element
	{ return(vec[i]); };
	const VECTOR2 & operator =(const VECTOR2 & v0)  // copy vector v0
	{ vec[0] = v0(0); vec[1] = v0(1); 
	return(*this); };
	void Zero()                                  // make zero vector
	{ vec[0] = vec[1] = 0.0; };
	void Set(const float x0, const float x1)
	{ vec[0] = x0; vec[1] = x1; };
};

// Line class
typedef struct line_2d {
	VECTOR2 start;
	VECTOR2 end;
}Line;

//////////////////////////////////////////////////////////////////////////
//	vector with 3 components
//////////////////////////////////////////////////////////////////////////
class VECTOR3
{
private :
	float vec[3];

public :
	// constructor
	VECTOR3() {Zero();}
	VECTOR3(const float x0, const float x1, const float x2) {vec[0] = x0; vec[1] = x1; vec[2] = x2;}
	// ADD-BY-LEETEN 02/07/2011-BEGIN
	VECTOR3(const VECTOR3& v)
	{
		for(int i = 0; i < 3; i++)
			vec[i] = v.vec[i];
	}
	// ADD-BY-LEETEN 02/07/2011-END

	void Set(const float x0, const float x1, const float x2){vec[0] = x0; vec[1] = x1; vec[2] = x2;}
	int Dimension() const {	return 3; }
	float x() {return vec[0];}
	float y() {return vec[1];}
	float z() {return vec[2];}
	float & operator [](const int i) {return(vec[i]);}				// index i'th element
	// ADD-BY-LEETEN 02/07/2011-BEGIN
	const float & operator [](const int i) const {return(vec[i]);}				// index i'th element
	// ADD-BY-LEETEN 02/07/2011-END

	float operator ()(const int i) const {return(vec[i]);}			// return i'th element
	bool operator ==(const VECTOR3& v) const
	{
		if((fabs(vec[0]-v(0))<EPS) && (fabs(vec[1]-v(1))<EPS) && (fabs(vec[2]-v(2))<EPS))
			return true;
		else
			return false;
	}
	void Zero() {vec[0] = vec[1] = vec[2] = 0.0;}					// make zero vector
	float GetMag() {return (float)(sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]));}	// get magnitude
	const VECTOR3 & operator =(const VECTOR3 & v0)					// copy vector v0
	{vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); return(*this);}
	float GetMax();													// get the maximum value								
	void Clamp();													// make sure all dimension <=1.0
	void Normalize();					                            // normalize vector
	bool IsSame(VECTOR3& a);							
	void scale(const float s);
	void minus(VECTOR3& v1, VECTOR3& v2);
	void add(float a, float b, float c) { vec[0]+= a; vec[1] += b; vec[2] += c; }
};

//////////////////////////////////////////////////////////////////////////
// 4d vector class
//////////////////////////////////////////////////////////////////////////
class VECTOR4
{
private :
	float vec[4];

public :
	// constructor
	VECTOR4() {Zero();}
	// MOD-BY-LEETEN 02/07/2011-FROM:
		// VECTOR4(VECTOR3 v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = 1.0;}
	// TO:
	VECTOR4(const VECTOR3& v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = 1.0;}
	VECTOR4(const VECTOR4& v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = v[3];}
	// MOD-BY-LEETEN 02/07/2011-END
	VECTOR4(const float x0, const float x1, const float x2, const float x3)
	{vec[0] = x0; vec[1] = x1; vec[2] = x2; vec[3] = x3;}

	void Set(const float x0, const float x1, const float x2, const float x3)
	{vec[0] = x0; vec[1] = x1; vec[2] = x2; vec[3] = x3;}

	int Dimension() const {	return 4;}
	float & operator [](const int i) { return(vec[i]);}			// index i'th element
	// ADD-BY-LEETEN 02/07/2011-BEGIN
	const float & operator [](const int i) const { return(vec[i]);}			// index i'th element
	// ADD-BY-LEETEN 02/07/2011-END
	float operator ()(const int i) const {return(vec[i]);}		// return i'th element
	const VECTOR4 & operator =(const VECTOR4 & v0)				// copy vector v0
	{vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = v0(3);return(*this);}
	const VECTOR4 & operator =(const VECTOR3 & v0)				// copy vector v0
	{vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = 1.0; return(*this);}
	VECTOR3 get_vector3() {VECTOR3 temp(vec[0]/vec[3], vec[1]/vec[3], vec[2]/vec[3]); return temp;}
	void Zero() {vec[0] = vec[1] = vec[2] = vec[3] = 0.0;}		// make zero vector

	void Normalize();											// normalize vector
};

//////////////////////////////////////////////////////////////////////////
// 3d matrix
//////////////////////////////////////////////////////////////////////////
class MATRIX3
{
private :
	VECTOR3 mat[3];       // a vector represents each matrix row

public :

	MATRIX3()                                    // constructor
	{	Identity(); };
	MATRIX3(const VECTOR3 & v0, const VECTOR3 & v1, const VECTOR3 & v2)
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; };  // constructor
	int Dimension() const
	{	return 3; };
	VECTOR3 & operator [](const int i)           // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX3 m0 as m0[i][j] (not m0[i,j])
	VECTOR3 operator()(const int i) const        // return row i
	{	return(mat[i]); };
	float operator ()(const int i, const int j) const   
	{	return(mat[i](j)); };                    // return element (i,j)
	MATRIX3 & operator =(const MATRIX3 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); 
	return(*this); };
	void Identity();                             // set to identity

	int inverse( MATRIX3& m) ;//added by lijie to handle curvilinear grid
	MATRIX3 transpose();//added by lijie to handle curvilinear grid

};

//////////////////////////////////////////////////////////////////////////
// 4d matrix
//////////////////////////////////////////////////////////////////////////
class MATRIX4
{

private :
	VECTOR4 mat[4];       // a vector represents each matrix row

public :		

	MATRIX4()                                    // constructor
	{	Identity(); };
	MATRIX4(const VECTOR4 & v0, const VECTOR4 & v1, 
		const VECTOR4 & v2, const VECTOR4 & v3) // constructor
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; mat[3] = v3; };  
	int Dimension() const
	{	return 4; };
	VECTOR4 & operator [](int i)                 // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX4 m0 as m0[i][j] (not m0[i,j])
	VECTOR4 operator()(const int i) const        // return row i
	{	return(mat[i]); };
	float operator ()(const int i, const int j) const
	{	return(mat[i](j)); };                    // return element (i,j)
	MATRIX4 & operator =(const MATRIX4 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); mat[3] = m0(3);
	return(*this); };
	MATRIX4 & operator =(const MATRIX3 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); 
	VECTOR4 temp(0.0,0.0,0.0,1.0);
	mat[3] = temp;
	return(*this); };
	void Identity();                             // set to identity
};

//************************
// VECTOR2 operations
//************************

inline VECTOR2 operator +(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 + v1
{	return(VECTOR2(v0(0) + v1(0), v0(1) + v1(1))); };

inline VECTOR2 operator -(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 - v1
{	return(VECTOR2(v0(0) - v1(0), v0(1) - v1(1))); };

inline VECTOR2 operator *(float x0, const VECTOR2 & v0)
// return x0*v0
{	return(VECTOR2(x0*v0(0), x0*v0(1))); };

inline VECTOR2 operator *(const VECTOR2 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };


//************************
// VECTOR3 operations
//************************

inline float dot(const VECTOR3 & v0, const VECTOR3 & v1)   
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2)); };

inline VECTOR3 cross(const VECTOR3 & v0, const VECTOR3 & v1)
// return cross product of v0 and v1
{	return(VECTOR3(v0(1)*v1(2) - v0(2)*v1(1), 
		   v0(2)*v1(0) - v0(0)*v1(2),
		   v0(0)*v1(1) - v0(1)*v1(0))); };

inline VECTOR3 operator +(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 + v1
{	return(VECTOR3(v0(0) + v1(0), v0(1) + v1(1), v0(2) + v1(2))); };

inline VECTOR3 operator -(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 - v1
{	return(VECTOR3(v0(0) - v1(0), v0(1) - v1(1), v0(2) - v1(2))); };

inline float operator *(VECTOR3& v0, VECTOR3& v1)
// return v0*v1T
{	return (v0(0) * v1(0) + v0(1) * v1(1) + v0(2) * v0(2)); };

inline VECTOR3 operator *(float x0, const VECTOR3 & v0)
// return x0*v0
{	return(VECTOR3(x0*v0(0), x0*v0(1), x0*v0(2))); };

inline VECTOR3 operator *(const VECTOR3 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

inline VECTOR4 get_vector4(VECTOR3 vec)
{	VECTOR4 temp(vec[0], vec[1], vec[2], 1.0); return temp;};

//************************
// VECTOR4 operations
//************************

inline float dot(const VECTOR4 & v0, const VECTOR4 & v1)   
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)); };

inline VECTOR4 operator +(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 + v1
{	return(VECTOR4(v0(0)+v1(0), v0(1)+v1(1), v0(2)+v1(2), v0(3)+v1(3))); };

inline VECTOR4 operator -(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 - v1
{	return(VECTOR4(v0(0)-v1(0), v0(1)-v1(1), v0(2)-v1(2), v0(3)-v1(3))); };

inline VECTOR4 operator *(float x0, const VECTOR4 & v0)
// return x0*v0
{	return(VECTOR4(x0*v0(0), x0*v0(1), x0*v0(2), x0*v0(3))); };

inline VECTOR4 operator *(const VECTOR4 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

//************************
// MATRIX3 operations
//************************

MATRIX3 operator +(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 + m1
MATRIX3 operator -(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 - m1
MATRIX3 operator *(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 * m1
MATRIX3 operator *(const float x0, const MATRIX3 & m0);    // return x0 * m0
MATRIX3 operator *(const MATRIX3 & m0, const float x0);    // return m0 * x0
VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX3 & m0); // return v0 * m0

//************************
// MATRIX4 operations
//************************

MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 + m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 - m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 * m1
MATRIX4 operator *(const float x0, const MATRIX4 & m0);    // return x0 * m0
MATRIX4 operator *(const MATRIX4 & m0, const float x0);    // return m0 * x0
VECTOR4 operator *(const MATRIX4 & m0, const VECTOR4 & v0); // return m0 * v0
VECTOR4 operator *(const VECTOR4 & v0, const MATRIX4 & m0); // return v0 * m0
VECTOR3 operator *(const MATRIX4 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX4 & m0); // return v0 * m0


MATRIX4 inverse(const MATRIX4 & m);  // return inverse of m; return 0 matrix if
// m is singular
MATRIX4 rotate_matrix(int type, float angle); // type: 1:x, 2:y, 3:z
MATRIX4 translate_matrix(float dx, float dy, float dz);
MATRIX4 scale_matrix(float sx, float sy, float sz);

#endif
