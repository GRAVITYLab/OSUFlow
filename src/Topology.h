
/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                
// ---------------------------------------------------
//	This is to implement a module to analyze the critical points and its 
//	type. The design and implementation is based on the FAST CFD package 
//	and the document "A tool for visualizing the topology of 3D vector 
//	fields" by A. Globus.											    
//
///////////////////////////////////////////////////////////////////////////////

// The ultimate goal is to first locate the candidate cells that may contain
// critical points (CP), and then find CP within these candidate cells, 
// classifies each critical point using the eigenvalues
// Currently, only implement: given critical points, classify them

#ifndef _VECTOR_TOPOLOGY_H_
#define _VECTOR_TOPOLOGY_H_

#include "header.h"
#include "VectorMatrix.h"
#include "Field.h"

#define ALL_REAL		3
#define TWO_COMPLEX		-3

typedef enum
{
	AttractingNode, AttractingSpiral, AttractingNodeSaddle, AttractingSpiralSaddle,
	RepellingNode, RepellingSpiral, RepellingNodeSaddle, RepellingSpiralSaddle,
	Degenerate3d, No3dCPtype 
}t3dCPtype;

// data structure for critical points
typedef struct
{
	VECTOR3 physicalLocation;
	VECTOR3 compOffset;
	float jacobian[3][3];
	float eigenValues[3];
	int eigenType;
	float eigenVectors[3][3];
	float eigenVectorPL[3][3];	/* physical space location */
	t3dCPtype cp3dType;
}tCriticalPoint;

class CPTopology
{
private:
	CVectorField* m_pField;

public:
	CPTopology();
	CPTopology(CVectorField* pField);
	~CPTopology();

	void CPClassify(tCriticalPoint& cp);
	void makeJacobian(tCriticalPoint& cp);
	void putOddManOutEigenVectorFirst(tCriticalPoint& cp);
	int howManyPositiveEigenValues(tCriticalPoint cp);
	t3dCPtype find3dType(tCriticalPoint cp);
};

/* file polynomials.c */
extern double cube_root ( double x );
extern int solve_cubic ( float, float, float, float, float*, float*, float* );
extern int solve_quadratic ( float, float, float, float*, float* );
extern int solve_linear ( float, float, float* );

/* file eigenvals.c */
extern int compute_eigenvalues( float m[3][3], float eigenvalues[3]);
extern void compute_real_eigenvectors (float m[3][3], float vals[3], float vecs[3][3]);
extern void compute_complex_eigenvectors(float m[3][3], float vals[3], float vecs[3][3]);

#endif

