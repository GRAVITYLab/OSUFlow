#ifndef _POLYNOMIALS_H_
#define _POLYNOMIALS_H_

#include <stdio.h>
#include <stdlib.h>
#include <draw_topo.h>
#include <fields.h>
#include "Data.h"

void Tinit( void );
void TmakeTopology( void );
void TsetGlyphSize( float );
void TfindCriticalPoints( void );
void TsurfaceTopology( void );
void TsurfaceVortexCores( void );
void TvortexCores( void );
void TnodesOnly( void );
void Tall( void );
void ToffSurface( void );
void Ttrace( void );

extern tTopology theTopology;	/* UGLY: shouldn't be global */
extern tGraphicObject theGraphicObject;	/* UGLY: shouldn't be global */

/* file ipc.c */

extern void IPCcheck( void );
extern void IPCinitialize( int argc, char* argv[], char* module );

/* file main.c */

extern void MustDraw( void );
extern void setSelectionUI( void );

/* file polynomials.c */

extern double cube_root ( double x );
extern int solve_cubic ( float, float, float, float, float*, float*, float* );
extern int solve_quadratic ( float, float, float, float*, float* );
extern int solve_linear ( float, float, float* );

/* file phys_to_comp.c */

extern void phys_to_comp( tData data );


/* file initialize.c */

extern void initialize( int x, int y );


/* file panels.c */

extern void initialize_panels( int x, int y );

/* file eigenvals.c */

extern int compute_eigenvalues( float m[3][3], float eigenvalues[3] );

/* file eigenvecs.c */
extern void compute_real_eigenvectors (float m[3][3], float vals[3], float vecs[3][3]);
extern void compute_complex_eigenvectors( float m[3][3], float vals[3], float vecs[3][3] );

/* file inv3x3.f */

extern void inv3x3( float jac[3][3], float jacinvf[3][3], int* status);
extern void inv3x3_( float jac[3][3], float jacinvf[3][3], int* status);


#endif



