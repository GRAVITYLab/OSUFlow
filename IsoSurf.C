////////////////////////////////////////////////////////////
//
//  UIsoSurface Class Definitions
//
//  Grid independent isosurface extraction
//  module using FEL. 
//
//  Han-Wei Shen
//  Oct 1, 1997
//
//
////////////////////////////////////////////////////////////

#include "IsoSurf.h"
#include "triangulator.h"

extern int vt_off[9]; 

/////////////////////////////////////////////////////////
//
//   Constructor and Destructor
//
vtIsoSurf::vtIsoSurf(float* infield, int dx, int dy, int dz):
field(infield), xdim(dx), ydim(dy), zdim(dz)

{
	isov = NULL; 
	nisov_allocated = nisov = 0; 

	// too speed up the data access 
	// the external array 
	vt_off[1] = 0; 
	vt_off[2] = 1; 
	vt_off[3] = xdim+1; 
	vt_off[4] = xdim; 
	vt_off[5] = xdim*ydim; 
	vt_off[6] = vt_off[5]+vt_off[2]; 
	vt_off[7] = vt_off[5]+vt_off[3];
	vt_off[8] = vt_off[5]+vt_off[4]; 
}


vtIsoSurf::~vtIsoSurf()
{
}


///////////////////////////////////////////////////
void vtIsoSurf::set_isov(float v) {

	if (nisov_allocated == 0) {
		isov = new float[1];   
		nisov_allocated = 1; 
	}
	nisov = 1; 
	isov[0] = v; 
}

void vtIsoSurf::set_isov(float* vlst, int nv) {
	if (nisov_allocated < nv) {
		if (nisov_allocated != 0) delete[] isov; 
		isov = new float[nv]; 
		nisov_allocated = nv; 
	}
	nisov = nv; 
	for (int i=0; i<nv; i++){
		isov[i] = vlst[i]; 
	}
}

//////////////////////////////////////////
//
//
void vtIsoSurf::execute(vector<VECTOR3*>& vTriangles)
{
	for (int i=0; i<nisov; i++) {
		iso_cell(isov[i], vTriangles); 
	}
}

/////////////////////////////////////////////////////////////
//
// Methods for Regular Marching Cubes Algorithm for regular grids
//
/////////////////////////////////////////////////////////////
//
void vtIsoSurf::iso_cell(float isoval, vector<VECTOR3*>& vTriangles)
{
	int total=0; 
	int index; 

	int dxdy = (xdim-1)*(ydim-1); 
	int dx = (xdim-1); 
	int isocell=0; 

	for (int i=0; i<xdim-1; i++)
		for (int j=0; j<ydim-1; j++) 
			for (int k=0; k<zdim-1; k++) {
				index = k*dxdy + j* dx + i; 
				total+= iso_hex(index, xdim, ydim, zdim, 
					isoval, field, vTriangles); 
			}
			printf("--->find %d triangles.\n", total); 
}
