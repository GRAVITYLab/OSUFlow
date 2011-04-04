////////////////////////////////////////////////////////////
//
//  Marching Cubes  IsoSurface Class
//
//  Han-Wei Shen   hwshen@nas.nasa.gov
//
////////////////////////////////////////////////////////////

#ifndef VT_ISOSURFACE_H 
#define VT_ISOSURFACE_H 1 

#include "header.h"
#include "VectorMatrix.h"

class vtIsoSurf {

protected: 
  
  float *isov; 
  int nisov; 
  int nisov_allocated; 

  float* field; 
  int xdim, ydim, zdim; 

  void iso_cell(float v, vector<VECTOR3*>& vTriangles); 
  
public:

  vtIsoSurf(float* field, int dx, int dy, int dz); 
  ~vtIsoSurf(); 

  void set_isov(float v); 
  void set_isov(float*, int); 
  
  virtual void execute(vector<VECTOR3*>& vTriangles); 

}; 

#endif










