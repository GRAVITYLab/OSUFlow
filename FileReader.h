#ifdef _MPI
#include <mpi.h>
#endif

float* ReadStaticDataRaw(char* fname, int* dimension); 

float* ReadStaticDataRaw(char *fname, int* dimension, 
			 float* sMin, float* sMax); 

float* ReadStaticDataRaw(char *fname, float *sMin, float* sMax, int* dim);

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension); 

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension, 
			       float *minB, float *maxB, 
			       int min_t, int max_t); 

float** ReadTimeVaryingDataRaw(char *fname, float* sMin, float* sMax, 
			       int* dim, int bt_max, int t_min, int t_max);

float** ReadTimeVaryingDataRaw(char *fname, float *dim, float *minB, 
			       float *maxB, int min_t, int max_t,
			       bool read_dims);
void swap4(char *n);
