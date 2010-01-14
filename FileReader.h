#ifdef _MPI
#include <mpi.h>
#endif

float* ReadStaticDataRaw(char* fname, int* dimension); 

float* ReadStaticDataRaw(char *fname, int* dimension, 
			 float* sMin, float* sMax); 

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension); 

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension, 
			       float *minB, float *maxB, 
			       int min_t, int max_t); 

float* ReadStaticDataRaw(char *fname, float *sMin, float* sMax, int* dim);

float** ReadTimeVaryingDataRaw(char *fname, float* sMin, float* sMax, 
			       int* dim, int bt_max, int t_min, int t_max);
float** IndepReadTimeVaryingDataRaw(char *flowName, float* sMin, float* sMax, 
				    int* dim, int bt_max, int t_min, int t_max);
