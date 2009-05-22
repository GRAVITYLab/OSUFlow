
float* ReadStaticDataRaw(char* fname, int* dimension); 

float* ReadStaticDataRaw(char *fname, int* dimension, 
		      float* sMin, float* sMax); 

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			    int *dimension); 

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			    int *dimension, 
			    float *minB, float *maxB, 
			    int min_t, int max_t); 
