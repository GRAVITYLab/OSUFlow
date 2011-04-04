#ifdef _MPI
#include <mpi.h>
#endif

#ifndef _FILEREADER_H
#define _FILEREADER_H

typedef enum {
  RAW = 0,
  RAW_HEADER,
  NETCDF,
  HDF_FLOAT,
  HDF_DOUBLE,
} DataMode;

void DatasetFiles(char **dataset_files, int num_dataset_files);

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

float** ReadData(char *fname, float *dim, float *minB, 
		 float *maxB, int min_t, int max_t,
		 DataMode data_mode);

#ifdef _MPI
void Mpi_ioReadDataRaw(MPI_File f, float *dim, float *minB, 
		       float *maxB, float *p, DataMode dm);
#endif

void PosixReadDataRaw(FILE *f, float *dim, float *minB, 
		      float *maxB, float *p, DataMode dm);
void swap4(char *n);

#endif
