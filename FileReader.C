
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

#ifdef _MPI
#include <mpi.h>
#endif

#include "FileReader.h"

/////////////////////////////////////////////////////////

float* ReadStaticDataRaw(char* fname, int* dimension) 
{
 
  FILE * fIn; 
  int totalNum;
  float *pData; 

  fIn = fopen(fname, "rb");
  assert(fIn != NULL);
  fread(dimension, sizeof(int), 3, fIn);
  totalNum = dimension[0] * dimension[1] * dimension[2];
  pData = new float[totalNum * 3];
  fread(pData, sizeof(float), totalNum*3, fIn);
  fclose(fIn);
  return(pData); 
}


//////////////////////////////////////////////////////////
//
//  Read a subset of the data file
//
float* ReadStaticDataRaw(char *fname, int* dimension, 
			       float* sMin, float* sMax) 
{
 
  FILE * fIn; 
  int totalNum;
  float *pData; 
  int lxdim, lydim, lzdim; 
  
  fIn = fopen(fname, "rb");
  assert(fIn != NULL);
  fread(dimension, sizeof(int), 3, fIn);

  lxdim = (int)(sMax[0]-sMin[0]+1); 
  lydim = (int)(sMax[1]-sMin[1]+1); 
  lzdim = (int)(sMax[2]-sMin[2]+1); 

  totalNum = lxdim*lydim*lzdim; 

  pData = new float[totalNum * 3];

  float *p = pData; 
  for (int z = (int)sMin[2]; z<=(int)sMax[2]; z++) {
    for (int y = (int)sMin[1]; y<=(int)sMax[1]; y++) {
      long offset = (long)(z*dimension[0]*dimension[1]+y*dimension[0]+sMin[0])*3*4;
      fseek(fIn, offset, SEEK_SET); 
      int size = (int)(sMax[0]-sMin[0]+1)*3; 
      fread(p, sizeof(float), size, fIn); 
      p+=size; 
    }
  }
  fclose(fIn);
  return(pData); 
}

///////////////////////////////////////////////////////////////

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
				     int *dimension)
{
  FILE *fIn;
  FILE *fVecIn;
  char *filename = new char[100];
  int  totalNum; 
  float* pData = NULL;
  float ** ppData = NULL; 

  fIn = fopen(fname, "r");
  assert(fIn != NULL);
  fscanf(fIn, "%d", &n_timesteps);

#ifdef DEBUG
  printf(" %d time steps to read. \n", n_timesteps); 
#endif

  ppData = new float*[n_timesteps]; 

  for (int i=0; i<n_timesteps; i++) {
    fscanf(fIn, "%s", filename);

#ifdef DEBUG
    printf(" to read %s ...\n", filename); 
#endif

    fVecIn = fopen(filename, "rb");
    if (fVecIn==NULL) {
      printf(" problem opening %s. skip\n", filename); 
      continue; 
    }
    fread(dimension, sizeof(int), 3, fVecIn);
    totalNum = dimension[0] * dimension[1] * dimension[2];
    pData = new float[totalNum * 3];
    fread(pData, sizeof(float), totalNum*3, fVecIn);
    ppData[i] = pData; 
    fclose(fVecIn);
  }
  return(ppData); 
}

//////////////////////////////////////////////////////////

float** ReadTimeVaryingDataRaw(char *fname, int& n_timesteps, 
			       int *dimension, 
			       float *minB, float *maxB, 
			       int min_t, int max_t) {

  FILE *fIn;
  FILE *fVecIn;
  char *filename = new char[100];
  int  totalNum; 
  float* pData = NULL;
  float ** ppData = NULL; 
  int lxdim, lydim, lzdim; 
  int i;

  fIn = fopen(fname, "r");
  assert(fIn != NULL);
  fscanf(fIn, "%d", &n_timesteps);

#ifdef DEBUG
  printf(" %d time steps to read. \n", n_timesteps); 
#endif

  if (min_t == -1 || max_t == -1)  {
    min_t = 0; 
    max_t = n_timesteps-1; 
  }

  int numTimesteps = max_t-min_t+1; 

  ppData = new float*[numTimesteps]; 

  lxdim = (int)(maxB[0]-minB[0]+1); 
  lydim = (int)(maxB[1]-minB[1]+1); 
  lzdim = (int)(maxB[2]-minB[2]+1); 

#ifdef DEBUG
  printf(" min max t = %d %d \n", min_t, max_t); 
#endif

  // for all timesteps
  for(int iFor = 0; iFor < n_timesteps; iFor++) {

    fscanf(fIn, "%s", filename);
    if (iFor <min_t || iFor >max_t) 
      continue; 

#ifdef DEBUG
    printf(" to read %s ...\n", filename); 
#endif

    fVecIn = fopen(filename, "rb");
    if (fVecIn==NULL) {
      printf(" problem opening %s. skip\n", filename); 
      continue; 
    }

    // read the size at the start of the file
    fread(dimension, sizeof(int), 3, fVecIn);
    // swap bytes
#ifdef BYTE_SWAP
    for (i = 0; i < 3; i++)
      swap4((char *)&(dimension[i]));
#endif

    totalNum = lxdim*lydim*lzdim; 
    pData = new float[totalNum * 3];

    float *p = pData; 
    for (int z = (int)minB[2]; z<=(int)maxB[2]; z++) {

      for (int y = (int)minB[1]; y<=(int)maxB[1]; y++) {
	long offset = (long)(z*dimension[0]*dimension[1]+y*dimension[0]+
			     minB[0])*3*4;
	fseek(fVecIn, offset, SEEK_SET);
	int size = (int)(maxB[0]-minB[0]+1)*3;
	fread(p, sizeof(float), size, fVecIn);

	// swap bytes
#ifdef BYTE_SWAP
	for (i = 0; i < size; i++)
	  swap4((char *)&(p[i]));
#endif
	p+=size;
      }

    }

    fclose(fVecIn);
    ppData[iFor-min_t] = pData; 

  } // for all timesteps

  fclose(fIn);
  return(ppData); 

}
//---------------------------------------------------------------------------
//
// ReadTimeVaryingDataRaw
//
// Read time-varying subdomain data using posix I/O
//
// fname: dataset name
// dim: size of total domain
// minB, maxB: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// t_min, t_max: time range of subdomain
// read_dims: whether to read dimensions from the start of the file
//            (if so, dim will contain the dimensions read and can
//            be passed in uninitialized)
//
float** ReadTimeVaryingDataRaw(char *fname, float *dim, float *minB, 
			       float *maxB, int min_t, int max_t,
			       bool read_dims) {

  FILE *fIn;
  FILE *fVecIn;
  char *filename = new char[256];
  float* pData = NULL;
  float ** ppData = NULL; 
  int n_timesteps;
  int i, j;

  assert((fIn = fopen(fname, "r")) != NULL);
  fscanf(fIn, "%d", &n_timesteps);

  if (min_t == -1 || max_t == -1)  {
    min_t = 0; 
    max_t = n_timesteps - 1; 
  }

  int numTimesteps = max_t - min_t + 1; 

  ppData = new float*[numTimesteps]; 

  // for all timesteps
  for(i = 0; i < n_timesteps; i++) {

    fscanf(fIn, "%s", filename);

    if (i < min_t || i > max_t) 
      continue; 

#ifdef DEBUG
    printf(" to read %s ...\n", filename); 
#endif

    assert((fVecIn = fopen(filename, "rb")) != NULL);

    // read the size at the start of the file
    if (read_dims) {
      fread(dim, sizeof(int), 3, fVecIn);
#ifdef BYTE_SWAP
      for (j = 0; j < 3; j++)
	swap4((char *)&(dim[j]));
#endif
    }

    pData = new float[(int)(maxB[0] - minB[0]) * (int)(maxB[1] - minB[1]) *
		      (int)(maxB[2] - minB[2]) * 3];
    float *p = pData; 

    // read contiguous rows of vectors
    for (int z = (int)minB[2]; z < (int)maxB[2]; z++) {

      for (int y = (int)minB[1]; y < (int)maxB[1]; y++) {

	long offset = (long)(z * dim[0] * dim[1] + y * dim[0] + minB[0]) *
	  3 * sizeof(float);
	int size = (int)(maxB[0] - minB[0]) * 3;

	fseek(fVecIn, offset, SEEK_SET);
	fread(p, sizeof(float), size, fVecIn);

#ifdef BYTE_SWAP
	for (j = 0; j < size; j++)
	  swap4((char *)&(p[j]));
#endif
	p += size;

      }

    }

    fclose(fVecIn);
    ppData[i - min_t] = pData; 

  } // for all timesteps

  fclose(fIn);
  return(ppData); 

}
//-----------------------------------------------------------------------

// MPI functions

#ifdef _MPI

//-----------------------------------------------------------------------
//
// ReadTimeVaryingDataRaw
//
// Read time-varying subdomain data collectively
// used MPI-IO collectives to perform the I/O
//
// flowName: dataset name
// sMin, sMax: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// dim: size of total domain
// bt_max: max number of time steps in any block
// t_min, t_max: time range of subdomain
//
float** ReadTimeVaryingDataRaw(char *flowName, float* sMin, float* sMax, 
			       int* dim, int bt_max, int t_min, int t_max) {

  float **Data; // the data
  MPI_File fd;
  MPI_Datatype filetype;
  int size[3]; // sizes of entire dataset
  int subsize[3]; // sizes of my subdomain
  int start[3]; // starting indices of my subdomain
  MPI_Status status;
  int rank;
  FILE *meta; // metadata file with file names of time step files
  int timesteps; // number of timesteps in the entire dataset
  char filename[256]; // individual timestep file name
  int nflt; // number of floats in my spatial domain (3 * number of points)
  int nt; // number of time steps in my time-space domain
  int null_reads; // number of null reads to fill out max number of time steps
  int i, t;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &rank);

  // read metadata
  meta = fopen(flowName, "r");
  assert(meta != NULL);
  fscanf(meta, "%d", &timesteps);

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // sMin and sMax are [x][y][z]
  for (i = 0; i < 3; i++) {
    size[2 - i] = dim[i];
    start[2 - i] = (int)sMin[i];
    subsize[2 - i] = (int)(sMax[i] - sMin[i]);
  }
  size[2] *= 3;
  start[2] *= 3;
  subsize[2] *= 3;
  nflt = subsize[0] * subsize[1] * subsize[2];

  // allocate data space
  nt = t_max - t_min + 1;
  assert((Data = new float*[nt]) != NULL);
  for (i = 0; i < nt; i++)
    assert((Data[i] = new float[nflt]) != NULL);

  // for all time steps
  for (t = 0; t < timesteps; t++) {

    fscanf(meta, "%s", filename);
    if (t < t_min || t > t_max)
      continue;

    // open the file
    assert(MPI_File_open(comm, filename, MPI_MODE_RDONLY,
			 MPI_INFO_NULL, &fd) == MPI_SUCCESS);

    // do the actual collective io
    MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
			     MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(fd, 0, MPI_FLOAT, filetype, (char *)"native", 
		      MPI_INFO_NULL);
    assert(MPI_File_read_all(fd, Data[t - t_min], nflt,
			     MPI_FLOAT, &status) == MPI_SUCCESS);
    assert(status.count == sizeof(float) * nflt);

    // clean up
    MPI_File_close(&fd);
    MPI_Type_free(&filetype);

    //swap bytes
#ifdef BYTE_SWAP
    for (i = 0; i < 3 * nflt; i++)
      swap4((char *)&(Data[i]));
#endif

  } // for all time steps

  // null reads to fill out the max number of time steps
  null_reads = bt_max - t_max + t_min - 1;
  assert(null_reads >= 0);
  if (null_reads)
    assert(MPI_File_open(comm, filename, MPI_MODE_RDONLY,
			 MPI_INFO_NULL, &fd) == MPI_SUCCESS);

  // for any null reads needed
  for (t = 0; t < null_reads; t++) {
    subsize[0] = subsize[1] = subsize[2] = 1;
    start[0] = start[1] = start[2] = 0;
    MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
			     MPI_FLOAT, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(fd, 0, MPI_FLOAT, filetype, (char *)"native", 
		      MPI_INFO_NULL);
    assert(MPI_File_read_all(fd, Data, 0, MPI_FLOAT, &status) == MPI_SUCCESS);
    MPI_Type_free(&filetype);

  }

  // clean up
  if (null_reads)
    MPI_File_close(&fd);
  fclose(meta);

  return Data;

}
//--------------------------------------------------------------------------

#endif

//---------------------------------------------------------------------------

// utility functions

//---------------------------------------------------------------------------
//
// swap4(n)
//
// Swaps 4 bytes from 1-2-3-4 to 4-3-2-1 order.
// cast the input as a char and use on any 4 byte variable
//
void swap4(char *n) {

  char *n1;
  char c;

  n1 = n + 3;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//----------------------------------------------------------------------------
