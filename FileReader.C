
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

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

  lxdim = sMax[0]-sMin[0]+1; 
  lydim = sMax[1]-sMin[1]+1; 
  lzdim = sMax[2]-sMin[2]+1; 

  totalNum = lxdim*lydim*lzdim; 

  pData = new float[totalNum * 3];

  float *p = pData; 
  for (int z = sMin[2]; z<=sMax[2]; z++) {
    for (int y = sMin[1]; y<=sMax[1]; y++) {
      long offset = (z*dimension[0]*dimension[1]+y*dimension[0]+sMin[0])*3*4; 
      fseek(fIn, offset, SEEK_SET); 
      int size = (sMax[0]-sMin[0]+1)*3; 
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
  printf(" %d time steps to read. \n", n_timesteps); 

  ppData = new float*[n_timesteps]; 

  for (int i=0; i<n_timesteps; i++) {
    fscanf(fIn, "%s", filename);
    printf(" to read %s ...\n", filename); 
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
				     int min_t, int max_t)
{
  FILE *fIn;
  FILE *fVecIn;
  char *filename = new char[100];
  int  totalNum; 
  float* pData = NULL;
  float ** ppData = NULL; 
  int lxdim, lydim, lzdim; 

  fIn = fopen(fname, "r");
  assert(fIn != NULL);
  fscanf(fIn, "%d", &n_timesteps);
  printf(" %d time steps to read. \n", n_timesteps); 

  if (min_t == -1 || max_t == -1)  {
    min_t = 0; 
    max_t = n_timesteps-1; 
  }

  int numTimesteps = max_t-min_t+1; 

  ppData = new float*[numTimesteps]; 

  lxdim = maxB[0]-minB[0]+1; 
  lydim = maxB[1]-minB[1]+1; 
  lzdim = maxB[2]-minB[2]+1; 

  printf(" min max t = %d %d \n", min_t, max_t); 

  for(int iFor = 0; iFor < n_timesteps; iFor++)
    {
      fscanf(fIn, "%s", filename);
      if (iFor <min_t || iFor >max_t) 
	continue; 
      printf(" to read %s ...\n", filename); 
      fVecIn = fopen(filename, "rb");
      if (fVecIn==NULL) {
	printf(" problem opening %s. skip\n", filename); 
	continue; 
      }
      fread(dimension, sizeof(int), 3, fVecIn);
      totalNum = lxdim*lydim*lzdim; 
      pData = new float[totalNum * 3];

      float *p = pData; 
      for (int z = minB[2]; z<=maxB[2]; z++) {
	for (int y = minB[1]; y<=maxB[1]; y++) {
	  long offset = (z*dimension[0]*dimension[1]+y*dimension[0]+minB[0])*3*4;
	  fseek(fVecIn, offset, SEEK_SET);
	  int size = (maxB[0]-minB[0]+1)*3;
	  fread(p, sizeof(float), size, fVecIn);
	  p+=size;
	}
      }
      fclose(fVecIn);
      ppData[iFor-min_t] = pData; 
    }
  return(ppData); 
}
