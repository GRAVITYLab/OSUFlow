
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>
#include <string.h>
#include <limits.h>

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef PNETCDF
#include "pnetcdf.h"
#endif

#include "FileReader.h"

static char **files = NULL; // individual file names of timesteps
static int num_files = 0; // number of timestep files

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
	// MOD-BY-LEETEN 02/07/2011-FROM:
		//	char *filename = new char[100];
	// TO:
  char filenameonly[1024];
  char filename[1024];
	// MOD-BY-LEETEN 02/07/2011-END
  int  totalNum; 
  float* pData = NULL;
  float ** ppData = NULL; 

	// ADD-BY-LEETEN 02/07/2011-BEGIN
	char szPath[1024];
	strcpy(szPath, fname);

	char *szSeparator;
	szSeparator = strrchr(szPath, '/');
	#ifdef WIN32
		if( !szSeparator )
			szSeparator = strrchr(szPath, '\\');
	#endif
	if( NULL == szSeparator )
		szSeparator = &szPath[0];

	*szSeparator = '\0';
	// ADD-BY-LEETEN 02/07/2011-END

  fIn = fopen(fname, "r");
  assert(fIn != NULL);
  fscanf(fIn, "%d", &n_timesteps);

#ifdef DEBUG
  printf(" %d time steps to read. \n", n_timesteps); 
#endif

  ppData = new float*[n_timesteps]; 

  for (int i=0; i<n_timesteps; i++) {
	// MOD-BY-LEETEN 02/07/2011-FROM:
		// fscanf(fIn, "%s", filename);
	// TO:
    fscanf(fIn, "%s", filenameonly);
	if(
		#ifdef WIN32
			filenameonly[1] != ':' &&
			filenameonly[0] != '\\' && 
		#endif
		filenameonly[0] != '/' )
	{
		sprintf(filename, "%s/%s", szPath, filenameonly);
	}
	else
	{
		strcpy(filename, filenameonly);
	}
	// MOD-BY-LEETEN 02/07/2011-END

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
// stores dataset file names when header has already been read
//
void DatasetFiles(char **names, int num_names) {

  files = names;
  num_files = num_names;

}
//---------------------------------------------------------------------------
//
// Read time-varying subdomain data using various data modes
//
// fname: dataset name
// dim: size of total domain
// minB, maxB: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// min_t, max_t: time range of subdomain
// data_mode: RAW, RAW_HEADER, NETCDF
//
float** ReadData(char *fname, float *dim, float *minB, 
		 float *maxB, int min_t, int max_t,
		 DataMode data_mode) { 

  FILE *fIn = NULL;
  FILE *fVecIn;
  char filename[256];
  float* pData = NULL;
  float ** ppData = NULL; 
  int n_timesteps;
  int i, j, k;

  if (files == NULL) {
    assert((fIn = fopen(fname, "r")) != NULL);
    assert(fscanf(fIn, "%d", &n_timesteps) > 0);
  }

#ifdef _MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#ifdef PNETCDF
  char *nc_var_names[3];
  float nc_fill_vals[3];
  char *nc_dim_names[3];
  
  // read the netcdf variable names from the file
  if (data_mode == NETCDF) {

    // dim names of x, y, and z come first line by line
    for (i = 0; i < 3; i++) {
      nc_dim_names[i] = new char[128];
      assert(nc_dim_names[i] != NULL);
      assert(fscanf(fIn, "%s", nc_dim_names[i]) < 128);
    }
    // variable names for u, v, and w come next. a fill value
    // may be specified by the name
    for (i = 0; i < 3; i++) {
      nc_var_names[i] = new char[128];
      assert(nc_var_names[i] != NULL);
      assert(fscanf(fIn, "%s", nc_var_names[i]) < 128);
      nc_fill_vals[i] = NC_FILL_FLOAT;
      assert(fscanf(fIn, "%f", &(nc_fill_vals[i])) > 0);
    }
  }
#endif
#endif

  int numTimesteps = max_t - min_t + 1; 

  ppData = new float*[numTimesteps]; 

  // for all timesteps
  for(i = 0; i < (files == NULL ? n_timesteps : num_files); i++) {

    if (files == NULL)
      assert(fscanf(fIn, "%s", filename) > 0);
    else 
      strcpy(filename, files[i]);

    if (i < min_t || i > max_t) 
      continue; 

#ifdef DEBUG
    printf(" to read %s ...\n", filename); 
#endif

    pData = new float[(int)(maxB[0] - minB[0] + 1) * 
		      (int)(maxB[1] - minB[1] + 1) *
		      (int)(maxB[2] - minB[2] + 1) * 3];


#ifdef _MPI

    if (data_mode == RAW || data_mode == RAW_HEADER) {

      // MPI_IO
      if (0) { // not ready yet
	MPI_File fd;
	assert(MPI_File_open(comm, filename, MPI_MODE_RDONLY,
			     MPI_INFO_NULL, &fd) == MPI_SUCCESS);
	Mpi_ioReadDataRaw(fd, dim, minB, maxB, pData, data_mode);
	MPI_File_close(&fd);
      }

      // Posix
      else {
	assert((fVecIn = fopen(filename, "rb")) != NULL);
	PosixReadDataRaw(fVecIn, dim, minB, maxB, pData, data_mode);
	fclose(fVecIn);
      }	


    } // data_mode == RAW || data_mode == RAW_HEADER


#ifdef PNETCDF

    // Parallel netCDF
    else if (data_mode == NETCDF) {

      int fd;
      assert(ncmpi_open(comm, filename, NC_NOWRITE, MPI_INFO_NULL, 
			&fd) == NC_NOERR);
			
      int nc_dim_ids[3];
      int min_dim_id = INT_MAX;
      for (j = 0; j < 3; j++) {		
	assert(ncmpi_inq_dimid(fd, nc_dim_names[j], 
			       &(nc_dim_ids[j])) == NC_NOERR);
	if (nc_dim_ids[j] < min_dim_id)
	  min_dim_id = nc_dim_ids[j];
      }

      // make the first dim id = 0
      for (j = 0; j < 3; j++)
	nc_dim_ids[j] -= min_dim_id;
			
      MPI_Offset size[3]; // sizes of entire dataset
      MPI_Offset subsize[3]; // sizes of my subdomain
      MPI_Offset start[3]; // starting indices of my subdomain
			
      for (j = 0; j < 3; j++) {
	size[nc_dim_ids[j]] = dim[j];
	start[nc_dim_ids[j]] = (MPI_Offset)minB[j];
	subsize[nc_dim_ids[j]] = (MPI_Offset)(maxB[j] - minB[j] + 1);
      }

      int nflt = subsize[0] * subsize[1] * subsize[2];
      float *var_buf = new float[nflt];
      assert(var_buf != NULL);

      for (j = 0; j < 3; j++) {

	// allow the user to specify null varnames. These will be
	// filled with zeroes
	if (strcmp(nc_var_names[j], "NULL") == 0) {
	  // stride the data buf to zero
	  for (k = 0; k < nflt; k++)
	    pData[k * 3 + j] = 0.0f;
	}

	else {

	  int nc_var_id;
	  // get the variable id
	  assert(ncmpi_inq_varid(fd, nc_var_names[j], 
				 &nc_var_id) == NC_NOERR);
					
	  // collectively read the data
	  assert(ncmpi_get_vara_float_all(fd, nc_var_id,
					  start, subsize, var_buf) == NC_NOERR);

	  // unpack the var buf
	  // handle data already in z y x order (for speed)
	  if (nc_dim_ids[0] == 2 && nc_dim_ids[1] == 1 && nc_dim_ids[2] == 0) {
	    for (k = 0; k < nflt; k++) {
	      if (var_buf[k] == nc_fill_vals[j])
		pData[k * 3 + j] = 0.0;
	      else
		pData[k * 3 + j] = var_buf[k];

	    }
	  }

	  // handle arbitarily ordered data
	  else {
	    int a, b, c;
	    int *abc[3] = {&a, &b, &c};
	    int reorder_k;
	    k = 0;
	    for (a = 0; a < subsize[nc_dim_ids[0]]; a++) {
	      for (b = 0; b < subsize[nc_dim_ids[1]]; b++) {
		for (c = 0; c < subsize[nc_dim_ids[2]]; c++, k++) {
		  reorder_k = 
		    (*(abc[nc_dim_ids[2]]) * subsize[0] * subsize[1])
		    + (*(abc[nc_dim_ids[1]]) * subsize[0]) 
		    + *(abc[nc_dim_ids[0]]);
		  if (var_buf[k] == nc_fill_vals[j])
		    pData[reorder_k * 3 + j] = 0.0;
		  else
		    pData[reorder_k * 3 + j] = var_buf[k];
		}
	      }
	    }
	  } // else

	} // else

      } // for

      // clean up
      delete []var_buf;
      ncmpi_close(fd);

    } // data_mode == NETCDF

#endif // PNETCDF    

    else {
      fprintf(stderr, "Data mode not supported using MPI\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

#endif // MPI

#ifndef _MPI

    // POSIX IO
    assert((fVecIn = fopen(filename, "rb")) != NULL);
    PosixReadDataRaw(fVecIn, dim, minB, maxB, pData, data_mode);
    fclose(fVecIn);

#endif

    // set the data, no matter which io method was used
    ppData[i - min_t] = pData; 

  } // for all timesteps

#ifdef PNETCDF
  if (data_mode == NETCDF) {
    for (i = 0; i < 3; i++) {
      delete [](nc_dim_names[i]);
      delete [](nc_var_names[i]);
    }
  }
#endif

  if (files == NULL)
    fclose(fIn);

  return(ppData); 

}
//-----------------------------------------------------------------------

#ifdef _MPI

//-----------------------------------------------------------------------
//
// Read subdomain data using MPI-IO
//
// f: opened file pointer
// dim: size of total domain
// minB, maxB: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// p: pointer to allocated space for data
// dm: data mode (RAW or RAW_HEADER)
//
void Mpi_ioReadDataRaw(MPI_File f, float *dim, float *minB, 
		       float *maxB, float *p, DataMode dm) {

  int size[3]; // sizes of entire dataset
  int subsize[3]; // sizes of my subdomain
  int start[3]; // starting indices of my subdomain
  MPI_Datatype filetype;
  MPI_Status status;
  int j;

  assert(dm == RAW || dm == RAW_HEADER);

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // sMin and sMax are [x][y][z]
  for (j = 0; j < 3; j++) {
    size[2 - j] = dim[j];
    start[2 - j] = (int)minB[j];
    subsize[2 - j] = (int)(maxB[j] - minB[j] + 1);
  }
				
  // multiply by three since we are reading in raw u,v,w coords	
  size[2] *= 3;
  start[2] *= 3;
  subsize[2] *= 3;
  int nflt = subsize[0] * subsize[1] * subsize[2];
		
  MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C,
			   MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);
  MPI_File_set_view(f, 0, MPI_FLOAT, filetype, (char *)"native",
		    MPI_INFO_NULL);
  assert(MPI_File_read_all(f, p, nflt,
			   MPI_FLOAT, &status) == MPI_SUCCESS);
  assert(status.count == (int)(sizeof(float)) * nflt);

#ifdef BYTE_SWAP
  for (j = 0; j < nflt; j++)
    swap4((char *)&(p[j]));
#endif

  MPI_Type_free(&filetype);

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// Read subdomain data using POSIX
//
// f: opened file pointer
// dim: size of total domain
// minB, maxB: corners of subdomain (inclusive, node-centered)
// ie, the range [0, 10] contains 10 voxels
// p: pointer to allocated space for data
// dm: data mode (RAW or RAW_HEADER)
//
void PosixReadDataRaw(FILE *f, float *dim, float *minB, 
			float *maxB, float *p, DataMode dm) {

  int y, z, size;
  long offset;

  assert(dm == RAW || dm == RAW_HEADER);

  // read the size at the start of the file
  int headeroffset = 0;
  int dimi[3];
  if (dm == RAW_HEADER) {
    assert(fread(dimi, sizeof(int), 3, f) > 0);
#ifdef BYTE_SWAP
    for (int j = 0; j < 3; j++)
      swap4((char *)&(dimi[j]));
#endif
    dim[0] = (float)dimi[0];
    dim[1] = (float)dimi[1];
    dim[2] = (float)dimi[2];
    headeroffset = sizeof(int) * 3;
  }

  // read contiguous rows of vectors
  for (z = (int)minB[2]; z <= (int)maxB[2]; z++) {

    for (y = (int)minB[1]; y <= (int)maxB[1]; y++) {

      offset = (long)(z * dim[0] * dim[1] + y * dim[0] + minB[0]) *
		3 * sizeof(float) + headeroffset;
      size = (int)(maxB[0] - minB[0] + 1) * 3;
      fseek(f, offset, SEEK_SET);
      assert(fread(p, sizeof(float), size, f) > 0);

#ifdef BYTE_SWAP
      for (int j = 0; j < size; j++)
	swap4((char *)&(p[j]));
#endif

      p += size;

    }

  }

}
//-----------------------------------------------------------------------

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
