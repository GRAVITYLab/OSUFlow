#ifdef HDF5

//------------------------------------------------------------------------------
//
// FLASH I/O library
//
// Used with permission from Mike Papka
// Modified by Tom Peterka
//
//------------------------------------------------------------------------------

// General Headers
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

// use HDF5 version 1.6 API
#define H5_USE_16_API

// Class Headers
#include <flashhdf5_float.h>
// extern "C" 
// {
// #include<H5Ppublic.h>
// }

// Utility Functions
#define CONSTRUCTOR 1
#define LEAF_NODE 1
#define HDF_RETURN_CHECK(i) {if(i == -1) return 0;else return 1;}

void
HDF_ERROR_CHECK(int errorflag, char *message)
{
    // HDF5 Returns a negative number on failure
    if(errorflag < 0)
    {
	cout << "ERROR: " << message << endl;
	exit(0);
    }// End of IF
}// End of HDF_ERROR_CHECK()
 
// Private Functions

void
FlashHDFFile::_ResetSettings(int constructor = 0)
{
    int i;

    numberOfLeafs = 0;
    numberOfDimensions = 0; // Number of dimensions in the dataset
    numberOfBlocks = 0;
    if(!constructor)
    {
	if(min != NULL)
	    delete [] min;
	if(max != NULL)
	    delete [] max;
    }// End of IF
    min = max = NULL;
    for(i=0;i<3;i++)
	cellDimensions[i] = -1;
    if (node_type != NULL)
      delete [] node_type;
    node_type = NULL;
    
}// End of _ResetSettings()

void
FlashHDFFile::_InitSettings()
{
    
    numberOfBlocks = GetNumberOfBlocks();
    //    cout << "Blocks finished " << endl;
    
    numberOfLeafs = GetNumberOfLeafs();
    //    cout << "leaves finished" << endl;
    
    numberOfDimensions = GetNumberOfDimensions();
    //    cout << "dims finished " << endl;
    
    _SetCellDimensions();
    //    cout << "leaving init" << endl;
    
}// End of _initSettings(

// FlashHDFFile Public Class Functions

FlashHDFFile::FlashHDFFile()
{
    _ResetSettings(CONSTRUCTOR);
}// End of FlashHDFFile()

#ifdef _MPI

//-----------------------------------------------------------------------
//
// Constructor for collective I/O
//
// filename: dataset file name
// comm: MPI communicator
//
// Tom Peterka, 11/2/09
//
FlashHDFFile::FlashHDFFile(char *filename, float scale, MPI_Comm comm) {

  MPI_Info info;

  node_type = NULL;
  _ResetSettings(CONSTRUCTOR);

  this->scale = scale;

  // open the file w/ collective I/O
  MPI_Info_create(&info); // todo: figure out appropriate hints for hdf5
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  assert((datasetId = H5Fopen(filename, H5F_ACC_RDONLY, plist_id)) >= 0);

  // create the proplist and set it for collective I/O
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  _InitSettings();

}
//-----------------------------------------------------------------------

#endif

//-----------------------------------------------------------------------
//
// Constructor for independent I/O
//
// filename: dataset file name
//
// Tom Peterka, 11/2/09
//
FlashHDFFile::FlashHDFFile(char *filename, float scale) {

  node_type = NULL;
  _ResetSettings(CONSTRUCTOR);

  this->scale = scale;

  assert((datasetId = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT)) >= 0);

  // create the proplist
  plist_id = H5Pcreate(H5P_DATASET_XFER);

  _InitSettings();

}
//-----------------------------------------------------------------------

#define FILE_SPLITTING 0

int
FlashHDFFile::Open(char *filename)
{
  _ResetSettings(CONSTRUCTOR);
  /* set the file access template for parallel IO access */
  hid_t acc_template = H5Pcreate(H5P_FILE_ACCESS);
  
  if(strchr(filename,'%'))
    {
      fprintf(stderr,"\nOpening a split file\n");
      hid_t gacc_template = H5Pcreate(H5P_FILE_ACCESS);
      hsize_t ms = 0;
      H5Pset_fapl_family(gacc_template,ms,acc_template);
      datasetId = H5Fopen(filename,H5F_ACC_RDONLY,gacc_template);
    }/* End of IF */
  else
    {
      fprintf(stderr,"\nOpening a single file\n");
      datasetId = H5Fopen(filename,H5F_ACC_RDONLY,acc_template);
    }/* End of ELSE */
  
  HDF_ERROR_CHECK(datasetId,(char *)"ERROR: Failed On File Open");

  _InitSettings();

  return 1;
    
}// End of Open(char *filename)

int
FlashHDFFile::Close()
{
    _ResetSettings();
    H5Fclose(datasetId);
    return 1;
}// End of Close()

////
int 
FlashHDFFile::GetNumberOfDimensions()
{ 

  typedef struct int_list_t {
  char name[80];
  int value;
} int_list_t;
                     
hsize_t dimens_1d, maxdimens_1d;
hid_t int_list_type;
hid_t string_type;
int_list_t *int_list;
char *string_index;
char testName[80];
char aname[80];
int len;
        
    if(numberOfDimensions <= 0) {
        
        string_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(string_type, MAX_STRING_LENGTH);
        hid_t dataset = H5Dopen(datasetId, "integer scalars");
        hid_t dataspace = H5Dget_space(dataset);
        /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
        H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
        /* malloc a pointer to a list of int_list_t's */
        int_list = new int_list_t[dimens_1d];
        /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
        hid_t memspace = H5Screate_simple(1, &dimens_1d, NULL);
        /* create an empty vessel sized to hold one int_list_t's worth of data */
        int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
        /* subdivide the empty vessel into its component sections (name and value) */
        H5Tinsert(int_list_type,
                "name",
                HOFFSET(int_list_t, name),
                string_type);
        H5Tinsert(int_list_type,
                "value", 
                HOFFSET(int_list_t, value),
                H5T_NATIVE_INT);
        /* read the data into 'int_list' */
        status = H5Dread(dataset, int_list_type, memspace, dataspace,
                        H5P_DEFAULT, int_list);
        for(unsigned int i=0; i<dimens_1d; i++) {

                strcpy(testName,int_list[i].name);
                string_index = testName;
                len = 0;

                while(*string_index != ' ' && *string_index != '\0') {
                        aname[len] = testName[len];
                        len++;
                        string_index++;
                        }
                *(aname+len) = '\0';

                if(strcmp(aname,"nxb") == 0) {
                        cellDimensions[0] = int_list[i].value;
                        }
                if(strcmp(aname,"nyb") == 0) {
                        cellDimensions[1] = int_list[i].value;
                        }
                if(strcmp(aname,"nzb") == 0) {
                        cellDimensions[2] = int_list[i].value;
                        }
                }
        delete [] int_list;
        H5Tclose(int_list_type);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);

        for(int i=0;i<3;i++)
        {
            numberOfDimensions = ((cellDimensions[i] > 0) ? (numberOfDimensions + 1) : numberOfDimensions);
        }// End of FOR

    }// End of IF

    return numberOfDimensions;

}// End of GetNumberOfDimensions()

//-----------------------------------------------------------------------
//
// GetNumBlocks
//
// reads the total number of blocks in the file
//
// Tom Peterka, 11/6/09
//
int FlashHDFFile::GetNumBlocks() {

  typedef struct int_list_t {
    char name[80];
    int value;
  } int_list_t;

  hsize_t dimens_1d, maxdimens_1d;
  hid_t int_list_type;
  hid_t string_type;
  int_list_t *int_list;
  char testName[80];
  char *string_index;
  char aname[80];
  int len;

  if(numberOfBlocks <= 0) {
                  
    string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, MAX_STRING_LENGTH);
    hid_t dataset = H5Dopen(datasetId, "integer scalars");
    hid_t dataspace = H5Dget_space(dataset);
    /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
    H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
    /* malloc a pointer to a list of int_list_t's */
    int_list = new int_list_t[dimens_1d];
    /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
    hid_t memspace = H5Screate_simple(1, &dimens_1d, NULL);
    /* create an empty vessel sized to hold one int_list_t's worth of data */
    int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
    /* subdivide the empty vessel into its component sections (name and value) */
    H5Tinsert(int_list_type,
	      "name", 
	      HOFFSET(int_list_t, name),
	      string_type);
    H5Tinsert(int_list_type,
	      "value",
	      HOFFSET(int_list_t, value),
	      H5T_NATIVE_INT);
    /* read the data into 'int_list' */
    status = H5Dread(dataset, int_list_type, memspace, dataspace,
		     H5P_DEFAULT, int_list);
    for(unsigned int i=0; i<dimens_1d; i++) {
      strcpy(testName,int_list[i].name);
      string_index = testName;

      len = 0;
      while(*string_index != ' ' && *string_index != '\0') {
	aname[len] = testName[len];
	len++;
	string_index++;
      }
      *(aname+len) = '\0';
      if(strcmp(aname,"globalnumblocks") == 0) {
	numberOfBlocks = int_list[i].value;
      }
    }
    delete [] int_list;
    H5Tclose(int_list_type);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);


  }

  return numberOfBlocks;

}
//-----------------------------------------------------------------------

int
FlashHDFFile::GetNumberOfBlocks()
{

  typedef struct int_list_t {
  char name[80];
  int value;
} int_list_t;

hsize_t dimens_1d, maxdimens_1d;
hid_t int_list_type;
hid_t string_type;
int_list_t *int_list;
char testName[80];
char *string_index;
char aname[80];
int len;

    if(numberOfBlocks <= 0) {
                  
        string_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(string_type, MAX_STRING_LENGTH);
        hid_t dataset = H5Dopen(datasetId, "integer scalars");
        hid_t dataspace = H5Dget_space(dataset);
        /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
        H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
        /* malloc a pointer to a list of int_list_t's */
        int_list = new int_list_t[dimens_1d];
        /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
        hid_t memspace = H5Screate_simple(1, &dimens_1d, NULL);
        /* create an empty vessel sized to hold one int_list_t's worth of data */
        int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
        /* subdivide the empty vessel into its component sections (name and value) */
        H5Tinsert(int_list_type,
                "name", 
                HOFFSET(int_list_t, name),
                string_type);
        H5Tinsert(int_list_type,
                "value",
                HOFFSET(int_list_t, value),
                H5T_NATIVE_INT);
        /* read the data into 'int_list' */
        status = H5Dread(dataset, int_list_type, memspace, dataspace,
                        H5P_DEFAULT, int_list);
        for(unsigned int i=0; i<dimens_1d; i++) {
                strcpy(testName,int_list[i].name);
                string_index = testName;

                len = 0;
                while(*string_index != ' ' && *string_index != '\0') {
                        aname[len] = testName[len];
                        len++;
                        string_index++;
                        }
                *(aname+len) = '\0';
                if(strcmp(aname,"globalnumblocks") == 0) {
                        numberOfBlocks = int_list[i].value;
                        }
                }
        delete [] int_list;
        H5Tclose(int_list_type);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);


    }// End of IF

    return numberOfBlocks;
}// End of GetNumberOfBlocks



int
FlashHDFFile::GetNumberOfLeafs()
{
    if(numberOfLeafs <= 0)
    {
	numberOfLeafs = 0;
	if(numberOfBlocks <= 0)
	    GetNumberOfBlocks();
	int i;
	for(i=0;i<numberOfBlocks;i++)
	    if(GetNodeType(i) == LEAF_NODE)
		numberOfLeafs++;
    }// End of IF
    
    return numberOfLeafs;
}// End of GetNumberOfLeafs()

int
FlashHDFFile::GetNumberOfGlobalIds()
{
    if(numberOfDimensions <= 0)
    {
	this->GetNumberOfDimensions();
	return ((2*numberOfDimensions) + (int)powf(2.0, (float)numberOfDimensions) + 1);
    }// End of IF
    else
	return ((2*numberOfDimensions) + (int)powf(2.0, (float)numberOfDimensions) + 1);
    
}// End of GetNumberOfGlobalIds()

int
FlashHDFFile::GetGlobalIds(int idx,int gids[])
{
    hid_t dataspace, dataset, memspace;
    int rank;
    hsize_t dimens_1d, dimens_2d[2];
    hid_t ierr;
    
    hsize_t start_2d[2];
    hsize_t stride_2d[2], count_2d[2];
    
    herr_t status;
    
    unsigned int i;
    for(i=0;i<15;i++)
	gids[i] = -2;

    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 15; // WRONG Hardcoded value;
    
    /* define the dataspace -- as described above */
    start_2d[0]  = (hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = 1;
    count_2d[1]  = 15; // WRONG Hardcoded value;
    
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_2d, stride_2d, count_2d, NULL);

    /* define the memory space */
    rank = 1;
    dimens_1d = 15; // WRONG Hardcoded value;
    memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    
    dataset = H5Dopen(datasetId, "gid");
    status  = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
		      H5P_DEFAULT, gids);
    

    for(i=0;i<count_2d[1];i++)
    {
	if(gids[i] >= 1)
	    gids[i] = gids[i]-1;
    }// End of FOR

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 1;
}// End of GetGlobalIds(int idx,int gids[])

// Gets the dimensions of an individual block, for viz its usually 9x9x9
//
void
FlashHDFFile::GetCellDimensions(int dimensions[3])
{
typedef struct int_list_t {
  char name[80];
  int value;
} int_list_t;

hsize_t dimens_1d, maxdimens_1d;
hid_t int_list_type;
hid_t string_type;
int_list_t *int_list;
char *string_index;
char testName[80];
char aname[80];
int len;

        string_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(string_type, MAX_STRING_LENGTH);
        hid_t dataset = H5Dopen(datasetId, "integer scalars");
        hid_t dataspace = H5Dget_space(dataset);
        /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
        H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
        /* malloc a pointer to a list of int_list_t's */
        int_list = new int_list_t[dimens_1d];
        /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
        hid_t memspace = H5Screate_simple(1, &dimens_1d, NULL);
        /* create an empty vessel sized to hold one int_list_t's worth of data */
        int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));
        /* subdivide the empty vessel into its component sections (name and value) */
        H5Tinsert(int_list_type,
                "name",
                HOFFSET(int_list_t, name),
                string_type);
        H5Tinsert(int_list_type,
                "value",
                HOFFSET(int_list_t, value),
                H5T_NATIVE_INT);
        /* read the data into 'int_list' */
        status = H5Dread(dataset, int_list_type, memspace, dataspace,
                        H5P_DEFAULT, int_list);
        for(unsigned int i=0; i<dimens_1d; i++) {

                strcpy(testName,int_list[i].name);
                string_index = testName;
                len = 0;

                while(*string_index != ' ' && *string_index != '\0') {
                        aname[len] = testName[len];
                        len++;
                        string_index++;
                        }
                *(aname+len) = '\0';

                if(strcmp(aname,"nxb") == 0) {
                        cellDimensions[0] = int_list[i].value;
                        }
                if(strcmp(aname,"nyb") == 0) {
                        cellDimensions[1] = int_list[i].value;
                        }
                if(strcmp(aname,"nzb") == 0) {
                        cellDimensions[2] = int_list[i].value;
                        }
                }
        delete [] int_list;
        H5Tclose(int_list_type);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        H5Dclose(dataset);

        dimensions[0] = cellDimensions[0];
        dimensions[1] = cellDimensions[1];
        dimensions[2] = cellDimensions[2];

    return;
}// End of GetNumberOfCellDimensions()

// This file is out of data 8/15/01 needs to be updated
void
FlashHDFFile::_SetCellDimensions()
{
    if((cellDimensions[0] == -1) || (cellDimensions[1] == -1)
       || (cellDimensions[2] == -1))
    {
	int dimsizes[3];
	// I don't think this should be hardcoded
	hid_t dataset = H5Dopen(datasetId, "number of zones per block");
	
	status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
			 H5P_DEFAULT, dimsizes);
	
	for(int i=0;i<3;i++)
	{
	    cellDimensions[i] = dimsizes[i];
	}// End of FOR
    }// End of IF

    return;
}

//----------------------------------------------------------------------------
//
// GetNodeType
// returns the node type of a block, 1 = leaf node
//
// reads the file only upon the first call and stores all node types in memory
// after the first time called reads the node type from a buffer in memory
//
// Tom Peterka, 11/7/09
//
int FlashHDFFile::GetNodeType(int idx) {

  if (node_type == NULL) {
    node_type = new int[numberOfBlocks];
    GetIntVecBlocks((char *)"node type", 0, numberOfBlocks, 1, 0, node_type);
  }

  return node_type[idx];

}
//----------------------------------------------------------------------------
int FlashHDFFile::GetNodeType(int idx,int runsize,int *nodearray) {
  
  int rank;
  hsize_t dimens_1d;
  hid_t ierr;
  hid_t dataspace, memspace, dataset;
  hsize_t start_1d;
  hsize_t stride_1d, count_1d;
  herr_t status;

  rank = 1;
  dimens_1d = numberOfBlocks;

  /* define the dataspace -- same as above */
  start_1d  = (hsize_t)idx;
  stride_1d = 1;
  count_1d  = runsize;
    
  dataspace = H5Screate_simple(rank, &dimens_1d, NULL);

  ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			     &start_1d, &stride_1d, &count_1d, NULL);
    
  /* define the memory space */
  dimens_1d = runsize;
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    
  /* create the dataset from scratch only if this is our first time */
  dataset = H5Dopen(datasetId, "node type");
  status  = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
		    H5P_DEFAULT, (void *)nodearray);

    
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  int leafcount = 0;
  for(int i=0;i<runsize;i++) {
    if(nodearray[i] == LEAF_NODE)
      leafcount++;
  }
    
 
  return leafcount;

}

int
FlashHDFFile::GetRefinementLevel(int idx)
{
    int refinelevel;
    
    int rank;
    hsize_t dimens_1d;
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    hsize_t start_1d; 
    //hsize_t start_2d[2];  //Not Used
    hsize_t stride_1d;
    //hsize_t stride_2d[2];  //Not Used
    hsize_t count_1d; 
    //hsize_t count_2d[2];  //Not Used
    
    herr_t status;

    rank = 1;
    dimens_1d = numberOfBlocks;

    /* define the dataspace -- same as above */
    start_1d  = (hsize_t)idx;
    stride_1d = 1;
    count_1d  = 1;
    
    dataspace = H5Screate_simple(rank, &dimens_1d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			     &start_1d, &stride_1d, &count_1d, NULL);
    
    /* define the memory space */
    dimens_1d = 1;
    memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    
    /* create the dataset from scratch only if this is our first time */
    dataset = H5Dopen(datasetId, "refine level");
    status  = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
		      H5P_DEFAULT, (void *)&refinelevel);
  
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    return refinelevel;
}

void
FlashHDFFile::GetMinMax(char varname[5],
          		float *varMin,
          		float *varMax)
{
hid_t attr, dataset,status;

 dataset = H5Dopen(datasetId, varname);
 //get min and max for this variable

  attr = H5Aopen_name(dataset, "minimum");
  status = H5Aread(attr, H5T_NATIVE_FLOAT, varMin);
        if(status < 0) {
                fprintf(stderr, "Error reading minimum attribute!\n");
                }
  H5Aclose(attr);

  attr = H5Aopen_name(dataset, "maximum");
  status = H5Aread(attr, H5T_NATIVE_FLOAT, varMax);
        if(status < 0) {
                fprintf(stderr, "Error reading maximum attribute!\n");
                }
  H5Aclose(attr);

  H5Dclose(dataset);

}

float
FlashHDFFile::GetMinimumAndMaximumValue(char variableName[5], int start, int run, float r[2])
{

    r[0] = FLT_MAX;
    r[1] = -FLT_MAX;
    
    int nob = run;

    int dims[3];
    GetCellDimensions(dims);

    int totalsize = nob*dims[0]*dims[1]*dims[2];
    //int blocksize = dims[0]*dims[1]*dims[2];  //Not Used
    int *lt = new int[totalsize];
    
    float *variable = new float[totalsize];
    
    //    GetNodeType(start,run,lt);
    
    GetScalarVariable(variableName,start,run,variable);

    //int count = 0;  //Not Used
    float average = 0.0;

    for(int i = 0; i < totalsize; i++)
    {
      average += variable[i];
	if(variable[i] > r[1])
	{
	    r[1] = variable[i];
	}
	if(variable[i] < r[0])
	{
	    r[0] = variable[i];
	}
    }// End of i FOR
    
    cerr << "Average Value: " << average << endl;

    delete [] variable;
    delete [] lt;
    
    return average/totalsize;
}


#define nguard (0)// 4

int
FlashHDFFile::GetScalarVariable(char variableName[5], int dataPointIndex, float *variable)
{
    int rank;
    hsize_t dimens_4d[4];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    hsize_t start_4d[5];
    hsize_t stride_4d[5], count_4d[5];
    
    herr_t status;
    
    //char record_label_new[5];  //Not Used
    
    /* ------------------------== unknowns ==--------------------------- */
    if((cellDimensions[0] <= -1) || (cellDimensions[1] <= -1) || 
       (cellDimensions[2] <= -1))
    {
	_SetCellDimensions(); // Cell dimensions have not been set.
    }// End of IF
    	
    rank = 4;
    dimens_4d[0] = numberOfBlocks;    
    dimens_4d[1] = cellDimensions[2];
    dimens_4d[2] = cellDimensions[1];
    dimens_4d[3] = cellDimensions[0];
  /* define the dataspace -- as described above */
    start_4d[0]  = (hsize_t)dataPointIndex;
    start_4d[1]  = 0;
    start_4d[2]  = 0;
    start_4d[3]  = 0;
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0]  = 1;
    count_4d[1]  = cellDimensions[2];
    count_4d[2]  = cellDimensions[1];
    count_4d[3]  = cellDimensions[0];
    
    dataspace = H5Screate_simple(rank, dimens_4d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);

    int k3d,k2d;
    k3d = k2d = 1;// Hard coded bad bad bad
    
    /* 
       define the memory space -- the unkt array that was passed includes all
       the variables and guardcells for a single block.  We want to cut out the
       portion that just has the interior cells and the variable specified by 
       index
    */
    rank = 4;
    dimens_4d[0] = cellDimensions[2]+(nguard)*2*k3d;
    dimens_4d[1] = cellDimensions[1]+(nguard)*2*k2d;
    dimens_4d[2] = cellDimensions[0]+(nguard)*2;
    dimens_4d[3] = 1; // Number of Variables 1

    memspace = H5Screate_simple(rank, dimens_4d, NULL);

    /* exclude the guardcells and take only the desired variable */
    start_4d[0] = (nguard)*k3d;
    start_4d[1] = (nguard)*k2d;  
    start_4d[2] = (nguard);
    start_4d[3] = 0;//variableIndex;  /* should be 0 *//* remember: 0 based indexing */
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;  
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0] = cellDimensions[2];
    count_4d[1] = cellDimensions[1]; 
    count_4d[2] = cellDimensions[0]; 
    count_4d[3] = 1; 
    
    //    double *unknowns = new double[cellDimensions[2]*cellDimensions[1]*cellDimensions[0]]; // This is bad
    
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);
    
    dataset = H5Dopen(datasetId,variableName);
    status  = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		      H5P_DEFAULT, variable);
    
//      int i;

//      // This seems bad
//      for(i=0;i<cellDimensions[0]*cellDimensions[1]*cellDimensions[2];i++)
//  	variable[i] = (float)unknowns[i];


//      delete [] unknowns;
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 1;
    
}// End of

int
FlashHDFFile::GetScalarVariable(char variableName[5], int dataPointIndex, float bounds[6], float *variable)
{

    int rank;
    hsize_t dimens_4d[4];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    hsize_t start_4d[5];
    hsize_t stride_4d[5], count_4d[5];
    
    herr_t status;
    
    // char record_label_new[5];  //Not Used
    
    /* ------------------------== unknowns ==--------------------------- */
    if((cellDimensions[0] <= -1) || (cellDimensions[1] <= -1) || 
       (cellDimensions[2] <= -1))
    {
	_SetCellDimensions(); // Cell dimensions have not been set.
    }// End of IF
    	
    rank = 4;
    dimens_4d[0] = numberOfBlocks;    
    dimens_4d[1] = cellDimensions[2];
    dimens_4d[2] = cellDimensions[1];
    dimens_4d[3] = cellDimensions[0];

  /* define the dataspace -- as described above */
    start_4d[0]  = (hsize_t)dataPointIndex;
    start_4d[1]  = 0;
    start_4d[2]  = 0;
    start_4d[3]  = 0;
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0]  = 1;
    count_4d[1]  = cellDimensions[2];
    count_4d[2]  = cellDimensions[1];
    count_4d[3]  = cellDimensions[0];

    //    cout << "GSV: Setting up Source" << endl;
    
    dataspace = H5Screate_simple(rank, dimens_4d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,start_4d, stride_4d, count_4d, NULL);

    //    cout << "GSV: Finished Source" << endl;
    
    int k3d,k2d;
    k3d = k2d = 1;// Hard coded bad bad bad
    
    /* 
       define the memory space -- the unkt array that was passed includes all
       the variables and guardcells for a single block.  We want to cut out the
       portion that just has the interior cells and the variable specified by 
       index
    */
    rank = 4;
    dimens_4d[0] = cellDimensions[2]+(nguard)*2*k3d;
    dimens_4d[1] = cellDimensions[1]+(nguard)*2*k2d;
    dimens_4d[2] = cellDimensions[0]+(nguard)*2;
    dimens_4d[3] = 1; // Number of Variables 1

    memspace = H5Screate_simple(rank, dimens_4d, NULL);

    /* exclude the guardcells and take only the desired variable */
    start_4d[0] = (nguard)*k3d;
    start_4d[1] = (nguard)*k2d;  
    start_4d[2] = (nguard);
    start_4d[3] = 0;//variableIndex;  /* should be 0 *//* remember: 0 based indexing */
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;  
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0] = cellDimensions[2];
    count_4d[1] = cellDimensions[1]; 
    count_4d[2] = cellDimensions[0]; 
    count_4d[3] = 1; 
    
    double *unknowns = new double[cellDimensions[0]*cellDimensions[1]*cellDimensions[2]]; 

    
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);
    
    dataset = H5Dopen(datasetId,variableName);
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
		      H5P_DEFAULT, unknowns);
    
    int i;
    
    for(i=0;i<cellDimensions[0]*cellDimensions[1]*cellDimensions[2];i++)
	variable[i] = (float)unknowns[i];


    double minbounds[3],maxbounds[3];
    Get3dMinimumBounds(dataPointIndex,(double *)&minbounds);
    Get3dMaximumBounds(dataPointIndex,(double *)&maxbounds);

    bounds[0] = (float)minbounds[0];
    bounds[1] = (float)maxbounds[0];
    bounds[2] = (float)minbounds[1];
    bounds[3] = (float)maxbounds[1];    
    bounds[4] = (float)minbounds[2];
    bounds[5] = (float)maxbounds[2];
    
    delete [] unknowns;
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 1;
    
}// End of

void
FlashHDFFile::GetScalarVariable(char variableName[5], int dataPointIndex, int slices[6], float *variable)
{
                                                                                                                                                                                                                                
    int rank;
    hsize_t dimens_4d[4];
    hid_t ierr;
    int i,k,j,idx,count;
                                                                                                                                                                                                                                
    hid_t dataspace, memspace, dataset;
                                                                                                                                                                                                                                
    hsize_t start_4d[5];
    hsize_t stride_4d[5], count_4d[5];
                                                                                                                                                                                                                                
    herr_t status;

    idx = -1;
                                                                                                                                                                                                                                
    // char record_label_new[5];  //Not Used
                                                                                                                                                                                                                                
    /* ------------------------== unknowns ==--------------------------- */
    if((cellDimensions[0] <= -1) || (cellDimensions[1] <= -1) ||
       (cellDimensions[2] <= -1))
    {
        _SetCellDimensions(); // Cell dimensions have not been set.
    }// End of IF
                                                                                                                                                                                                                                
    rank = 4;
    dimens_4d[0] = numberOfBlocks;
    dimens_4d[1] = cellDimensions[2];
    dimens_4d[2] = cellDimensions[1];
    dimens_4d[3] = cellDimensions[0];
                                                                                                                                                                                                                                
  /* define the dataspace -- as described above */
    start_4d[0]  = (hsize_t)dataPointIndex;
    start_4d[1]  = 0;
    start_4d[2]  = 0;
    start_4d[3]  = 0;
                                                                                                                                                                                                                                
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
                                                                                                                                                                                                                                
    count_4d[0]  = 1;
    count_4d[1]  = cellDimensions[2];
    count_4d[2]  = cellDimensions[1];
    count_4d[3]  = cellDimensions[0];
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                
    dataspace = H5Screate_simple(rank, dimens_4d, NULL);
                                                                                                                                                                                                                                
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,start_4d, stride_4d, count_4d, NULL);
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                
    /*
       define the memory space -- the unkt array that was passed includes all
       the variables and guardcells for a single block.  We want to cut out the
       portion that just has the interior cells and the variable specified by
       index
    */
    rank = 4;
    dimens_4d[0] = cellDimensions[2];
    dimens_4d[1] = cellDimensions[1];
    dimens_4d[2] = cellDimensions[0];
    dimens_4d[3] = 1; // Number of Variables 1
                                                                                                                                                                                                                                
    memspace = H5Screate_simple(rank, dimens_4d, NULL);
                                                                                                                                                                                                                                
    /* exclude the guardcells and take only the desired variable */
    start_4d[0] = 0;
    start_4d[1] = 0;
    start_4d[2] = 0;
    start_4d[3] = 0;//variableIndex;  /* should be 0 *//* remember: 0 based indexing */
                                                                                                                                                                                                                                
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
                                                                                                                                                                                                                                
    count_4d[0] = cellDimensions[2];
    count_4d[1] = cellDimensions[1];
    count_4d[2] = cellDimensions[0];
    count_4d[3] = 1;
                                                                                                                                                                                                                                
    double *unknowns = new double[cellDimensions[0]*cellDimensions[1]*cellDimensions[2]];
                                                                                                                                                                                                                                
                                                                                                                                                                                                                                
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                               start_4d, stride_4d, count_4d, NULL);
                                                                                                                                                                                                                                
    dataset = H5Dopen(datasetId,variableName);
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                      H5P_DEFAULT, unknowns);
                                                                                                                                                                                                                                
    for(i=0;i<cellDimensions[0]*cellDimensions[1]*cellDimensions[2];i++)
        unknowns[i] = (float)unknowns[i];

   count = 0;
   for(i=slices[4]; i<(slices[4] + slices[5]); i++) {
	for(k = slices[2]; k<(slices[2] + slices[3]); k++) {
		for(j = slices[0]; j<(slices[0] + slices[1]); j++) {
			idx = j + (k * cellDimensions[1]) + (i * (cellDimensions[2]*cellDimensions[2]));
			variable[count] = unknowns[idx];
			count++;
			}
		}
	}

delete [] unknowns;
}



int
FlashHDFFile::GetScalarVariable(char variableName[5], int dataPointIndex, int sizeofrun, float *variable)
{
    //    cout << "In GSV" << endl;
    
    int rank;
    hsize_t dimens_4d[4];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    hsize_t start_4d[5];
    hsize_t stride_4d[5], count_4d[5];
    
    herr_t status;
    
    //char record_label_new[5];  //Not Used
    int k3d,k2d;
    k3d = k2d = 1;// Hard coded bad bad bad
    
    /* ------------------------== unknowns ==--------------------------- */
    if((cellDimensions[0] <= -1) || (cellDimensions[1] <= -1) || 
       (cellDimensions[2] <= -1))
    {
	_SetCellDimensions(); // Cell dimensions have not been set.
    }// End of IF
    	
    rank = 4;
    dimens_4d[0] = numberOfBlocks;    
    dimens_4d[1] = cellDimensions[2]; 
    dimens_4d[2] = cellDimensions[1];
    dimens_4d[3] = cellDimensions[0];
    
  /* define the dataspace -- as described above */
    start_4d[0]  = (hsize_t)dataPointIndex;
    start_4d[1]  = 0;
    start_4d[2]  = 0;
    start_4d[3]  = 0;
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0]  = sizeofrun;// was 1, now the number of blocks grabbed
    count_4d[1]  = cellDimensions[2];
    count_4d[2]  = cellDimensions[1];
    count_4d[3]  = cellDimensions[0];

    //    cout << "GSV: Source" << endl;
    
    dataspace = H5Screate_simple(rank, dimens_4d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);

    /* 
       define the memory space -- the unkt array that was passed includes all
       the variables and guardcells for a single block.  We want to cut out the
       portion that just has the interior cells and the variable specified by 
       index
    */
    rank = 4;
    dimens_4d[0] = numberOfBlocks;
    dimens_4d[1] = cellDimensions[2]+(nguard)*2*k3d;
    dimens_4d[2] = cellDimensions[1]+(nguard)*2*k2d;
    dimens_4d[3] = cellDimensions[0]+(nguard)*2;
    //    dimens_4d[3] = sizeofrun; // Number of Variables 1

    memspace = H5Screate_simple(rank, dimens_4d, NULL);

    /* exclude the guardcells and take only the desired variable */
    start_4d[0] = (nguard)*k3d;
    start_4d[1] = (nguard)*k2d;  
    start_4d[2] = (nguard);
    start_4d[3] = 0;//variableIndex;  /* should be 0 *//* remember: 0 based indexing */
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;  
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[1] = cellDimensions[2];
    count_4d[2] = cellDimensions[1]; 
    count_4d[3] = cellDimensions[0]; 
    count_4d[0] = sizeofrun; 
    
    // This is bad why hard coded ...
    //    double *unknowns = new double[sizeofrun*cellDimensions[2]*cellDimensions[1]*cellDimensions[0]];

    //    cout << "GSV: Destination" << endl;
    
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);

    //    cout << "GSV: Error " << ierr << endl;
    
    dataset = H5Dopen(datasetId,variableName);
    if(dataset < 0)
    {
	cout << "Failed to Open Dataset, Using Variable " << variableName << endl;
    }// End of IF
    
    status  = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
		      H5P_DEFAULT, variable);

    //    cout << "GSV: Read" << endl;
    
//      int i;
    
//      // This is bad why recopy 
//      for(i=0;i<sizeofrun*cellDimensions[0]*cellDimensions[1]*cellDimensions[2];i++)
//      {
//  	//	cout << unknowns[i] << endl;
//  	variable[i] = (float)unknowns[i];
//      }

//      delete [] unknowns;
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 1;
    
}// End of

int
FlashHDFFile::GetScalarVariable(char variableName[5], int dataPointIndex, int sizeofrun,float bounds[6], float *variable)
{
//    cout << "Loading Variable " << variableName << endl;
    int rank;
    hsize_t dimens_4d[4];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    hsize_t start_4d[5];
    hsize_t stride_4d[5], count_4d[5];
    
    herr_t status;
    
    //char record_label_new[5];  //Not Used
    
    /* ------------------------== unknowns ==--------------------------- */
    if((cellDimensions[0] <= -1) || (cellDimensions[1] <= -1) || 
       (cellDimensions[2] <= -1))
    {
	_SetCellDimensions(); // Cell dimensions have not been set.
    }// End of IF
    	
    rank = 4;
    dimens_4d[0] = numberOfBlocks;    
    dimens_4d[1] = cellDimensions[2]; 
    dimens_4d[2] = cellDimensions[1];
    dimens_4d[3] = cellDimensions[0];

  /* define the dataspace -- as described above */
    start_4d[0]  = (hsize_t)dataPointIndex;
    start_4d[1]  = 0;
    start_4d[2]  = 0;
    start_4d[3]  = 0;
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0]  = sizeofrun;// was 1, now the number of blocks grabbed
    count_4d[1]  = cellDimensions[2];
    count_4d[2]  = cellDimensions[1];
    count_4d[3]  = cellDimensions[0];
    
    dataspace = H5Screate_simple(rank, dimens_4d, NULL);

    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);

    int k3d,k2d;
    k3d = k2d = 1;// Hard coded bad bad bad
    
    /* 
       define the memory space -- the unkt array that was passed includes all
       the variables and guardcells for a single block.  We want to cut out the
       portion that just has the interior cells and the variable specified by 
       index
    */
    rank = 5;
    dimens_4d[0] = numberOfBlocks;
    dimens_4d[1] = cellDimensions[2]+(nguard)*2*k3d;
    dimens_4d[2] = cellDimensions[1]+(nguard)*2*k2d;
    dimens_4d[3] = cellDimensions[0]+(nguard)*2;

    memspace = H5Screate_simple(rank, dimens_4d, NULL);

    /* exclude the guardcells and take only the desired variable */
    start_4d[0] = 0;
    start_4d[1] = (nguard)*k3d;
    start_4d[2] = (nguard)*k2d;  
    start_4d[3] = (nguard);
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;  
    stride_4d[2] = 1;
    stride_4d[3] = 1;
    
    count_4d[0] = sizeofrun;
    
    count_4d[1] = cellDimensions[2];
    count_4d[2] = cellDimensions[1]; 
    count_4d[3] = cellDimensions[0]; 
    count_4d[4] = 1; 
    
    // This is bad why hard coded ...
    double *unknowns = new double[sizeofrun*cellDimensions[2]*cellDimensions[1]*cellDimensions[0]];
    
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
			       start_4d, stride_4d, count_4d, NULL);
    
    dataset = H5Dopen(datasetId,variableName);
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
		      H5P_DEFAULT, unknowns);
    

    int i;
    
    // This is bad why recopy 
    for(i=0;i<sizeofrun*cellDimensions[0]*cellDimensions[1]*cellDimensions[2];i++)
	variable[i] = (float)unknowns[i];


    delete [] unknowns;
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    double minbounds[3],maxbounds[3];
    //    Get3dMinimumBounds(dataPointIndex, sizeofrun, (double *)&minbounds);
    //    Get3dMaximumBounds(dataPointIndex, sizeofrun, (double *)&maxbounds);

    bounds[0] = (float)minbounds[0];
    bounds[1] = (float)maxbounds[0];
    bounds[2] = (float)minbounds[1];
    bounds[3] = (float)maxbounds[1];    
    bounds[4] = (float)minbounds[2];
    bounds[5] = (float)maxbounds[2];    

    return 1;
    
}// End of

int
FlashHDFFile::Get3dMinimumBounds(int idx, double coords[3])
{
    float coord[3]; 
    float bound[3];
    
    Get3dCoordinate(idx,coord);
    Get3dBlockSize(idx,bound);

    coords[0] = coord[0] - (bound[0]/2.0);
    coords[1] = coord[1] - (bound[1]/2.0);
    coords[2] = coord[2] - (bound[2]/2.0);

    return 1;  //Return some default value
}// End of

int
FlashHDFFile::Get3dBoundingBox(int idx, float bnd_boxes[6])
{

  float bnd_box[6];
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_3d[3];

  hsize_t start_3d[3];
  hsize_t stride_3d[3], count_3d[3];

 dataset = H5Dopen(datasetId, "bounding box");
 dataspace = H5Dget_space(dataset);

  start_3d[0] = (hssize_t) (idx);
  start_3d[1] = 0;
  start_3d[2] = 0;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  count_3d[0] = (hssize_t) (1);
  count_3d[1] = 3;
  count_3d[2] = 2;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_3d,
                              stride_3d, count_3d, NULL);



  rank = 3;
  dimens_3d[0] = numberOfBlocks;
  dimens_3d[1] = 3;
  dimens_3d[2] = 2;

  memspace = H5Screate_simple(rank, dimens_3d, NULL);

  start_3d[0] = 0;
  start_3d[1] = 0;
  start_3d[2] = 0;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  count_3d[0] = 1;
  count_3d[1] = 3;
  count_3d[2] = 2;

   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                             start_3d, stride_3d, count_3d, NULL);

   status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                   H5P_DEFAULT, bnd_box);

    bnd_boxes[0] = (float)bnd_box[0];
    bnd_boxes[1] = (float)bnd_box[1];
    bnd_boxes[2] = (float)bnd_box[2];
    bnd_boxes[3] = (float)bnd_box[3];
    bnd_boxes[4] = (float)bnd_box[4];
    bnd_boxes[5] = (float)bnd_box[5];

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return 1;
}

int
FlashHDFFile::Get3dMaximumBounds(int idx, double coords[3])
{
    float coord[3];
    float bound[3];
    
    Get3dCoordinate(idx,coord);
    Get3dBlockSize(idx,bound);

    coords[0] = coord[0] + (bound[0]/2.0);
    coords[1] = coord[1] + (bound[1]/2.0);
    coords[2] = coord[2] + (bound[2]/2.0);
    
    return 1;  //Return Some default value
}// End of


	
int
FlashHDFFile::GetCoordinateRangeEntireDataset(float ranges[6])
{
    int i;
    float coord[3];
    float bounds[3];
    float emin[3];
    float emax[3];
    
    ranges[0] = ranges[2] = ranges[4] = FLT_MAX;
    ranges[1] = ranges[3] = ranges[5] = -FLT_MAX;
    
    for(i=0;i<numberOfBlocks;i++)
    {
	Get3dCoordinate(i,coord);
	Get3dBlockSize(i,bounds);
	
	emin[0] = coord[0] - bounds[0]/2.0;
  	emax[0] = coord[0] + bounds[0]/2.0;
	
  	emin[1] = coord[1] - bounds[1]/2.0;
  	emax[1] = coord[1] + bounds[1]/2.0;
	
  	emin[2] = coord[2] - bounds[2]/2.0;
  	emax[2] = coord[2] + bounds[2]/2.0;
	
  	if(emin[0] < ranges[0])
  	    ranges[0] = (float)emin[0];
  	if(emin[1] < ranges[2])
  	    ranges[2] = (float)emin[1];
  	if(emin[2] < ranges[4])
  	    ranges[4] = (float)emin[2];
	
  	if(emax[0] > ranges[1])
  	    ranges[1] = (float)emax[0];
  	if(emax[1] > ranges[3])
  	    ranges[3] = (float)emax[1];
  	if(emax[2] > ranges[5])
  	    ranges[5] = (float)emax[2];
    }// End of FOR

    return 1;
      
}// End of GetCoordinateRangeEntireDataset(float range[6])

int
FlashHDFFile::Get3dCoordinate(int idx, float coords[3])
{
    double coord[3];
    int rank;
    hsize_t dimens_1d, dimens_2d[2];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    // hsize_t start_1d;  //Not Used
    hsize_t start_2d[2];
    //hsize_t stride_1d; //Not Used
    hsize_t stride_2d[2];
    //hsize_t count_1d; //Not Used
    hsize_t count_2d[2];
    herr_t status;
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;
    
    /* define the dataspace -- as described above */
    start_2d[0]  = (hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    count_2d[0]  = 1;
    count_2d[1]  = 3;//NDIM;
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
  			       start_2d, stride_2d, count_2d, NULL);
    
    /* define the memory space */
    rank = 1;
    dimens_1d = 3;//NDIM;
    memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    
    dataset = H5Dopen(datasetId, "coordinates");
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
  		      H5P_DEFAULT, coord);
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    coords[0] = (float)coord[0];
    coords[1] = (float)coord[1];
    coords[2] = (float)coord[2];
    
    return 1;
    
}// End of Get3dCoordinate(int, float[3])

    
int
FlashHDFFile::Get3dCoordinate(int idx, int sizeofrun, float *coords)
{
    double *coord = new double[3*sizeofrun];
    int rank;
    //hsize_t dimens_1d; //Not Used
    hsize_t dimens_2d[2];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;
    
    //hsize_t start_1d; //Not Used
    hsize_t start_2d[2];
    // hsize_t stride_1d; //Not Used
    hsize_t stride_2d[2];
    //hsize_t count_1d;  //Not Used
    hsize_t count_2d[2];
    
    herr_t status;
    
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;
    
    /* define the dataspace -- as described above */
    start_2d[0]  = (hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = sizeofrun;
    count_2d[1]  = 3;//NDIM;
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
  			       start_2d, stride_2d, count_2d, NULL);
    
    /* define the memory space */
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;

    memspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    /* define the dataspace -- as described above */
    start_2d[0]  = 0;//(hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = sizeofrun;
    count_2d[1]  = 3;//NDIM;

    
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
  			       start_2d, stride_2d, count_2d, NULL);
    
    dataset = H5Dopen(datasetId, "coordinates");
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
  		      H5P_DEFAULT, coord);

    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    for(int i=0;i<(sizeofrun*3);i++)
      coords[i] = (float)coord[i];

    delete [] coord;

    return 1;
    
}// End of Get3dCoordinate(int, float[3])

int
FlashHDFFile::Get3dBlockSize(int idx, float coords[3])
{
    double size[3];
    
    int rank;
    hsize_t dimens_1d, dimens_2d[2];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;

    //hsize_t start_1d;  //Not Used
    hsize_t start_2d[2];
    hsize_t stride_2d[2], count_2d[2];
    
    herr_t status;
    
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;
    
    /* define the dataspace -- as described above */
    start_2d[0]  = (hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = 1;
    count_2d[1]  = 3;//NDIM;
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_2d, stride_2d, count_2d, NULL);
    
    /* define the memory space */
    rank = 1;
    dimens_1d = 3;//NDIM;
    memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    
    dataset = H5Dopen(datasetId, "block size");
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
  		      H5P_DEFAULT, size);
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    coords[0] = (float)size[0];
    coords[1] = (float)size[1];
    coords[2] = (float)size[2];
    
    return 1;  // Return some default value
}// End of Get3dBlockSize(int idx, float coords[3])

int
FlashHDFFile::Get3dBlockSize(int idx, int sizeofrun, float *coords)
{
    double *size = new double[3*sizeofrun];
    
    int rank;
    //hsize_t dimens_1d;  //Not Used
    hsize_t dimens_2d[2];
    hid_t ierr;
    
    hid_t dataspace, memspace, dataset;

    //hsize_t start_1d;  //Not Used
    hsize_t start_2d[2];
    hsize_t stride_2d[2], count_2d[2];
    
    herr_t status;
    
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;
    
    /* define the dataspace -- as described above */
    start_2d[0]  = (hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = sizeofrun;
    count_2d[1]  = 3;//NDIM;
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			       start_2d, stride_2d, count_2d, NULL);
    
    /* define the memory space */
    rank = 2;
    dimens_2d[0] = numberOfBlocks;
    dimens_2d[1] = 3;//NDIM;

    memspace = H5Screate_simple(rank, dimens_2d, NULL);
    
    /* define the dataspace -- as described above */
    start_2d[0]  = 0;//(hsize_t)idx;
    start_2d[1]  = 0;
    
    stride_2d[0] = 1;
    stride_2d[1] = 1;
    
    count_2d[0]  = sizeofrun;
    count_2d[1]  = 3;//NDIM;    rank = 1;

    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
  			       start_2d, stride_2d, count_2d, NULL);
    
    dataset = H5Dopen(datasetId, "block size");
    status  = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
  		      H5P_DEFAULT, size);
    
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    for(int i=0;i<(sizeofrun*3);i++)
      coords[i] = (float)size[i];

    delete [] size;
    
    return 1;  //Return Some Default Value
}// End of Get3dBlockSize(int idx, int sizeofrun, float *coords)

int
FlashHDFFile::GetNumberOfVariables()
{
    hid_t dataspace, dataset;

    dataset = H5Dopen(datasetId, "unknown names");
    dataspace = H5Dget_space(dataset);
    
    hsize_t fields[2];
    hsize_t maxfields[2];
    
    H5Sget_simple_extent_dims (dataspace, fields, maxfields);

    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    return fields[0];
    
}// End GetNumberOfVariables()

int
FlashHDFFile::GetVariableNames(FlashVariableNames *name)
{
    int numberOfVariables = GetNumberOfVariables();
    
    char *tname = new char[4*numberOfVariables]; 

    //int rank = 2;  //Not Used
    int dimens_2d[2];
    dimens_2d[0] = numberOfVariables;
    dimens_2d[1] = 1;
    
    /* manually set the string size */
    //int string_size = 4;  //Not Used
    
    hid_t dataset = H5Dopen(datasetId, "unknown names");
    hid_t dataspace = H5Dget_space(dataset);

    /* setup the datatype for this string length */
    //    hid_t string_type = H5Tcopy(H5T_C_S1);
    //    H5Tset_size(string_type, string_size);
    hid_t string_type = H5Dget_type(dataset);


    /*hid_t status    = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, 
      H5P_DEFAULT, tname);*/

    H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, 
	    H5P_DEFAULT, tname);

    //name = new FlashVariableNames[numberOfVariables];

    for(int i=0;i<numberOfVariables;i++)
    {	
	strncpy(name[i],&tname[i*4],4);
	name[i][4] = '\0';
    }
    
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return numberOfVariables;
}

// Particle Data
// Returns negative number on failure, if successful returns the number of particles.
int
FlashHDFFile::GetParticleData(FlashParticleData **particleData)
{
  //int ierr;  //Not Used
  //int FAIL = -1;  //Not Used
  //hsize_t dimens_1d; //Not Used

	hid_t sp_type;
	
	hid_t dataset = H5Dopen(datasetId, "particle tracers");
	if(dataset < 0)
	{
	  debugprint((char *)"Particle Dataset Failed to Open");
	    return FAILED;
	}// End of IF
	
	hid_t dataspace = H5Dget_space(dataset);
	
	hsize_t fields[2];
	hsize_t maxfields[2];
	
	H5Sget_simple_extent_dims (dataspace, fields, maxfields);

	*particleData = new FlashParticleData[fields[0]];
	
	/* create the HDF 5 compound data type to describe the record */
	sp_type = H5Tcreate(H5T_COMPOUND, sizeof(FlashParticleData));
	
	H5Tinsert(sp_type,"tag",offsetof(FlashParticleData,id),H5T_NATIVE_INT);
	H5Tinsert(sp_type,"block_no",offsetof(FlashParticleData,blockLocation),H5T_NATIVE_INT);
	H5Tinsert(sp_type,"x pos",offsetof(FlashParticleData,x),H5T_NATIVE_DOUBLE);
	H5Tinsert(sp_type,"y pos",offsetof(FlashParticleData,y),H5T_NATIVE_DOUBLE);
	H5Tinsert(sp_type,"z pos",offsetof(FlashParticleData,z),H5T_NATIVE_DOUBLE);
	H5Tinsert(sp_type,"x vel",offsetof(FlashParticleData,vx),H5T_NATIVE_DOUBLE);
	H5Tinsert(sp_type,"y vel",offsetof(FlashParticleData,vy),H5T_NATIVE_DOUBLE);
	H5Tinsert(sp_type,"z vel",offsetof(FlashParticleData,vz),H5T_NATIVE_DOUBLE);

	status = H5Dread(dataset, sp_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, *particleData);

	H5Tclose(sp_type);
	H5Dclose(dataset);

	return fields[0]; // Number Of Particles
	
}// End of GetParticleData

//-----------------------------------------------------------------------
//
// GetFloatVecBlocks
//
// reads a (2D) vector of floats for selected blocks in the file
// vector can have size d = 1 (scalar)
//
// var: variable (dataset) name
//
// sblock, nblocks: starting block and number of blocks
// if leaf_only = 1, sblock and nblocks are in terms of leaf blocks
// else, sblock and nblocks are in terms of all file blocks
//
// d: size of vector in each block
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// Tom Peterka, 11/6/09
//
int FlashHDFFile::GetFloatVecBlocks(char *var, int sblock, int nblocks,
					int d, int leaf_only, float *data) {

  int rank = 2;
  hsize_t dims[2];
  hid_t dataspace, memspace, dataset;
  hsize_t start[2];
  int skip = 0; // number of nonleaf blocks skipped
  int i, j;
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d;
  start[1] = 0;

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id,
		 data) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // scale data
  if (scale != 1.0f) {
    for (i = 0; i < (int)dims[0]; i++) {
      if (!leaf_only || GetNodeType(start[0] + i) == 1) {
	for (j = 0; j < d; j++)
	  data[i * d + j] *= scale;
      }
    }
  }
 
  // remove nonleaf blocks (in-place data copy)
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++) {
	if (skip)
	  data[(i - skip) * d + j] = data[i * d + j];
      }
    }
    return(dims[0] - skip);

  }

  else
    return dims[0];
    
}
//-----------------------------------------------------------------------
//
// GetDoubleVecBlocks
//
// reads a (2D) vector of doubles for selected blocks in the file
// vector can have size d = 1 (scalar)
// output is converted to floats
//
// var: variable (dataset) name
//
// sblock, nblocks: starting block and number of blocks
// if leaf_only = 1, sblock and nblocks are in terms of leaf blocks
// else, sblock and nblocks are in terms of all file blocks
//
// d: size of vector in each block
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// Tom Peterka, 11/6/09
//
int FlashHDFFile::GetDoubleVecBlocks(char *var, int sblock, int nblocks,
					int d, int leaf_only, float *data) {

  int rank = 2;
  hsize_t dims[2];
  hid_t dataspace, memspace, dataset;
  hsize_t start[2];
  int skip = 0; // number of nonleaf blocks skipped
  int i, j;
  double *dd; // double version of data
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d;
  start[1] = 0;

  assert((dd = new double[dims[0] * dims[1]]) != NULL);

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id,
		 dd) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // partial copy from double to float, removing nonleaf blocks
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++)
	  data[(i - skip) * d + j] = dd[i * d + j] * scale;
    }
    delete[] dd;
    return(dims[0] - skip);

  }

  // full copy from double to float
  else {
    for (i = 0; i < (int)dims[0]; i++) {
      for (j = 0; j < d; j++)
	data[i * d + j] = dd[i * d + j] * scale;
    }
    delete[] dd;
    return dims[0];
  }
    
}
//-----------------------------------------------------------------------
//
// GetFloatMatBlocks
//
// reads a (2D) matrix of  floats for selected blocks in the file
// eg. 3x2 block bounds
//
// var: variable (dataset) name
// sblock, nblocks: starting block and number of blocks
// d1, d2: matrix size in the order listed in the file
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d1 * d2 * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// added by Tom Peterka, 11/2/09
//
int FlashHDFFile::GetFloatMatBlocks(char *var, int sblock, int nblocks, 
				    int d1, int d2, int leaf_only, 
				    float *data) {

  int rank = 3;
  hsize_t dims[3];
  hid_t dataspace, memspace, dataset;
  hsize_t start[3];
  int skip = 0; // number of nonleaf blocks skipped
  int d = d1 * d2; // number of elements per block
  int i, j;
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d1;
  dims[2] = d2;
  start[1] = start[2] = 0;

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id,
		 data) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // scale data
  if (scale != 1.0f) {
    for (i = 0; i < (int)dims[0]; i++) {
      if (!leaf_only || GetNodeType(start[0] + i) == 1) {
	for (j = 0; j < d; j++)
	  data[i * d + j] *= scale;
      }
    }
  }
 
  // remove nonleaf blocks (in-place data copy)
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++) {
	if (skip)
	  data[(i - skip) * d + j] = data[i * d + j];
      }
    }
    return(dims[0] - skip);

  }

  else
    return dims[0];
    
}
//-----------------------------------------------------------------------
//
// GetDoubleMatBlocks
//
// reads a (2D) matrix of doubles for selected blocks in the file
// eg. 3x2 block bounds
// output is converted to floats
//
// var: variable (dataset) name
// sblock, nblocks: starting block and number of blocks
// d1, d2: matrix size in the order listed in the file
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d1 * d2 * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// added by Tom Peterka, 11/2/09
//
int FlashHDFFile::GetDoubleMatBlocks(char *var, int sblock, int nblocks, 
				     int d1, int d2, int leaf_only, 
				     float *data) {

  int rank = 3;
  hsize_t dims[3];
  hid_t dataspace, memspace, dataset;
  hsize_t start[3];
  int skip = 0; // number of nonleaf blocks skipped
  int d = d1 * d2; // number of elements per block
  int i, j;
  double *dd; // double version of data
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d1;
  dims[2] = d2;
  start[1] = start[2] = 0;

  assert((dd = new double[dims[0] * dims[1] * dims[2]]) != NULL);

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id,
		 dd) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // partial copy from double to float, removing nonleaf blocks
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++)
	  data[(i - skip) * d + j] = dd[i * d + j] * scale;
    }
    delete[] dd;
    return(dims[0] - skip);

  }

  // full copy from double to float
  else {
    for (i = 0; i < (int)dims[0]; i++) {
      for (j = 0; j < d; j++)
	data[i * d + j] = dd[i * d + j] * scale;
    }
    delete[] dd;
    return dims[0];
  }
    
}
//-----------------------------------------------------------------------
//
// GetFloatVolBlocks
//
// reads a (3D) volume of floats for selected blocks in the file
//
// var: variable (dataset) name
//
// sblock, nblocks: starting block and number of blocks
// if leaf_only = 1, sblock and nblocks are in terms of leaf blocks
// else, sblock and nblocks are in terms of all file blocks
//
// d1, d2, d3: size of volume in each block
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// Tom Peterka, 11/6/09
//
int FlashHDFFile::GetFloatVolBlocks(char *var, int sblock, int nblocks,
					int d1, int d2, int d3, int leaf_only,
					float *data) {

  int rank = 4;
  hsize_t dims[4];
  hid_t dataspace, memspace, dataset;
  hsize_t start[4];
  int skip = 0; // number of nonleaf blocks skipped
  int d = d1 * d2 * d3; // total size of the block
  int i, j;
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d1;
  dims[2] = d2;
  dims[3] = d3;
  start[1] = start[2] = start[3] = 0;

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id,
		 data) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // scale data
  if (scale != 1.0f) {
    for (i = 0; i < (int)dims[0]; i++) {
      if (!leaf_only || GetNodeType(start[0] + i) == 1) {
	for (j = 0; j < d; j++)
	  data[i * d + j] *= scale;
      }
    }
  }
 
  // remove nonleaf blocks (in-place data copy)
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++) {
	if (skip)
	  data[(i - skip) * d + j] = data[i * d + j];
      }
    }
    return(dims[0] - skip);

  }

  else
    return dims[0];
    
}
//-----------------------------------------------------------------------
//
// GetDoubleVolBlocks
//
// reads a (3D) volume of doubles for selected blocks in the file
// output is converted to floats
//
// var: variable (dataset) name
//
// sblock, nblocks: starting block and number of blocks
// if leaf_only = 1, sblock and nblocks are in terms of leaf blocks
// else, sblock and nblocks are in terms of all file blocks
//
// d1, d2, d3: size of volume in each block
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
//
int FlashHDFFile::GetDoubleVolBlocks(char *var, int sblock, int nblocks,
				     int d1, int d2, int d3, int leaf_only,
				     float *data) {

  int rank = 4;
  hsize_t dims[4];
  hid_t dataspace, memspace, dataset;
  hsize_t start[4];
  int skip = 0; // number of nonleaf blocks skipped
  int d = d1 * d2 * d3; // total size of the block
  int i, j;
  double *dd; // double version of data
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d1;
  dims[2] = d2;
  dims[3] = d3;
  start[1] = start[2] = start[3] = 0;

  assert((dd = new double[dims[0] * dims[1] * dims[2] * dims[3]]) != NULL);

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id,
		 dd) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // partial copy from double to float, removing nonleaf blocks
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++)
	  data[(i - skip) * d + j] = dd[i * d + j] * scale;
    }

    delete[] dd;
    return(dims[0] - skip);

  }

  // full copy from double to float
  else {
    for (i = 0; i < (int)dims[0]; i++) {
      for (j = 0; j < d; j++)
	data[i * d + j] = dd[i * d + j] * scale;
    }
    delete[] dd;
    return dims[0];
  }
    
}
//-----------------------------------------------------------------------
//
// GetIntVecBlocks
//
// reads a (1D) vector of ints for all blocks in the file
// eg. block centers or block sizes
//
// var: variable (dataset) name
// sblock, nblocks: starting block and number of blocks
// d: size of vector in each block
// leaf_only: 0 = keep all nodes, 1 = keep only leaf nodes
// data: enough space to hold all the result (d * total number blocks)
//
// returns: number of blocks stored (leaf + nonleaf or only leaf)
//
// added by Tom Peterka, 11/2/09
//
int FlashHDFFile::GetIntVecBlocks(char *var, int sblock, int nblocks, 
				     int d, int leaf_only, int *data) {

  int rank = 2;
  hsize_t dims[2];
  hid_t dataspace, memspace, dataset;
  hsize_t start[2];
  int skip = 0; // number of nonleaf blocks skipped
  int i, j;
    
  // bracket the file blocks to read (start[0], dims[0])
  if (!leaf_only) {
    dims[0] = nblocks;
    start[0] = sblock;
  }
  else {
    j = -1; // count of leaf blocks
    for (i = 0; j < sblock && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    start[0] = i - 1;
    j = 0;
    for ( ; j < nblocks && i < numberOfBlocks; i++) {
      if (GetNodeType(i) != 1)
	continue;
      j++;
    }
    dims[0] = i - start[0];
  }

  dims[1] = d;
  start[1] = 0;

  // read the data    
  memspace = H5Screate_simple(rank, dims, NULL);    
  assert((dataset = H5Dopen(datasetId, var)) >= 0);
  dataspace = H5Dget_space(dataset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
		      start, NULL, dims, NULL);
  assert(H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, plist_id,
		 data) >= 0);

  // cleanup
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  // remove nonleaf blocks (in-place data copy)
  if (leaf_only) {

    for (i = 0; i < (int)dims[0]; i++) {
      if (GetNodeType(start[0] + i) != 1) {
	skip++;
	continue;
      }
      for (j = 0; j < d; j++) {
	if (skip)
	  data[(i - skip) * d + j] = data[i * d + j];
      }
    }
    return(dims[0] - skip);

  }

  else
    return dims[0];
    
}
//-----------------------------------------------------------------------

#if 0
int
FlashHDFFile::GetVectorVariable(int variableIndex, int dataPointIndex, float *variable)
{
    int index = SDnametoindex(fileId, "unknowns");
    int id = SDselect(fileId, index);
    
    char sds_name[64];
    int32 rank, dimsizes[MAX_VAR_DIMS], dataType, num_attrs;
	    
    // Retrieves the name, rank, dimension sizes, data type and number of attributes for a data set.
    istat = SDgetinfo(id, sds_name, &rank, dimsizes, &dataType, &num_attrs);
    HDF_ERROR_CHECK(istat,"flashGetCoordinateRangeEntireDataset:1");
    
    numberOfBlocks = dimsizes[0];
    int zRes = dimsizes[1];
    int yRes = dimsizes[2];
    int xRes = dimsizes[3];
    int numberOfVariables = dimsizes[4];
    
#if 0
    cout << "Rank: " << rank << endl;
    cout << "Number of blocks: " << numberOfBlocks << endl;
    cout << "X " << xRes << " Y " << yRes << " Z " << zRes << endl;
    cout << "Number of Variables: " << numberOfVariables << endl;
#endif
    int numberOfDimensions = 0;
    if(xRes > 1)
	numberOfDimensions++;
    if(yRes > 1)
	numberOfDimensions++;
    if(zRes > 1)
	numberOfDimensions++;
    
//    cout << "Number of Dimensions: " << numberOfDimensions << endl;
    
    int k2d = numberOfDimensions/2;
    int k3d = (numberOfDimensions-1)/2;
    
    void *data;
    int cellSizeX = (xRes+nguard-(1+nguard)+1);
    int cellSizeY = (yRes+nguard*k2d-(1+nguard*k2d)+1);
    int cellSizeZ = (zRes+nguard*k3d-(1+nguard*k3d)+1);

//    cout << "Cell X " << cellSizeX << " ";
//    cout << "Cell Y " << cellSizeY << " ";
//    cout << "Cell Z " << cellSizeZ << endl;

    if (dataType == DFNT_FLOAT32) 
	data = calloc(numberOfVariables * cellSizeX * cellSizeY * cellSizeZ, sizeof(float));
    else if (dataType == DFNT_FLOAT64) 
	data = calloc(numberOfVariables * cellSizeX * cellSizeY * cellSizeZ, sizeof(double));
    else 
    {
	cout << "ERROR: Unknown datatype" << endl;
	exit(1);
    }// End of ELSE

    int32 start_uk[5], stride_uk[5], edges_uk[5];
    
    start_uk[1] = start_uk[2] = start_uk[3] = start_uk[4] = 0;
    start_uk[0] = dataPointIndex;
    stride_uk[0] = stride_uk[1] = stride_uk[2] = stride_uk[3] = stride_uk[4] = 1;
	
    //	unknowns
    edges_uk[4] = numberOfVariables;
    
    edges_uk[3] = xRes;
    edges_uk[2] = yRes;
    edges_uk[1] = zRes;
    edges_uk[0] = 1;//numberOfBlocks;

    
    istat = SDreaddata(id, (int32 *)start_uk, (int32 *)stride_uk, (int32 *)edges_uk, data);
    HDF_ERROR_CHECK(istat,"Variable Reader");

    int i,j,k,n;
    int internalVariableCount = 0;
    int vectorComponent = 0;
    
    for(k=0; k<cellSizeZ; k++) 
    {
	for(j=0; j<cellSizeY; j++) 
	    for(i=0; i<cellSizeX; i++) 
		for(n=0; n<numberOfVariables; n++) 
		{
		    if((n >= variableIndex) && (n < (variableIndex+numberOfDimensions))) // Variable of Choice
		    {	
			vectorComponent = abs(n-variableIndex);
			
			if(dataType == DFNT_FLOAT32)
			    variable[(k*cellSizeX*cellSizeY*numberOfDimensions)+(j*cellSizeX*numberOfDimensions)+(i*numberOfDimensions)+vectorComponent]
				= ((float*)data)[internalVariableCount];
			else if(dataType == DFNT_FLOAT64)
			    variable[(k*cellSizeX*cellSizeY*numberOfDimensions)+(j*cellSizeX*numberOfDimensions)+(i*numberOfDimensions)+vectorComponent]
				= (float)((double*)data)[internalVariableCount];
			else 
			{
			    cout << "datatype_unknown is invalid" << endl;
			    exit(1);
			}// End of ELSE
		    }// End of IF
		    internalVariableCount++;
		}// End of n FOR
    }// End of k FOR
    
    return 1;
    
}// End of


/////////////////
int32 
sdPrintInformation(int32 id)
{
    char sds_name[64];
    int32 rank, dimsizes[MAX_VAR_DIMS], data_type, num_attrs;
    int32 istat;
    
    // Retrieves the name, rank, dimension sizes, data type and number of attributes for a data set.
    istat = SDgetinfo(id, sds_name, &rank, dimsizes, &data_type, &num_attrs);
    HDF_ERROR_CHECK(istat,"sdPrintInformation");

#if PRINT_HDF_INFO
    cout << endl << "sds_name: " << sds_name << endl;
    cout << "rank: " << rank << endl;
    for(int i = 0;i<rank;i++)
        cout << "dimsizes[" << i << "]: " << dimsizes[i] << endl;
    cout << "data_type: " << data_type << endl;
    cout << "num_attrs: " << num_attrs << endl;
#endif

    return data_type;
}// End of sdPrintInformation() 


#endif

#endif // HDF5
