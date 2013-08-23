/*  Wrappers for FLASH HDF Files 
 */
#ifndef __FLASH_HDF__
#define __FLASH_HDF__

#define COORD_SIZE 3
#define BND_SIZE 6

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _MPI
#include <mpi.h>
#endif

#define MAX_STRING_LENGTH 80

#ifdef HDF5
#include <hdf5.h>
#endif

using namespace std;

const int FAILED = -1;

inline void debugprint(char *message)
{
#if DEBUG
  cout << message << endl;
#endif
}

void HDF_ERROR_CHECK(int errorflag, char *message);

typedef char FlashVariableNames[5];

typedef struct FLASHPARTICLEDATA
{
  int id;
  int blockLocation;
    
  double x,y,z;
  double vx,vy,vz;
}FlashParticleData;


#ifdef HDF5

class FlashHDFFile
{
 public:
  // Default Constructor
  FlashHDFFile();
    
  // Constructor that opens a given HDF file
  FlashHDFFile(char *filename, float scale); 

#ifdef _MPI
  // constructor that opens a given HDF file for collective I/O
  FlashHDFFile(char *filename, float scale, MPI_Comm comm);
#endif

  // Open given HDF File
  int Open(char *filename);
    
  // Close currently opened HDF File
  int Close();

  // Get the number of dimensions in the current dataset {i.e. 1, 2, 3}
  int GetNumberOfDimensions();
    
  // Get the number of blocks in the current dataset
  int GetNumberOfBlocks();
  int GetNumBlocks();
    
  // Get the number of leaf blocks in the current dataset
  int GetNumberOfLeafs();

  // Get the number of cells per block !! This is bad, needs another function and this should be a pointer.
  void GetCellDimensions(int dimensions[3]);
    
  // Get the node type for cell at index
  int GetNodeType(int index);
  int GetNodeType(int index,int runsize,int *nodearray);

  // Get corner data for cell at index
  int Get3dMinimumBounds(int index, double mbbox[3]);
  int Get3dMaximumBounds(int index, double mbbox[3]);
    
  // Spatial extent of dataset (Boundaries)
  int GetCoordinateRangeEntireDataset(float ranges[6]);

  // Spatial location of a given block (index)
  int Get3dCoordinate(int index, float coords[3]);
  int Get3dCoordinate(int index, int sizeofrun, float *coords);

  // Size of a given block (index)
  int Get3dBlockSize(int index, float coords[3]);
  int Get3dBlockSize(int index, int sizeofrun, float *coords);

  //bnd box
  int Get3dBoundingBox(int idx, float bounds[6]);

  // Gets the minimum and maximum value starting at cell "start" over "run" blocks returned in "r", float return value is the average !! this seems bad
  float GetMinimumAndMaximumValue(char variableName[5], int start, int run, float r[2]);

  // !! This is not even implemented, what the ??
  int GetMinimumAndMaximumValue(char variableName[5], int start, int run, float r[2], int b[2]);
    
  // Spatial location of a given block (index)
  // int Get2dCoordinate(int index, float coords[2]);

  // Size of a given block (index)
  // int Get2dBlockSize(int index, float coords[2]);

  // How many variables are available
  int GetNumberOfVariables();
  // Name of variables
  int GetVariableNames(FlashVariableNames *names);
    
  // Vector data is returned in an interleaved format (x,y,z,x,y,z) and 
  // the the dimensions are automatically generated based on the 
  // dimensionallity of the dataset
  // int GetVectorVariable(int variableIndex, int index, float *variable);

  // Get scalar data -> variable name, index to block, pointer variable, runlength is for grabbing multiple blocks, bounds is spatial
  int GetScalarVariable(char variableName[5], int index, float *variable);
  void GetScalarVariable(char variableName[5], int index, int slices[6], float *variable);
  int GetScalarVariable(char variableName[5], int index, int runlength, float *variable);
  int GetScalarVariable(char variableName[5], int index, float bounds[6], float *variable);
  int GetScalarVariable(char variableName[5], int index, int runlength, float bounds[6], float *variable);
    
  int GetRefinementLevel(int index);

  int GetGlobalIds(int index,int globalIds[]);
  int GetNumberOfGlobalIds();

  // Particle Stuff
  int GetParticleData(FlashParticleData **particleData);

  //min max of scalar
  void GetMinMax(char variableName[5], float *min, float * max);

  // get (1D) vector of floats for selected blocks
  // data in the file can be float (GetFloat...) or double (GetDouble...)
  // output is float either way
  int GetFloatVecBlocks(char *var, int sblock, int nblocks,
			    int d, int leaf_only, float *data);
  int GetDoubleVecBlocks(char *var, int sblock, int nblocks,
			 int d, int leaf_only, float *data);

  // get (2D) matrix of floats for selected blocks
  // data in the file can be float (GetFloat...) or double (GetDouble...)
  // output is float either way
  int GetFloatMatBlocks(char *var, int sblock, int nblocks, 
			int d1, int d2, int leaf_only, float *data);
  int GetDoubleMatBlocks(char *var, int sblock, int nblocks, 
			int d1, int d2, int leaf_only, float *data);

  // get (3D) volume of floats for selected blocks
  // data in the file can be float (GetFloat...) or double (GetDouble...)
  // output is float either way
  int GetFloatVolBlocks(char *var, int sblock, int nblocks,
			    int d1, int d2, int d3, int leaf_only,
			    float *data);
  int GetDoubleVolBlocks(char *var, int sblock, int nblocks,
			 int d1, int d2, int d3, int leaf_only,
			 float *data);

  // get (2D) vector of ints for selected blocks
  int GetIntVecBlocks(char *var, int sblock, int nblocks, int d,
		      int leaf_only, int *data);



 private:

  void _SetCellDimensions();
  void _InitSettings();
  void _ResetSettings(int constructor);

  hid_t datasetId;  // Pointer to opened HDF5 File
  herr_t status;    // Commonly used return value

  int numberOfDimensions; // Number of dimensions in the dataset
  int numberOfBlocks;
  int numberOfLeafs;
    
  float *min,*max; // Minimum and maximum directions
  int *node_type; // node type (leaf = 1)

  int cellDimensions[3];
    
  hid_t plist_id; // HDF property list for collective or independent I/O
  float scale; // scaling factor for data

};

#endif // HDF5

#ifndef HDF5

// class with empty versions of all functions
// makes no sense to call this without HDF defined, but needed to compile
// OSUFLOW if there is not HDF installed on the system

class FlashHDFFile
{
 public:
  FlashHDFFile(){}
  FlashHDFFile(char *filename, float scale){}
#ifdef _MPI
  FlashHDFFile(char *filename, float scale, MPI_Comm comm){}
#endif
  int Open(char *filename){ return 0; }
  int Close(){ return 0; }
  int GetNumberOfDimensions(){ return 0; }
  int GetNumberOfBlocks(){ return 0; }
  int GetNumBlocks(){ return 0; }
  int GetNumberOfLeafs(){ return 0; }
  void GetCellDimensions(int dimensions[3]){}
  int GetNodeType(int index){ return 0; }
  int GetNodeType(int index,int runsize,int *nodearray){ return 0; }
  int Get3dMinimumBounds(int index, double mbbox[3]){ return 0; }
  int Get3dMaximumBounds(int index, double mbbox[3]){ return 0; }
  int GetCoordinateRangeEntireDataset(float ranges[6]){ return 0; }
  int Get3dCoordinate(int index, float coords[3]){ return 0; }
  int Get3dCoordinate(int index, int sizeofrun, float *coords){ return 0; }
  int Get3dBlockSize(int index, float coords[3]){ return 0; }
  int Get3dBlockSize(int index, int sizeofrun, float *coords){ return 0; }
  int Get3dBoundingBox(int idx, float bounds[6]){ return 0; }
  float GetMinimumAndMaximumValue(char variableName[5], int start, int run, float r[2]){ return 0.0; }
  int GetMinimumAndMaximumValue(char variableName[5], int start, int run, float r[2], int b[2]){ return 0; }
  int GetNumberOfVariables(){ return 0; }
  int GetVariableNames(FlashVariableNames *names){ return 0; }
  int GetScalarVariable(char variableName[5], int index, float *variable){ return 0; }
  void GetScalarVariable(char variableName[5], int index, int slices[6], float *variable){}
  int GetScalarVariable(char variableName[5], int index, int runlength, float *variable){ return 0; }
  int GetScalarVariable(char variableName[5], int index, float bounds[6], float *variable){ return 0; }
  int GetScalarVariable(char variableName[5], int index, int runlength, float bounds[6], float *variable){ return 0; }
  int GetRefinementLevel(int index){ return 0; }
  int GetGlobalIds(int index,int globalIds[]){ return 0; }
  int GetNumberOfGlobalIds(){ return 0; }
  int GetParticleData(FlashParticleData **particleData){ return 0; }
  void GetMinMax(char variableName[5], float *min, float * max){}
  int GetFloatVecBlocks(char *var, int sblock, int nblocks,
			int d, int leaf_only, float *data){ return 0; }
  int GetFloatMatBlocks(char *var, int sblock, int nblocks, 
			int d1, int d2, int leaf_only, float *data){ return 0; }
  int GetFloatVolBlocks(char *var, int sblock, int nblocks,
			int d1, int d2, int d3, int leaf_only,
			float *data){ return 0; }
  int GetDoubleVecBlocks(char *var, int sblock, int nblocks,
			int d, int leaf_only, float *data){ return 0; }
  int GetDoubleMatBlocks(char *var, int sblock, int nblocks, 
			int d1, int d2, int leaf_only, float *data){ return 0; }
  int GetDoubleVolBlocks(char *var, int sblock, int nblocks,
			int d1, int d2, int d3, int leaf_only,
			float *data){ return 0; }
  int GetIntVecBlocks(char *var, int sblock, int nblocks, int d,
		      int leaf_only, int *data){ return 0; }

 private:

  void _SetCellDimensions(){}
  void _InitSettings(){}
  void _ResetSettings(int constructor){}
  int numberOfDimensions;
  int numberOfBlocks;
  int numberOfLeafs;
  float *min,*max;
  int *node_type;
  int cellDimensions[3];
  float scale;

};

#endif // not HDF5

#endif
