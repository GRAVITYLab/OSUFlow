	#ifndef OSUFLOW_VTK
#define OSUFLOW_VTK

#include <assert.h>

#include <vtkSmartPointer.h>

#include "OSUFlow.h"
#include "VectorFieldVTK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

class vtkDataSet;
class vtkImageData;

class OSUFlowVTK: public OSUFlow
{
protected:
	vtkSmartPointer<vtkImageData> sImageData;
public:
	virtual ~OSUFlowVTK() { }



	// openmp support
	inline void initOpenMP(int nproc) {
#ifdef _OPENMP
		omp_set_num_threads(nproc);
		printf("Openmp max threads = %d\n", omp_get_max_threads());
#endif
	}

	void setData(vtkDataSet* input);

	inline bool getHasData() { return this->has_data; }

#if 0
	virtual void LoadData(const char* fname, bool bStatic, bool deferred = false)
	{	OSUFlow::LoadData(fname, bStatic, deferred);}

	virtual void LoadData(const char* fname, bool bStatic, VECTOR3 pMin, VECTOR3 pMax, bool deferred = false)
	{ assert(false); OSUFlow::LoadData(fname, bStatic, pMin, pMax, deferred); }

	virtual void LoadData(const char* fname, bool bStatic, VECTOR3 pMin,
		VECTOR3 pMax, int min_t, int max_t, bool deferred = false)
	{ assert(false); OSUFlow::LoadData(fname, bStatic, pMin, pMax, min_t, max_t, deferred); }

	virtual void LoadData(const char* fname, float *sMin, float *sMax,
		float *dim, int min_t, int max_t, DataMode mode,
		float **data = NULL)
	{ assert(false); OSUFlow::LoadData(fname, sMin, sMax, dim, min_t, max_t, mode, data); }

	virtual void LoadData(char **dataset_files, int num_dataset_files,
		float *sMin, float *sMax, float *dim, int min_t, int max_t,
		DataMode mode, float **data = NULL)
	{ assert(false); OSUFlow::LoadData(dataset_files, num_dataset_files, sMin, sMax, dim, min_t, max_t, mode, data); }

	//////////////////////////////////////////////////////////////////////////
	// sReal : does not include ghost cells
	//////////////////////////////////////////////////////////////////////////
	virtual void LoadData(char **dataset_files, int num_dataset_files,
				       float *sMin, float *sMax, int* sRealMin, int* sRealMax,
				       float *dim, int min_t, int max_t, DataMode mode,
				       float **data = NULL);
#endif

};

#endif //OSUFLOW_VTK
