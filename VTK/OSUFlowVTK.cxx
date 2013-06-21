#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>
#if 0
#include <vtkXMLFileReadTester.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkGenericDataObjectReader.h>
#endif

#include "OSUFlowVTK.h"

void OSUFlowVTK::setData(vtkDataSet *input)
{
	vtkImageData *imageData = vtkImageData::SafeDownCast(input);
	if (NULL == imageData) {  // not Regular-grid data

		this->flowField = new VectorFieldVTK(input);

	} else
	{
		// this is a RAW data format.  Use OSUFlow's CFlowField interpolator
		int numPoints = imageData->GetNumberOfPoints();

		// copy data
		float **ppData = new float*[1]; // time varying
		ppData[0] = new float[numPoints*3];
		assert(ppData[0]);
		memcpy(ppData[0], imageData->GetScalarPointer(), numPoints*12);

		// get info
		float sMin[3], sMax[3], dim[3];
		int sRealMin[3], sRealMax[3];
		int *bound = imageData->GetExtent();

		sMin[0] = bound[0]; sMax[0] = bound[1];
		sMin[1] = bound[2]; sMax[1] = bound[3];
		sMin[2] = bound[4]; sMax[2] = bound[5];
		sRealMin[0] = sMin[0]; sRealMax[0] = sMax[0];
		sRealMin[1] = sMin[1]; sRealMax[1] = sMax[1];
		sRealMin[2] = sMin[2]; sRealMax[2] = sMax[2];
		dim[0] = sRealMax[0]-sRealMin[0]+1;
		dim[1] = sRealMax[1]-sRealMin[1]+1;
		dim[2] = sRealMax[2]-sRealMin[2]+1;
		int t_min =0;
		int t_max = 0;

		this->InitFlowField(sMin, sMax, sRealMin, sRealMax, dim, t_min, t_max, RAW, ppData);

	}

	// entire data, static
	this->flowField->Boundary(gMin, gMax);
	printf("boundary: %f %f %f %f %f %f\n", gMin[0], gMax[0], gMin[1], gMax[1], gMin[2], gMax[2]);
	lMin = gMin;
	lMax = gMax;
	MinT = 0;
	MaxT = 0;

	has_data = true;
}

void OSUFlowVTK::LoadData(char **dataset_files, int num_dataset_files,
			       float *sMin, float *sMax, int* sRealMin, int* sRealMax,
			       float *dim, int min_t, int max_t, DataMode mode,
			       float **data )
{
	if (has_data)
		this->DeleteData();

	// currently: load in VTI data format
	printf("sMin: %f %f %f\n", sMin[0], sMin[1], sMin[2]);
	printf("sMax: %f %f %f\n", sMax[0], sMax[1], sMax[2]);
	printf("sRealMin: %d %d %d\n", sRealMin[0], sRealMin[1], sRealMin[2]);
	printf("sRealMax: %d %d %d\n", sRealMax[0], sRealMax[1], sRealMax[2]);

	this->dataset_files = dataset_files;
	this->num_dataset_files = num_dataset_files;
	bStaticFlow = false;


	lMin.Set(sMin[0], sMin[1], sMin[2]);
	lMax.Set(sMax[0], sMax[1], sMax[2]);

	if (max_t >= min_t) {
		numTimesteps = max_t - min_t + 1;
		MinT = min_t; MaxT = max_t;
	}
	else { // defaults to 1 time step
		numTimesteps = 1;
		MinT = MaxT = min_t;
	}

	if (numTimesteps > num_dataset_files) {
		fprintf(stderr, "[OSUFlowVTK] Error: time step range larger than number of files\n");
		exit(1);
	}

	float **ppData = new float *[numTimesteps];
	for (int t=MinT; t<=MaxT; t++)
	{
		vtkXMLImageDataReader *reader = vtkXMLImageDataReader::New();
		reader->SetFileName(dataset_files[t]); //TODO: now assuming timesteps always start from 0
		reader->UpdateInformation();

		// set local extent
		int extent[6];
		extent[0] = sMin[0]; extent[1] = sMax[0];
		extent[2] = sMin[1]; extent[3] = sMax[1];
		extent[4] = sMin[2]; extent[5] = sMax[2];
		printf("Set Extent: %d %d %d %d %d %d\n", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
		reader->SetUpdateExtent(0, extent);
		reader->Update();

		int *ext = reader->GetOutput()->GetExtent();
		printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

		vtkImageData *image = reader->GetOutput() ;

		// show content
		//data->PrintSelf(std::cout, vtkIndent(2));
		int numPoints = image->GetNumberOfPoints();

		// copy data
		float *pData = new float[numPoints*3];
		assert(pData);
		memcpy(pData, image->GetScalarPointer(), numPoints*12);

		ppData[t-MinT] = pData;

		//image->Delete();
		//reader->Delete();
	}

	InitFlowField(sMin, sMax, sRealMin, sRealMax, dim, min_t, max_t, mode, ppData);

	printf("File read\n");

#if 0 // unstructured
	//InitFlowField(sMin, sMax, sRealMin, sRealMax, dim, min_t, max_t, mode,data);

	vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
	reader->SetFileName(dataset_files[0]); //TODO
	reader->UpdateInformation();

	// set local extent
	int extent[6];
	extent[0] = sMin[0]; extent[1] = sMax[0];
	extent[2] = sMin[1]; extent[3] = sMax[1];
	extent[4] = sMin[2]; extent[5] = sMax[2];
	printf("Set Extent: %d %d %d %d %d %d\n", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
	reader->SetUpdateExtent(0, extent);
	reader->Update();

	int *ext = reader->GetOutput()->GetExtent();
	printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

	// set data
	vtkStructuredGrid *data = reader->GetOutput();
	this->flowField = new VectorFieldVTK(data);

	reader->Delete();

	VECTOR3 minB, maxB;
	this->flowField->Boundary(minB, maxB);

	break;
#endif

	has_data = true;

}
#if 0
vtkMultiBlockDataSet *OSUFlowVTK::loadVTKData(const char *fileName)
{
	vtkXMLFileReadTester *vtkXMLFormatFileChecker = vtkXMLFileReadTester::New();

	vtkXMLFormatFileChecker->SetFileName(fileName);

	if (vtkXMLFormatFileChecker->TestReadFile() > 0) {
		//Logger.getLogger(getClass().getName()).log(Level.INFO, "[{0}] is of type [{1}]",
		//		new Object[] { fileName, vtkXMLFormatFileChecker.GetFileDataType() });

		vtkXMLGenericDataObjectReader *vtkXMLFileReader = vtkXMLGenericDataObjectReader::New();

		vtkXMLFileReader->SetFileName(fileName);
		vtkXMLFileReader->Update();

		if (vtkXMLFileReader->GetOutput()->IsA("vtkMultiBlockDataSet")) {
			return vtkXMLFileReader->GetOutput();
		} else {
			vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

			makeMultiblock->SetInputConnection(vtkXMLFileReader->GetOutputPort());

			return makeMultiblock;
		}
	} else // legacy format
	{
		vtkGenericDataObjectReader *legacyVTKFileReader = vtkGenericDataObjectReader::New();

		legacyVTKFileReader->SetFileName(fileName);
		legacyVTKFileReader->Update();

		vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

		makeMultiblock->SetInput(legacyVTKFileReader->GetOutput());

		return makeMultiblock;
	}
}
#endif
