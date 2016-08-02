// convert .vec files to vti format
#include <stdio.h>
#include <stdlib.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include "vtkZLibDataCompressor.h"


int main ( int argc, char *argv[] )
{
	printf("Usage: vec2vti file.raw w h d\n");
	printf("Output: output.vti\n");
	FILE *fp = fopen(argv[1], "rb");
	if (!fp) {
		perror(argv[1]);
		exit(1);
	}

	int i, j, k;
	float x[3], v[3];
	int dims[3];

    dims[0] = atoi(argv[2]);
    dims[1] = atoi(argv[3]);
    dims[2] = atoi(argv[4]);

	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->SetDimensions(dims);

#if 0
	// We also create the points and vectors.
	vtkFloatArray *vectors = vtkFloatArray::New();
	vectors->SetNumberOfComponents(3);
	vectors->SetNumberOfTuples(dims[0]*dims[1]*dims[2]); // allocation done here
	fread(vectors->GetPointer(0), 12, dims[0]*dims[1]*dims[2], fp);

	vectors->SetName("Velocity");
	imageData->GetPointData()->AddArray(vectors);
	vectors->Delete();
#else
	#if VTK_MAJOR_VERSION <= 5
	imageData->SetNumberOfScalarComponents(1);
	imageData->SetScalarTypeToFloat();
	#else
	imageData->AllocateScalars(VTK_FLOAT, 1);
	#endif

	fread(imageData->GetScalarPointer(0,0,0), 4, dims[0]*dims[1]*dims[2], fp);
#endif

    // compressor
    vtkSmartPointer<vtkDataCompressor> compressor = vtkZLibDataCompressor::New();


	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("output.vti");
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(imageData->GetProducerPort());
	#else
	writer->SetInputData(imageData);
	#endif
    writer->SetCompressor(compressor);
    writer->EncodeAppendedDataOff ();

	writer->Write();

	fclose(fp);

    printf("Done.  Output: output.vti\n");
    return EXIT_SUCCESS;
}
