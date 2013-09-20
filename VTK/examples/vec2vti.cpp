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

int main ( int argc, char *argv[] )
{
	printf("Usage: vec2vti file.vec\n");
	printf("Output: output.vti\n");
	FILE *fp = fopen(argv[1], "rb");
	if (!fp) {
		perror(argv[1]);
		exit(1);
	}

	int i, j, k;
	float x[3], v[3];
	int dims[3];

	fread(dims, 3, 4, fp);
	printf("dim: %d %d %d\n", dims[0], dims[1], dims[2]);


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
	imageData->SetNumberOfScalarComponents(3);
	imageData->SetScalarTypeToFloat();
	#else
	imageData->AllocateScalars(VTK_FLOAT, 3);
	#endif

	fread(imageData->GetScalarPointer(0,0,0), 12, dims[0]*dims[1]*dims[2], fp);
#endif

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("output.vti");
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(imageData->GetProducerPort());
	#else
	writer->SetInputData(imageData);
	#endif
	writer->Write();

	fclose(fp);

  printf("Done.  Output: output.vti\n");
  return EXIT_SUCCESS;
}
