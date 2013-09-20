// convert .vec files to VTS format
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
	printf("Usage: vec2vts file.vec\n");
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

	// Create the structured grid.
	vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
	sgrid->SetDimensions(dims);

	// We also create the points and vectors.
	vtkFloatArray *vectors = vtkFloatArray::New();
	vectors->SetNumberOfComponents(3);
	vectors->SetNumberOfTuples(dims[0]*dims[1]*dims[2]);
	vtkPoints *points = vtkPoints::New();
	points->Allocate(dims[0]*dims[1]*dims[2]);

	vtkIdType id=0;
	for ( k=0; k<dims[2]; k++)
	{
	for (j=0; j<dims[1]; j++)
	  {
	  for (i=0; i<dims[0]; i++)
		{
		  x[0] = i;
		  x[1] = j;
		  x[2] = k;
		  fread(v, 1, 3, fp);
		  points->InsertPoint(id, x);
		  vectors->InsertTuple(id,v);
		  id++;
		}
	  }
	}
	//sgrid->SetPoints(points);
	points->Delete();

	vectors->SetName("velocity");
	//sgrid->GetPointData()->SetVectors(vectors);
	vectors->Delete();

	  // Write file
	  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
	  writer->SetFileName("output.vts");
	#if VTK_MAJOR_VERSION <= 5
	  writer->SetInput(structuredGrid);
	#else
	  writer->SetInputData(sgrid);
	#endif
	  writer->Write();
	fclose(fp);
 
  printf("Done\n");
  return EXIT_SUCCESS;
}
