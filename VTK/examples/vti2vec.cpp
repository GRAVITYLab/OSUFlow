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
#include <vtkXMLImageDataReader.h>

typedef struct {double x,y,z;} double3;

// VTK image format
vtkSmartPointer<vtkImageData> getData_vti(const char *fname)
{
	vtkXMLImageDataReader *reader = vtkXMLImageDataReader::New();
	reader->SetFileName(fname); //TODO
	reader->UpdateInformation();
	reader->Update();

	int *ext = reader->GetOutput()->GetExtent();
	printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

	vtkSmartPointer<vtkImageData > data = reader->GetOutput() ;

	reader->Delete();

	return data;
}

int main ( int argc, char *argv[] )
{
	printf("Usage: vti2vec file.vti\n");
	printf("Output: output.vec\n\n");

	vtkSmartPointer<vtkImageData> imageData = getData_vti(argv[1]);

	int dims[3];
	imageData->GetDimensions(dims);
	printf("dim: %d %d %d\n", dims[0], dims[1], dims[2]);
	printf("Number of scalar components: %d\n", imageData->GetNumberOfScalarComponents());
	if (imageData->GetNumberOfScalarComponents() != 3) {
		printf("Program does not support scalar components != 3.\n");
		exit(1);
	}

	// output file
	FILE *fp = fopen("output.vec", "wb");
	if (!fp) {
		perror("output.vec");
		exit(1);
	}

	int i, j, k;
	int size = dims[0]*dims[1]*dims[2];
	if (size != imageData->GetNumberOfPoints()) {
		printf("Data size not matching dimension!\n");
	}

	fwrite(dims, 3, 4, fp);
	if (imageData->GetScalarType() == VTK_FLOAT) {
		fwrite(imageData->GetScalarPointer(), size, 4*3	, fp);
	} else {
		double3 *p = (double3 *)imageData->GetScalarPointer();
		for (i=0; i<size; i++) {
			float floatp[3];
			floatp[0] = p[i].x;
			floatp[1] = p[i].y;
			floatp[2] = p[i].z;
			fwrite(floatp, 3, 4, fp);
		}
	}

	fclose(fp);
	printf("Done\n");

	return 0;
}
