#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

#include "OSUFlow.h"
#include "VectorFieldVTK.h"
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include <vtkInterpolatedVelocityField.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkFunctionSet.h>
#include <vtkImageInterpolator.h>
#include <vtkPointData.h>

using namespace std;

vtkSmartPointer<vtkDataSet> getData_plot3d()
{
	// Start by loading some data.
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	pl3dReader->SetXYZFileName(	SAMPLE_DATA_DIR "/curvilinear/combxyz.bin");
	pl3dReader->SetQFileName(SAMPLE_DATA_DIR "/curvilinear/combq.bin");
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();

	// random points
	//vtkStructuredGrid *grid = pl3dReader->GetOutput();
	//int *dim = grid->GetDimensions();
	printf("Number of blocks=%d\n", pl3dReader->GetOutput()->GetNumberOfBlocks());
	vtkSmartPointer<vtkDataSet> data = vtkDataSet::SafeDownCast( pl3dReader->GetOutput()->GetBlock(0) );
	pl3dReader->Delete();


	return data;
}

// structured VTK format
vtkSmartPointer<vtkDataSet> getData_vts()
{
	vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
	reader->SetFileName(SAMPLE_DATA_DIR "/curvilinear/comb.vts"); //TODO

#if 0
	int extent[6];
	extent[0] = sMin[0]; extent[1] = sMax[0];
	extent[2] = sMin[1]; extent[3] = sMax[1];
	extent[4] = sMin[2]; extent[5] = sMax[2];
	printf("Set Extent: %d %d %d %d %d %d\n", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);
	reader->SetUpdateExtent(0, extent);
#endif
	//reader->UpdateInformation();
	reader->Update();

	//reader->PrintSelf(std::cout, vtkIndent(2));

	int *ext = reader->GetOutput()->GetExtent();
	printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

	vtkSmartPointer<vtkDataSet> data = vtkDataSet::SafeDownCast( reader->GetOutput() );
	//data->GetPointData()->SetActiveAttribute("Momentum", vtkDataSetAttributes::VECTORS);

	printf("File read\n");
	// show content
	data->PrintSelf(std::cout, vtkIndent(2));

	reader->Delete();


	// test
	vtkInterpolatedVelocityField *interpolator = vtkInterpolatedVelocityField::New();
	interpolator->AddDataSet(data);
	double vel[3];
	double coords[4] = {2,2,2,0};
	for (int k=0; k<3; k++)
	{
		bool success = ((vtkFunctionSet *)interpolator)->FunctionValues(coords, vel);
		printf("querry success = %d, vec = %lf %lf %lf\n", success, vel[0], vel[1], vel[2]);
	}

	return data;
}


// VTK image format
vtkSmartPointer<vtkDataSet> getData_vti()
{
	vtkXMLImageDataReader *reader = vtkXMLImageDataReader::New();
	reader->SetFileName(SAMPLE_DATA_DIR "/regular/tornado/1.vti"); //TODO
	reader->UpdateInformation();
	reader->Update();

	int *ext = reader->GetOutput()->GetExtent();
	printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);

	vtkSmartPointer<vtkImageData > data = reader->GetOutput() ;

	printf("File read\n");
	// show content
	data->PrintSelf(std::cout, vtkIndent(2));
	// debug
	printf("pointer: %p\n", data->GetScalarPointer());


	reader->Delete();



	return data;
}



int main()
{
	// choose one of the following:
	vtkSmartPointer<vtkDataSet> data = getData_plot3d();
	//vtkSmartPointer<vtkDataSet> data = getData_vts();
	//vtkSmartPointer<vtkDataSet> data = getData_vti(); // use CVectorField


	OSUFlow *osuflow = new OSUFlow;
	CVectorField *field = new VectorFieldVTK( data );
	osuflow->SetFlowField( field );

	// gen seeds
	VECTOR3 minB, maxB;
	osuflow->Boundary(minB, maxB);
	printf("bounds: %lf %lf %lf %lf %lf %lf\n", minB[0], maxB[0], minB[1], maxB[1], minB[2], maxB[2]);

	osuflow->SetRandomSeedPoints(&minB[0], &maxB[0], 100);


	// openmp
#ifdef _OPENMP
	osuflow->initOpenMP(8);
#endif

	// init results
	list<vtListSeedTrace*> list;
	osuflow->SetIntegrationParams(1, 5);
	osuflow->GenStreamLines(list , BACKWARD_AND_FORWARD, 50, 0);
	printf(" done integrations\n");
	printf("list size = %d\n", (int)list.size());

	std::list<vtListSeedTrace*>::iterator pIter;


	pIter = list.begin();
	for (; pIter!=list.end(); pIter++) {
		vtListSeedTrace *trace = *pIter;
		std::list<VECTOR3*>::iterator pnIter;
		pnIter = trace->begin();
		for (; pnIter!=trace->end(); pnIter++) {
		  VECTOR3 p = **pnIter;
		  printf(" %f %f %f ", p[0], p[1], p[2]);
		}
		printf("\n");
	}

	delete osuflow;
}



