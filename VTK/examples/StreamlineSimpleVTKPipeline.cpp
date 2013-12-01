// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py

// Compatibility with Paraview
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL) 
// Compatibility with Paraview

#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

#include "vtkOSUFlow.h"
#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkSmartPointer.h"
#include "vtkLineSource.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkStructuredGridOutlineFilter.h"
#include "vtkProperty.h"
#include "vtkXMLStructuredGridReader.h"
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkImageData.h>
// streamline
#include "vtkStreamLine.h"

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
	//vtkSmartPointer<vtkDataSet> data = getData_plot3d();
	//vtkSmartPointer<vtkDataSet> data = getData_vts();
	vtkSmartPointer<vtkDataSet> data = getData_vti();

	//
	// Determine seeds
	//
	// use a static rake to generate streamlines
	vtkLineSource *rake = vtkLineSource::New();
	rake->SetPoint1(15, -5, 32);
	rake->SetPoint2(15, 5, 32);
	rake->SetResolution(21);
	vtkPolyDataMapper *rakeMapper = vtkPolyDataMapper::New();
	rakeMapper->SetInputConnection(rake->GetOutputPort());
	vtkActor *rakeActor = vtkActor::New();
	rakeActor->SetMapper(rakeMapper);
	rakeActor->VisibilityOn();

#if 1
	// use vtkOSUFlow
	vtkOSUFlow *streamer = vtkOSUFlow::New();
	streamer->SetInputData(data);
	streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationStepLength(.001);
	streamer->SetIntegrationDirectionToBackward();
	streamer->SetMaximumPropagationTime(200);
	streamer->SetNumberOfThreads(1);
	streamer->SetIntegrationDirectionToBackward();//  IntegrateBothDirections();
	//streamer->VorticityOn();
#else
	// use vtkStreamLine
	vtkStreamLine *streamer = vtkStreamLine::New();
	streamer->SetInputData(data);
	streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetMaximumPropagationTime(200);
	streamer->SetIntegrationStepLength(.001);
	streamer->SetStepLength(.001);
	streamer->SetNumberOfThreads(1);
	streamer->SetIntegrationDirectionToBackward();//  IntegrateBothDirections();
	streamer->VorticityOff();
#endif
	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(streamer->GetOutputPort());
	mapper->SetScalarRange(data->GetScalarRange());
	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);


	//
	// outline
	//
	vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
	outline->SetInputData(data);

	vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
	outlineMapper->SetInputConnection(outline->GetOutputPort());

	vtkActor *outlineActor = vtkActor::New();
	outlineActor->SetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(0,0,0);

	//
	// renderer
	//
	vtkRenderer *ren = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	ren->AddActor(rakeActor);
	ren->AddActor(actor);
	ren->AddActor(outlineActor);
	ren->SetBackground(.5,.5,.5);

	renWin->SetSize(500,500);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}



