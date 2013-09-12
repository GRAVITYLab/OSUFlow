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
// streamline
#include "vtkStreamLine.h"

using namespace std;


int main()
{
	// Start by loading some data.
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	// set data
	{
		char path[256];
		sprintf(path, "%s/curvilinear/combxyz.bin", SAMPLE_DATA_DIR);
		printf("%s\n", path);
		pl3dReader->SetXYZFileName(path);
		sprintf(path, "%s/curvilinear/combq.bin", SAMPLE_DATA_DIR);
		pl3dReader->SetQFileName(path);
	}
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();
	vtkDataSet *data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));
	data->Register(NULL);
	pl3dReader->Delete();

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

#if 0
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



