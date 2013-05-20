
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
#include "vtkTesting.h"
#include "vtkProperty.h"
#include "vtkLineWidget.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
// streamline
#include "vtkStreamLine.h"

using namespace std;

vtkOSUFlow *streamer;
void beginInteraction(vtkObject *caller, long unsigned int eventId, void *callData )
{
}

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("compute\n");
	streamer->Update();
}


int main()
{
	// Start by loading some data.
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	// set data
	{
		char path[256];
		vtkTesting *t = vtkTesting::New();
		sprintf(path, "%s/Data/combxyz.bin", "/home/jchen/project/VTKData"); //t->GetDataRoot());
		printf("%s\n", path);
		pl3dReader->SetXYZFileName(path);
		sprintf(path, "%s/Data/combq.bin", "/home/jchen/project/VTKData"); //t->GetDataRoot());
		pl3dReader->SetQFileName(path);
		t->Delete();
	}
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();
	vtkDataSet *data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));

	//
	// Determine seeds
	//
#if 0
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
#else
	// user can change the rake
	vtkLineWidget *lineWidget = vtkLineWidget::New();
	lineWidget->SetInputData(data);
	lineWidget->SetAlignToYAxis();
	lineWidget->PlaceWidget();
	lineWidget->ClampToBoundsOn();
	vtkPolyData *rake = vtkPolyData::New();
	lineWidget->GetPolyData(rake);
#endif

#if 1
	// vtkOSUFlow
	streamer = vtkOSUFlow::New();
	streamer->SetInputData(data);
	streamer->SetSourceData(rake);
//	/streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetStepLength(.001);
	streamer->SetIntegrationDirectionToBackward();
	streamer->SetMaximumPropagationTime(200);
	streamer->SetNumberOfThreads(1);
	streamer->SetIntegrationDirectionToBackward();//  IntegrateBothDirections();
	streamer->VorticityOn();
#else
	vtkStreamLine *streamer = vtkStreamLine::New();
	streamer->SetInputData(data);
	streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetMaximumPropagationTime(200);
	streamer->SetIntegrationStepLength(.2);
	streamer->SetStepLength(.001);
	streamer->SetNumberOfThreads(1);
	streamer->SetIntegrationDirectionToBackward();//  IntegrateBothDirections();
	streamer->VorticityOn();
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
	// rake interactor
	//



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

	// line widget
	lineWidget->SetInteractor(iren);
	vtkCallbackCommand *callback = vtkCallbackCommand::New();
	callback->SetCallback(computeStreamlines);
	lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);

	//ren->AddActor(rakeActor);
	ren->AddActor(actor);
	ren->AddActor(outlineActor);
	ren->SetBackground(.5,.5,.5);

	renWin->SetSize(500,500);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}



