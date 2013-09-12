
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
#include "vtkLineWidget.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkStructuredGrid.h"
#include "vtkOutlineFilter.h"
#include "vtkPlaneWidget.h"
// streamline
#include "vtkStreamLine.h"

using namespace std;

vtkLineWidget *lineWidget;
vtkPlaneWidget *planeWidget;
vtkOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds, *seeds2 ;

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	double *point1 = lineWidget->GetPoint1();
	double *point2 = lineWidget->GetPoint2();
	printf("LineWidget Point1: %lf %lf %lf, Point2: %lf %lf %lf\n", point1[0], point1[1], point1[2], point2[0], point2[1], point2[2]);

	lineWidget->GetPolyData(seeds);
	renWin->Render();
}
void computeStreamlines2(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("2\n");
	planeWidget->GetPolyData(seeds);
	renWin->Render();
	planeWidget->On();
}


int main(int argc, char **argv)
{
	printf("Press 'i' to change the rake\n");

	streamer = vtkOSUFlow::New();

	// read data
	OSUFlow *osuflow = streamer->getOSUFlow();
	const char *filename;
	if (argc>1)
		filename = argv[1];
	else
		filename = SAMPLE_DATA_DIR "/regular/tornado/1.vec";
	osuflow->LoadData(filename, true); //true: static dataset

#ifdef _OPENMP
	osuflow->initOpenMP(8);
#endif

	// dummy dataset for boundary
	VECTOR3 minB, maxB;
	osuflow->Boundary(minB, maxB);

	vtkPolyData *data = vtkPolyData::New();
	vtkPoints *points = vtkPoints::New();
	float p[3];
	p[0] = minB[0]; p[1] = minB[1]; p[2] = minB[2];	points->InsertNextPoint(p);
	p[0] = minB[0]; p[1] = minB[1]; p[2] = maxB[2];	points->InsertNextPoint(p);
	p[0] = minB[0]; p[1] = maxB[1]; p[2] = minB[2];	points->InsertNextPoint(p);
	p[0] = minB[0]; p[1] = maxB[1]; p[2] = maxB[2];	points->InsertNextPoint(p);
	p[0] = maxB[0]; p[1] = minB[1]; p[2] = minB[2];	points->InsertNextPoint(p);
	p[0] = maxB[0]; p[1] = minB[1]; p[2] = maxB[2];	points->InsertNextPoint(p);
	p[0] = maxB[0]; p[1] = maxB[1]; p[2] = minB[2];	points->InsertNextPoint(p);
	p[0] = maxB[0]; p[1] = maxB[1]; p[2] = maxB[2];	points->InsertNextPoint(p);
	data->SetPoints(points);
	points->Delete();


	//
	// Determine seeds
	//
	// rake
	lineWidget = vtkLineWidget::New();
	lineWidget->SetInputData(data);
	lineWidget->SetResolution(21); // 22 seeds along the line
	lineWidget->SetAlignToYAxis();
	lineWidget->PlaceWidget();
	lineWidget->ClampToBoundsOn();
	seeds = vtkPolyData::New();
	lineWidget->GetPolyData(seeds);

#if 0
	// plane
	planeWidget = vtkPlaneWidget::New();
	planeWidget->SetInputData(data);
	planeWidget->SetResolution(5);
	planeWidget->PlaceWidget();
	planeWidget->SetKeyPressActivationValue('j');
	seeds2 = vtkPolyData::New();
	planeWidget->GetPolyData(seeds2);
#endif

	//
	// vtkOSUFlow
	//
	streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumPropagationTime(1000);
	streamer->SetNumberOfThreads(1);
	streamer->VorticityOn();

	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(streamer->GetOutputPort());
	mapper->SetScalarRange(data->GetScalarRange());
	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);


	//
	// outline
	//
	vtkOutlineFilter *outline = vtkOutlineFilter ::New();
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
	renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);
	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	// line widget interactor
	lineWidget->SetInteractor(iren);
	lineWidget->SetDefaultRenderer(ren);
	vtkCallbackCommand *callback = vtkCallbackCommand::New();
	callback->SetCallback(computeStreamlines);
	lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);

#if 0
	// plane widget interactor
	planeWidget->SetInteractor(iren);
	planeWidget->SetDefaultRenderer(ren);
	vtkCallbackCommand *callback2 = vtkCallbackCommand::New();
	callback2->SetCallback(computeStreamlines2);
	planeWidget->AddObserver(vtkCommand::EndInteractionEvent, callback2);
#endif


	ren->AddActor(actor);
	ren->AddActor(outlineActor);
	ren->SetBackground(.5,.5,.5);

	renWin->SetSize(500,500);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}



