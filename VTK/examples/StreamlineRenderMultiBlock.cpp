// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <iostream>
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
#include "vtkXMLMultiBlockDataReader.h"
// streamline
#include "vtkStreamLine.h"

using namespace std;

vtkLineWidget *lineWidget;
vtkOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("compute\n");
	lineWidget->GetPolyData(seeds);
	renWin->Render();
	streamer->Update();
}


int main(int argc, char **argv)
{
	printf("Press 'i' to change the rake\n");

	// read data
	char file1[256];
	strcpy(file1, argv[1]);

	// Start by loading some data.
	vtkXMLMultiBlockDataReader *reader = vtkXMLMultiBlockDataReader::New();
	// set data
    reader->SetFileName(file1);
	reader->Update();
	vtkMultiBlockDataSet *data;

    vtkMultiBlockDataSet *mbdata = vtkMultiBlockDataSet::SafeDownCast(reader->GetOutput());
    if (mbdata==NULL)
        printf("First block not extracted\n");

    vtkDataSet *data0 = vtkDataSet::SafeDownCast(mbdata->GetBlock(0));

    //reader->PrintSelf(cout, vtkIndent(0));

	//
	// Determine seeds
	//
	// user can change the rake
	lineWidget = vtkLineWidget::New();
    lineWidget->SetInputData(data0);
	lineWidget->SetResolution(21); // 22 seeds along the line
	lineWidget->SetAlignToYAxis();
	lineWidget->PlaceWidget();
	lineWidget->ClampToBoundsOn();
	seeds = vtkPolyData::New();
	lineWidget->GetPolyData(seeds);


	//
	// vtkOSUFlow
	//
	streamer = vtkOSUFlow::New();
    streamer->SetInputData(mbdata);
	streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumPropagationTime(200);
	streamer->SetNumberOfThreads(1);
	streamer->VorticityOn();
	//streamer->getOSUFlow()->initOpenMP(8); // set number of processes to 8

	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(streamer->GetOutputPort());
	mapper->SetScalarRange(data0->GetScalarRange());
	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);


	//
	// outline
	//
	vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
    outline->SetInputData(data0);

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



