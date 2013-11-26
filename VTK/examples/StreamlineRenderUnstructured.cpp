// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <list>
#include <iterator>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

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
#include "vtkScalarBarActor.h"
#include "vtkLookupTable.h"
// streamline
#include "vtkStreamLine.h"
// vtu
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocatorInterpolatedVelocityField.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkDataSetMapper.h>

using namespace std;

vtkLineWidget *lineWidget;
vtkOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;
vtkPolyDataMapper *mapper;
vtkSmartPointer<vtkScalarBarActor> scalarBar;

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("compute\n");
	lineWidget->GetPolyData(seeds);
	streamer->Update();
	//mapper->SetScalarRange(streamer->GetOutput()->GetScalarRange());
	mapper->SetInputData(streamer->GetOutput());
	scalarBar->SetLookupTable(mapper->GetLookupTable());
	renWin->Render();

	double *range = streamer->GetOutput()->GetScalarRange();
	printf("Range: %lf %lf size=%lld\n", range[0], range[1], streamer->GetOutput()->GetPointData()->GetArray(0)->GetSize());
}


vtkDataSet *getData_vtu()
{
	vtkXMLUnstructuredGridReader *reader = vtkXMLUnstructuredGridReader::New();
	reader->SetFileName(SAMPLE_DATA_DIR "/unstructured/disk.vtu");
	reader->Update();

	vtkUnstructuredGrid *data = reader->GetOutput();
	if (data==0) {
		printf("data error\n");
	}
	data->Print(cout);

	data->GetPointData()->Print(cout);
	int r = data->GetPointData()->SetActiveAttribute("V", vtkDataSetAttributes::VECTORS);
	printf("index=%d\n", r);

	//reader->Delete();
	// test
	vtkCellLocatorInterpolatedVelocityField *interpolator = vtkCellLocatorInterpolatedVelocityField::New();
	interpolator->AddDataSet(data);
	double vel[3];
	double coords[4] = {0,0,0,0};
	for (int k=0; k<3; k++)
	{
		bool success = ((vtkFunctionSet *)interpolator)->FunctionValues(coords, vel);
		printf("querry success = %d, vec = %lf %lf %lf\n", success, vel[0], vel[1], vel[2]);
	}

	reader->Delete();
	return data;
}

int main(int argc, char **argv)
{
	printf("Press 'i' to change the rake\n");

	// read data
	vtkDataSet *data = getData_vtu();

	//
	// Determine seeds
	//
	// user can change the rake
	lineWidget = vtkLineWidget::New();
	lineWidget->SetInputData(data);
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
	streamer->SetInputData(data);
	streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumPropagationTime(200);
	streamer->SetNumberOfThreads(1);
	streamer->VorticityOn();
	streamer->Update();

	mapper = vtkPolyDataMapper::New();
	mapper->SetInputConnection(streamer->GetOutputPort());
	vtkActor *actor = vtkActor::New();
	actor->SetMapper(mapper);


	//
	// outline
	//
	vtkActor *outlineActor;
	if (data->IsA("vtkUnstructuredGrid")) {
		vtkUnstructuredGridGeometryFilter *outline = vtkUnstructuredGridGeometryFilter::New();
		outline->SetInputData(data);

		vtkDataSetMapper *outlineMapper = vtkDataSetMapper::New();
		outlineMapper->SetInputConnection(outline->GetOutputPort());

		outlineActor = vtkActor::New();
		outlineActor->SetMapper(outlineMapper);
		outlineActor->GetProperty()->SetColor(0,0,1.0);
		outlineActor->GetProperty()->SetOpacity(.2);
	} else
	{
		vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
		outline->SetInputData(data);

		vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
		outlineMapper->SetInputConnection(outline->GetOutputPort());

		outlineActor = vtkActor::New();
		outlineActor->SetMapper(outlineMapper);
		outlineActor->GetProperty()->SetColor(0,0,0);
	}



	// scalar bar
	  scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();
	  scalarBar->SetLookupTable(mapper->GetLookupTable());
	  scalarBar->SetTitle("Curvature");
	  scalarBar->SetNumberOfLabels(4);


	  // Create a lookup table to share between the mapper and the scalarbar
	  vtkSmartPointer<vtkLookupTable> hueLut =
	    vtkSmartPointer<vtkLookupTable>::New();
	  hueLut->SetTableRange (streamer->GetOutput()->GetScalarRange());
	  hueLut->SetHueRange (.7, 0);
	  hueLut->SetSaturationRange (1, 1);
	  hueLut->SetValueRange (1, 1);
	  hueLut->Build();

	  mapper->SetLookupTable( hueLut );
	  scalarBar->SetLookupTable( hueLut );


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
	//ren->AddActor(scalarBar);
	ren->SetBackground(.5,.5,.5);
	//ren->SetBackground(1,1,1);

	renWin->SetSize(500,500);

	iren->Initialize();
	renWin->Render();
	iren->Start();

	return 0;
}



