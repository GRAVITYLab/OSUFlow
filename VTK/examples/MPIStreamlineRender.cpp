// Reference: VTK/Examples/VisualizationAlgorithms/Python/StreamlinesWithLineWidget.py
#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

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
#include <vtkXMLImageDataReader.h>
#include <vtkOutlineFilter.h>
#include <vtkPStreamTracer.h>
#include <vtkImageData.h>
// MPI
#include <vtkXMLPImageDataReader.h>
#include <vtkMPIController.h>
#include <vtkCompositeRenderManager.h>
#include <vtkStreamingDemandDrivenPipeline.h>

// streamline
#include "vtkPOSUFlow.h"
//#include "vtkPStreamTracer.h"

using namespace std;

vtkLineWidget *lineWidget;
vtkPOSUFlow *streamer;
vtkRenderWindow *renWin;
vtkPolyData *seeds ;


struct Args
{
  int argc;
  char** argv;
};

void computeStreamlines(vtkObject* caller, unsigned long eventId, void *clientdata, void *calldata)
{
	printf("compute\n");
	lineWidget->GetPolyData(seeds);
	renWin->Render();
	streamer->Update();
}



// This will be called by all processes
void MyMain( vtkMultiProcessController *controller, void *arg )
{
	int rank = controller->GetLocalProcessId();
	int nproc = controller->GetNumberOfProcesses();

	if (rank==0)
		printf("Press 'i' to change the rake\n");

	// read data
	Args *args = (Args *)arg;

	// Start by loading some data.
#if 1
	vtkXMLImageDataReader *reader = vtkXMLImageDataReader::New();
	// setdata
	reader->SetFileName(SAMPLE_DATA_DIR "/regular/tornado/1.vti");

	//reader->Update();
	//int *extent = reader->GetOutput()->GetExtent(); // vtkImageData::SafeDownCast(reader->GetOutputDataObject(0))->GetExtent();
	//printf("Extent: %d %d %d %d %d %d\n", extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);

#else
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	// set data
	pl3dReader->SetXYZFileName(file1);
	pl3dReader->SetQFileName(file2);
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();
	vtkDataSet *data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));
#endif

	//
	// Determine seeds
	//
#if 0
	// user can change the rake
	lineWidget = vtkLineWidget::New();
	lineWidget->SetInputConnection(reader->GetOutputPort());
	lineWidget->SetResolution(21); // 22 seeds along the line
	lineWidget->SetAlignToYAxis();
	lineWidget->PlaceWidget();
	lineWidget->ClampToBoundsOn();
	seeds = vtkPolyData::New();
	lineWidget->GetPolyData(seeds);
#else
	// use a static rake to generate streamlines
	vtkLineSource *rake = vtkLineSource::New();
	rake->SetPoint1(0,0,0); //extent[0], extent[2], extent[4]);
	rake->SetPoint2(47,47,47); //extent[1], extent[3], extent[5]);
	rake->SetResolution(100);
	vtkPolyDataMapper *rakeMapper = vtkPolyDataMapper::New();
	rakeMapper->SetInputConnection(rake->GetOutputPort());
	vtkActor *rakeActor = vtkActor::New();
	rakeActor->SetMapper(rakeMapper);
	rakeActor->VisibilityOn();
#endif


#if 0
	//
	// vtkOSUFlow
	//
	streamer = vtkPOSUFlow::New();
	streamer->SetInputConnection(0, reader->GetOutputPort());
	streamer->SetInputConnection(1, rake->GetOutputPort());
	//streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumIntegrationStep(100);
#else
	//
	// vtk StreamTracer
	//
	vtkPStreamTracer *streamer = vtkPStreamTracer::New();
	streamer->SetInputConnection(0, reader->GetOutputPort());
	streamer->SetInputConnection(1, rake->GetOutputPort());
	//streamer->SetSourceData(seeds);	//streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetIntegrationDirectionToForward();
	streamer->SetMaximumIntegrationStep(100);
#endif


	streamer->Update();
#if 0 // run more times
	rake->SetPoint1(0,0,0); //extent[0], extent[2], extent[4]);
	rake->SetPoint2(37,37,37); //extent[1], extent[3], extent[5]);
	streamer->Update();
	rake->SetPoint1(0,0,0); //extent[0], extent[2], extent[4]);
	rake->SetPoint2(47,47,47); //extent[1], extent[3], extent[5]);
	streamer->Update();
#endif

	if (rank==0)
	{
		vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
		mapper->SetInputConnection(streamer->GetOutputPort());
		vtkActor *actor = vtkActor::New();
		actor->SetMapper(mapper);


		//streamer->Update();

		// Tell the pipeline which piece we want to update.
		//vtkStreamingDemandDrivenPipeline* exec =
				//vtkStreamingDemandDrivenPipeline::SafeDownCast(streamer->GetExecutive());
		//exec->SetUpdateNumberOfPieces(exec->GetOutputInformation(0), nproc);
		//exec->SetUpdatePiece(exec->GetOutputInformation(0), rank);

		//
		// outline
		//

	#if 1
		vtkOutlineFilter *outline = vtkOutlineFilter::New();
		outline->SetInputConnection(0, reader->GetOutputPort());
	#else
		vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
		outline->SetInputConnection(0, reader->GetOutputPort());
	#endif

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

		// This class allows all processes to composite their images.
		// The root process then displays it in it's render window.
		//vtkCompositeRenderManager* tc = vtkCompositeRenderManager::New();
		//tc->SetRenderWindow(renWin);

		// line widget interactor
		//lineWidget->SetInteractor(iren);
		//lineWidget->SetDefaultRenderer(ren);
		//vtkCallbackCommand *callback = vtkCallbackCommand::New();
		//callback->SetCallback(computeStreamlines);
		//lineWidget->AddObserver(vtkCommand::EndInteractionEvent, callback);

		//ren->AddActor(rakeActor);
		ren->AddActor(actor);
		ren->AddActor(outlineActor);
		ren->SetBackground(.5,.5,.5);

		renWin->SetSize(500,500);
		iren->Initialize();
		renWin->Render();
		iren->Start();

		// Only the root process will have an active in	teractor. All
		// the other render windows will be slaved to the root.
		//tc->StartInteractor();
	} // rank==0

}



int main( int argc, char* argv[] )
{
  // This is here to avoid false leak messages from vtkDebugLeaks when
  // using mpich. It appears that the root process which spawns all the
  // main processes waits in MPI_Init() and calls exit() when
  // the others are done, causing apparent memory leaks for any objects
  // created before MPI_Init().
  MPI_Init(&argc, &argv);

  // Note that this will create a vtkMPIController if MPI
  // is configured, vtkThreadedController otherwise.
  vtkMPIController* controller = vtkMPIController::New();

  controller->Initialize(&argc, &argv, 1);

  Args args;
  args.argc = argc;
  args.argv = argv;

  controller->SetSingleMethod(MyMain, &args);
  controller->SingleMethodExecute();

  controller->Finalize();
  controller->Delete();

  return 0;
}


