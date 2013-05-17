
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
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkStructuredGridOutlineFilter.h"

using namespace std;


int main()
{
	// Start by loading some data.
	vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
	// set data
	{
		char path[256];
		vtkTesting *t = vtkTesting::New();
		sprintf(path, "%s/Data/combxyz.bin", t->GetDataRoot());
		printf("%s\n", path);
		pl3dReader->SetXYZFileName(path);
		sprintf(path, "%s/Data/combq.bin", t->GetDataRoot());
		pl3dReader->SetQFileName(path);
		t->Delete();
	}
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();
	vtkDataSet *data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));


	// use a rake to generate streamlines
	vtkLineSource *rake = vtkLineSource::New();
	rake->SetPoint1(15, -5, 32);
	rake->SetPoint2(15, 5, 32);
	rake->SetResolution(21);

	vtkPolyDataMapper *rakeMapper = vtkPolyDataMapper::New();
	rakeMapper->SetInputConnection(rake->GetOutputPort());

	vtkActor *rakeActor = vtkActor::New();
	rakeActor->SetMapper(rakeMapper);

	// vtkOSUFlow
	vtkOSUFlow *streamer = vtkOSUFlow::New();
	streamer->SetInputData(data);
	streamer->SetSourceConnection(rake->GetOutputPort());
	streamer->SetStepLength(0.001);
	streamer->SetIntegrationDirectionToIntegrateBothDirections();

#if 0
	// random points
	//vtkStructuredGrid *grid = pl3dReader->GetOutput();
	//int *dim = grid->GetDimensions();
	printf("Number of blocks=%d\n", pl3dReader->GetOutput()->GetNumberOfBlocks());
	vtkDataSet *data = vtkDataSet::SafeDownCast(pl3dReader->GetOutput()->GetBlock(0));
	double *bounds = data->GetBounds();
	printf("bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	// gen seeds
	float from[3], to[3];
	from[0]= bounds[0]; to[0] = bounds[1];
	from[1]= bounds[2]; to[1] = bounds[3];
	from[2]= bounds[4]; to[2] = bounds[5];
	osuflow->SetRandomSeedPoints(from, to, 100);
#endif

	vtkStructuredGridOutlineFilter *outline = vtkStructuredGridOutlineFilter::New();
	outline->SetInputData(data);

	vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
	outlineMapper->SetInputConnection(outline->GetOutputPort());

	vtkActor *outlineActor = vtkActor::New();
	outlineActor->GetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(0,0,0);

	vtkRenderer *ren = vtkRenderer::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren);
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	ren->AddActor(rakeActor);
	ren->AddActor(actor);
	ren->AddActor(outlineActor);
	ren->SetBackground(1,1,1);

	renWin->SetSize(300,300);

	iren->Initialize();
	renWin->Render();
	iren->Start();


}



