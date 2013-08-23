
#include <stdio.h>
#include <stdlib.h>

#include <list>
#include <iterator>

#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkPLOT3DReader.h"
#include "vtkSmartPointer.h"
#include "OSUFlowVTK.h"

using namespace std;

#define VTK_DATA_ROOT "/home/local/KHQ/chunming.chen/project/VTKData"


int main()
{
	// Start by loading some data.
	vtkPLOT3DReader *pl3dReader = vtkPLOT3DReader::New();
	pl3dReader->SetXYZFileName(VTK_DATA_ROOT  "/Data/combxyz.bin");
	pl3dReader->SetQFileName(VTK_DATA_ROOT  "/Data/combq.bin");
	pl3dReader->SetScalarFunctionNumber(100);
	pl3dReader->SetVectorFunctionNumber(202);
	pl3dReader->Update();

	// random points
	//vtkStructuredGrid *grid = pl3dReader->GetOutput();
	//int *dim = grid->GetDimensions();
	double *bounds = pl3dReader->GetOutput()->GetBounds();
	printf("bounds: %lf %lf %lf %lf %lf %lf\n", bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	OSUFlowVTK *osuflow = new OSUFlowVTK;
	vtkSmartPointer<vtkDataSet> sData = pl3dReader->GetOutput();
	osuflow->setData(sData);

	// gen seeds
	float from[3], to[3];
	from[0]= bounds[0]; to[0] = bounds[1];
	from[1]= bounds[2]; to[1] = bounds[3];
	from[2]= bounds[4]; to[2] = bounds[5];
	osuflow->SetRandomSeedPoints(from, to, 100);

	// init results
	vtk
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



