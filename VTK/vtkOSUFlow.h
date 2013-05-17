#ifndef VTK_OSUFLOW
#define VTK_OSUFLOW

#include <assert.h>

#include <vtkStreamer.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>
#include <vtkInformationVector.h>

#include <OSUFlowVTK.h>
#include <Field.h>


class vtkDataSet;

class vtkOSUFlow: public vtkStreamer   //vtkPolyDataAlgorithm
{
	vtkTypeMacro(vtkOSUFlow,vtkStreamer);
	OSUFlowVTK *osuflow;

public:

	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// Construct object with step size set to 1.0.
	static vtkOSUFlow *New();


	// Description:
	// Specify the length of a line segment. The length is expressed in terms of
	// elapsed time. Smaller values result in smoother appearing streamlines, but
	// greater numbers of line primitives.
	vtkSetClampMacro(StepLength,double,0.000001,VTK_DOUBLE_MAX);
	vtkGetMacro(StepLength,double);

	inline void SetRandomSeedPoints(float bMin[3], float bMax[3], int num_seeds) {
		this->SetNumberOfInputPorts(0);
		osuflow->SetRandomSeedPoints(bMin, bMax, num_seeds);
	}
protected:
	vtkOSUFlow();
	~vtkOSUFlow() {};

	// Convert streamer array into vtkPolyData
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	// the length of line primitives
	double StepLength;

private:
	// Not implemented:
	//vtkOSUFlow(const vtkOSUFlow&);
	//void operator=(const vtkOSUFlow&);
};


#endif //VTK_OSUFLOW
