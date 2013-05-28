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


	inline void SetRandomSeedPoints(float bMin[3], float bMax[3], int num_seeds) {
		this->SetNumberOfInputPorts(0);
		osuflow->SetRandomSeedPoints(bMin, bMax, num_seeds);
	}

	// compatibility
	inline OSUFlowVTK *getOSUFlow() { return osuflow; }

	virtual int FillInputPortInformation(int port, vtkInformation *info);

protected:
	vtkOSUFlow();
	~vtkOSUFlow() ;

	// Convert streamer array into vtkPolyData
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	// Not implemented:
	//vtkOSUFlow(const vtkOSUFlow&);
	//void operator=(const vtkOSUFlow&);
};


#endif //VTK_OSUFLOW
