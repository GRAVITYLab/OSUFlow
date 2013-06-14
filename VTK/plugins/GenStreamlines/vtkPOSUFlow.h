#ifndef VTK_P_OSUFLOW
#define VTK_P_OSUFLOW

#include <assert.h>

#include <vtkStreamer.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInterpolatedVelocityField.h>
#include <vtkCellType.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>

#include "vtkOSUFlow.h"


class vtkPOSUFlow: public vtkOSUFlow //vtkPolyDataAlgorithm
{
	vtkTypeMacro(vtkPOSUFlow,vtkOSUFlow);

public:
	static vtkPOSUFlow *New();

	void PrintSelf(ostream& os, vtkIndent indent);

protected:

	vtkPOSUFlow();
	~vtkPOSUFlow() ;

	virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	// Not implemented:
	//vtkPOSUFlow(const vtkPOSUFlow&);
	//void operator=(const vtkPOSUFlow&);
};


#endif //VTK_OSUFLOW
