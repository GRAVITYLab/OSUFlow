#ifndef VTK_OSUFLOW
#define VTK_OSUFLOW

#include <assert.h>

#include <vtkStreamer.h>
#include <vtkSmartPointer.h>

#include <OSUFlow.h>
#include <Field.h>

class vtkDataSet;
class vtcInformationVector;
class vtkInformation;

class vtkOSUFlow: public vtkStreamer   //vtkPolyDataAlgorithm
{
	vtkTypeMacro(vtkOSUFlow,vtkStreamer);

public:
	static vtkOSUFlow *New();

	// Description
	// Use OSUFlow's seeding function
	inline void SetRandomSeedPoints(float bMin[3], float bMax[3], int num_seeds) {
		osuflow->SetRandomSeedPoints(bMin, bMax, num_seeds);
	}

	inline OSUFlow *getOSUFlow() { return osuflow; }

	// Description
	// Tells VTK pipeline that both input ports are optional (Data can be assigned to OSUFlow in advance)
	// The number of ports (2) are assigned in the super class vtkStreamer
	virtual int FillInputPortInformation(int port, vtkInformation *info);

	void PrintSelf(ostream& os, vtkIndent indent);


	// Set/get variables
	vtkSetMacro(IntegratorOrder, int);
	vtkGetMacro(IntegratorOrder, int);

	vtkSetClampMacro(MinimumIntegrationStep,double,0.0000001,VTK_DOUBLE_MAX);
	vtkGetMacro(MinimumIntegrationStep,double);

	vtkSetClampMacro(MaximumIntegrationStep,double,0.0000001,VTK_DOUBLE_MAX);
	vtkGetMacro(MaximumIntegrationStep,double);

	vtkSetMacro(MaximumError, double);
	vtkGetMacro(MaximumError, double);

	vtkSetMacro(MaximumNumberOfSteps, vtkIdType);
	vtkGetMacro(MaximumNumberOfSteps, vtkIdType);

	vtkSetMacro(scale, double);
	vtkGetMacro(scale, double);

protected:
	double MinimumIntegrationStep;
	double MaximumIntegrationStep;

	// Description
	// Used in RK45
	double MaximumError;

	// Description
	// Specify the maximum number of steps for integrating a streamline.
	int MaximumNumberOfSteps;

	// Description
	// RK2, RK4, RK45 (INTEG_ORD)
	int IntegratorOrder;

	double scale;

	vtkOSUFlow();
	~vtkOSUFlow() ;

	// Description
	// Convert streamer array into vtkPolyData
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	OSUFlow *osuflow;

	// Not implemented:
	//vtkOSUFlow(const vtkOSUFlow&);
	//void operator=(const vtkOSUFlow&);
};


#endif //VTK_OSUFLOW
