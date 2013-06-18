#ifndef VTK_P_OSUFLOW
#define VTK_P_OSUFLOW

#include <vector>
#include <vtkSmartPointer.h>

#include "diy.h"
#include "vtkOSUFlow.h"


class vtkInformation;
class vtkInformationVector;
class vtkTableExtentTranslator;

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

	vtkTableExtentTranslator *extentTable;
	void getNeighborIds(std::vector<gb_t> &neighborIdAry, vtkExtentTranslator *translator, int rank);
	void initExtentTable(vtkExtentTranslator *translator);

private:

	// Not implemented:
	//vtkPOSUFlow(const vtkPOSUFlow&);
	//void operator=(const vtkPOSUFlow&);
};


#endif //VTK_OSUFLOW
