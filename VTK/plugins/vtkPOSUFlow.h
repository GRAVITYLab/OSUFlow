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

	vtkSetMacro(UseDIYPartition, bool);
	vtkGetMacro(UseDIYPartition, bool);

	//vtkSetClampMacro(MaxRounds, int, 0, VTK_INT_MAX);
	vtkSetMacro(MaxRounds, int);
	vtkGetMacro(MaxRounds, int);
protected:
	vtkPOSUFlow();
	~vtkPOSUFlow() ;

	virtual int RequestUpdateExtent(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	vtkTableExtentTranslator *extentTable;
	void getNeighborIds(std::vector<gb_t> &neighborIdAry, vtkExtentTranslator *translator, int rank);
	void initExtentTable(vtkExtentTranslator *translator);
	void initExtentTableByDIY(vtkExtentTranslator *translator);

protected:
	int MaxRounds;	// max number of rounds (seed interchanges)
	bool UseDIYPartition; // partition the data by diy or by the input extent
	float WaitFactor; // wait factor for nonblocking communication
	                // wait for this portion of messages to arrive each round
	bool diy_initialized;

	double totTime;
	double totInTime;
	double totCompCommTime;	// computation and communication time
	double totOutTime; 		// gathering time
private:

	// Not implemented:
	//vtkPOSUFlow(const vtkPOSUFlow&);
	//void operator=(const vtkPOSUFlow&);
};


#endif //VTK_OSUFLOW
