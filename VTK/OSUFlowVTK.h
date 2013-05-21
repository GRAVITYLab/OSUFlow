#ifndef OSUFLOW_VTK
#define OSUFLOW_VTK

#include <assert.h>

#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#if 0
#include <vtkXMLFileReadTester.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkMultiBlockDataGroupFilter.h>
#include <vtkGenericDataObjectReader.h>
#endif

#include <Field.h>
#include "OSUFlow.h"
#include "VectorFieldVTK.h"



class OSUFlowVTK: public OSUFlow
{
public:
	inline void setData(vtkDataSet* input)
	{
		double *bounds = input->GetBounds();

		gMin.Set(bounds[0], bounds[2], bounds[4]);
		gMax.Set(bounds[1], bounds[3], bounds[5]);
		lMin = gMin;
		lMax = gMax;
		MinT = 0;
		MaxT = 0;

		CVectorField* field = new VectorFieldVTK(input);
		flowField = field;
		has_data = true;
	}

	inline bool getHasData() { return this->has_data; }

#if 0
	vtkObject *loadVTKData(const char *fileName)
	{
	    vtkXMLFileReadTester *vtkXMLFormatFileChecker = vtkXMLFileReadTester::New();

	    vtkXMLFormatFileChecker->SetFileName(fileName);

		if (vtkXMLFormatFileChecker->TestReadFile() > 0) {
			//Logger.getLogger(getClass().getName()).log(Level.INFO, "[{0}] is of type [{1}]",
			//		new Object[] { fileName, vtkXMLFormatFileChecker.GetFileDataType() });

			vtkXMLGenericDataObjectReader *vtkXMLFileReader = vtkXMLGenericDataObjectReader::New();

			vtkXMLFileReader->SetFileName(fileName);
			vtkXMLFileReader->Update();

			if (vtkXMLFileReader->GetOutput()->IsA("vtkMultiBlockDataSet")) {
				return vtkXMLFileReader->GetOutput();
			} else {
				vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

				makeMultiblock->SetInputConnection(vtkXMLFileReader->GetOutputPort());

				return makeMultiblock;
			}
		} else // legacy format
		{
			vtkGenericDataObjectReader *legacyVTKFileReader = vtkGenericDataObjectReader::New();

			legacyVTKFileReader->SetFileName(fileName);
			legacyVTKFileReader->Update();

			vtkMultiBlockDataGroupFilter *makeMultiblock = vtkMultiBlockDataGroupFilter::New();

			makeMultiblock->SetInput(legacyVTKFileReader->GetOutput());

			return makeMultiblock;
		}
	}
#endif
};

#endif //OSUFLOW_VTK
