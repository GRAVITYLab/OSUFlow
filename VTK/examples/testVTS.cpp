#include <iostream>

#include <assert.h>

#include <vtkInformation.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <Field.h>
#include "OSUFlow.h"
#include "VectorFieldVTK.h"


int main()
{
vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
reader->SetFileName(SAMPLE_DATA_DIR "/curvilinear/comb.vts");
reader->UpdateInformation();

//int extent[6] = {14, 28, 16, 32, 12, 24};
//reader->SetUpdateExtent(0, extent);
reader->Update();

reader->PrintSelf(std::cout, vtkIndent(2));

int *ext = reader->GetOutput()->GetExtent();
printf("Extent: %d %d %d %d %d %d\n", ext[0], ext[1], ext[2], ext[3], ext[4], ext[5]);
}


