// An example detecting vortices in curvilinear grid data
// Samples regularly and compute the four criteria
// Stores into raw files with nrrd headers
// By Chun-Ming Chen

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vtkDataSet.h"
#include "vtkStructuredGrid.h"
#include "vtkMultiBlockPLOT3DReader.h"
#include "vtkMultiBlockDataSet.h"
#include <vtkInterpolatedVelocityField.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkImageData.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkFunctionSet.h>
#include <vtkImageInterpolator.h>
#include <vtkPointData.h>
#include <vtkImageData.h>

#include "VectorFieldVTK.h"
#include "OSUFlow.h"

int main(int argc, char ** argv)
{
#if 0
    // read PLOT3D data
    char file1[256], file2[256];
    int files;
    if (argc<=1) { // load default data
        sprintf(file1, "%s/curvilinear/combxyz.bin", SAMPLE_DATA_DIR); //t->GetDataRoot());
        printf("%s\n", file1);
        sprintf(file2, "%s/curvilinear/combq.bin", SAMPLE_DATA_DIR); //t->GetDataRoot());
        files = 2;
    } else {
        strcpy(file1, argv[1]);
        strcpy(file2, argv[2]);
    }

    // Start by loading some data.
    vtkMultiBlockPLOT3DReader *pl3dReader = vtkMultiBlockPLOT3DReader::New();
    // set data
    pl3dReader->SetXYZFileName(file1);
    pl3dReader->SetQFileName(file2);
    pl3dReader->SetAutoDetectFormat(1);  // should be on for loading binary file
    //pl3dReader->SetScalarFunctionNumber(100);
    pl3dReader->SetVectorFunctionNumber(200); // load velocity
    pl3dReader->Update();
    vtkDataSet *data = vtkDataSet::SafeDownCast( pl3dReader->GetOutput()->GetBlock(0) );

    OSUFlow *osuflow = new OSUFlow;
    CVectorField *field = new VectorFieldVTK( data );
    osuflow->SetFlowField( field );

#else
    char file[256];
    if (argc<=1) {
      sprintf(file, "%s/regular/tornado/1.vec", SAMPLE_DATA_DIR);
    } else {
      strcpy(file, argv[1]);
    }

    // debug with regular grids
    OSUFlow *osuflow = new OSUFlow;
    osuflow->LoadData(file, true); //true: static dataset
    CVectorField *field = osuflow->GetFlowField();
#endif

    //field->NormalizeField(true);


    VECTOR3 minLen, maxLen;
    osuflow->Boundary(minLen, maxLen);
    printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n",
                                minLen[0], maxLen[0], minLen[1], maxLen[1],
                                minLen[2], maxLen[2]);


    int dim[3]	;
    field->getDimension(dim[0], dim[1], dim[2]);
    printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);

    // define range
    int i1, i2, j1, j2, k1, k2;
    int i,j,k;
    i1 = 0; i2 = dim[0];
    j1 = 0; j2 = dim[1];
    k1 = 0; k2 = dim[2];

    // init file
    char out_fname[4][256];
    char prefix[4][10]={"lambda2", "q", "delta", "gamma"};
    for (i=0; i<4; i++)
    {
        sprintf(out_fname[i], "vortex_%s.raw", prefix[i]);
        printf("Output file: %s\n", out_fname[i]);
    }

    // determine unit
    float x,y,z;
    float unit = std::min(std::min(maxLen[0]-minLen[0], maxLen[1]-minLen[1]), maxLen[2]-minLen[2]) / 50;
    float delta = unit * .1f; // for Jacobian computation
    printf("Sampling unit: %f\n", unit);

    // open file
    int count=0;
    FILE *fp[4];
    for (i=0; i<4; i++)
        fp[i] = fopen(out_fname[i], "wb");

    // detect vortices
    for (z=minLen[2]; z< maxLen[2]; z+=unit) {
        for (y=minLen[1]; y< maxLen[1]; y+=unit)
            for (x=minLen[0]; x< maxLen[0]; x+=unit)
            {
#if 0
                //VECTOR3 v;
                //osuflow->GetFlowField()->at_phys(VECTOR3(x,y,z), 0, v);
                //printf("%f %f %f\n", v[0], v[1], v[2]);
#endif
                float f[4]; //lambda2, q, delta, gamma;
                field->GenerateVortexMetrics(VECTOR3(x,y,z), f[0], f[1], f[2], f[3], delta);
                for (i=0; i<4; i++)
                    fwrite((char *)&f[i], 1, 4, fp[i]);
                count++;
            }
        printf("z=%f\n", z);
    }
    for (i=0; i<4; i++)
        fclose(fp[i]);


    // get sampling dim
    int bdim[3];
    bdim[0] = bdim[1] = bdim[2] = 0;
    {
        int d;
        for (d=0; d<3; d++)
            for (x=minLen[d]; x<=maxLen[d]; x+=unit)
                bdim[d] ++;
    }

    // output
    for (i=0; i<4; i++)
    {
        // get out_fname filename only (no path)
        char *out_fname_no_path = strrchr(out_fname[i], '/');
        if (out_fname_no_path==NULL) out_fname_no_path = out_fname[i]; else out_fname_no_path++;

        char out_desc_fname[256];
        sprintf(out_desc_fname, "%s_%s.nhdr", "vortex", prefix[i]);
        FILE *fp = fopen(out_desc_fname, "wt");
        fprintf(fp,
                "NRRD0001\n"
                "type: float\n"
                "dimension: 3\n"
                "sizes: %d %d %d\n"
                "encoding: raw\n"
                "data file: %s\n"
                "space origin: (%f,%f,%f)\n"
                "space directions: (%f,0,0) (0,%f,0) (0,0,%f)\n"
                "# sampling distance: %f\n"
                "# grid range: %d %d %d - %d %d %d\n"
                "# physical range: %f %f %f - %f %f %f\n",
                bdim[0], bdim[1], bdim[2], out_fname_no_path,
                minLen[0], minLen[1], minLen[2], unit, unit, unit,
                unit, i1, j1, k1, i2, j2, k2, minLen[0], minLen[1], minLen[2], maxLen[0], maxLen[1], maxLen[2]);
        fclose(fp);
    }

    printf("Done (%d elems)\n", count);
    return 0;
}
