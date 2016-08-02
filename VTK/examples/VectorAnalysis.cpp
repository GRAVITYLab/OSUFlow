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

void compute_vorticity(CVectorField *field)
{
  int dim[3]	;
  field->getDimension(dim[0], dim[1], dim[2]);
  printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);

  vector<VECTOR3> vort(dim[0]*dim[1]*dim[2]);

  field->GenerateVortField(0, false, &vort[0]);

  // init file
  char out_fname[256];
  sprintf(out_fname, "output_vortmag.raw");
  printf("Output file: %s\n", out_fname);

  // open file
  FILE *fp;
  fp = fopen(out_fname, "wb");

  int i,j,k;
  int count=0;
  for (k=0; k<dim[2]; k++)
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
      {
        VECTOR3 vec = vort[count++];
        float f = vec.GetMag();
        fwrite(&f, 1, 4, fp);
      }
  fclose(fp);

}

void vortex_analysis(CVectorField *field)
{
  int i,j,k;
  int dim[3]	;
  field->getDimension(dim[0], dim[1], dim[2]);
  printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);

  // init file
  const int MEASURES = 8;
  char out_fname[MEASURES][256];
  char prefix[MEASURES][20]={"lambda2", "q", "delta", "gamma", "divergence", "mag", "curvature", "JacDet"};
  for (i=0; i<MEASURES; i++)
  {
      sprintf(out_fname[i], "output_%s.raw", prefix[i]);
      printf("Output file: %s\n", out_fname[i]);
  }

  // open file
  FILE *fp[MEASURES];
  for (i=0; i<MEASURES; i++)
      fp[i] = fopen(out_fname[i], "wb");

  // detect vortices
  for (k=0; k<dim[2]; k++)
  {
    for (j=0; j<dim[1]; j++)
      for (i=0; i<dim[0]; i++)
      {
        float f[MEASURES]; //lambda2, q, delta, gamma, divergence, mag;
        field->GenerateVortexMetrics(i,j,k, f[0], f[1], f[2], f[3]);

        MATRIX3 J = field->JacobianStructuredGrid(i,j,k);

        // divergence
        f[4] = J(0,0) + J(1,1) + J(2,2);

        // velocity mag
        VECTOR3 vec;
        field->at_vert(i,j,k,0, vec);
        f[5] = vec.GetMag();

        // curvature
        f[6] = field->Curvature(i,j,k);

        // Jacobian det
        f[7] = J.det();

        for (int q=0; q<MEASURES; q++)
          fwrite(&f[q], 1, 4, fp[q]);

      }
    //printf("=%f\n", k);
  }
  for (i=0; i<MEASURES; i++)
      fclose(fp[i]);

  // output
  for (i=0; i<MEASURES; i++)
  {
    // get out_fname filename only (no path)
    char *out_fname_no_path = strrchr(out_fname[i], '/');
    if (out_fname_no_path==NULL) out_fname_no_path = out_fname[i]; else out_fname_no_path++;

    char out_desc_fname[256];
    sprintf(out_desc_fname, "%s_%s.nhdr", "output", prefix[i]);
    FILE *fp = fopen(out_desc_fname, "wt");
    fprintf(fp,
            "NRRD0001\n"
            "type: float\n"
            "dimension: 3\n"
            "sizes: %d %d %d\n"
            "encoding: raw\n"
            "data file: %s\n",
            dim[0], dim[1], dim[2], out_fname_no_path);
    fclose(fp);
  }

}

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

    compute_vorticity(field);

    vortex_analysis(field);



    printf("Done \n");
    return 0;
}
