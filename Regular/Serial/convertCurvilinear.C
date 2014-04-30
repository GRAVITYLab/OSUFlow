#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "OSUFlow.h"
#include "CurvilinearGrid.h"
#include "vtk.h"

#define unit 0.01f
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>

using namespace std;

float *read_txt(const char *filename)
{
	int startX,endX,startY,endY,startZ,endZ;
	FILE *fp = fopen(filename, "rb");

	fscanf(fp, "%d %d %d %d %d %d", &endX, &startX, &endY, &startY, &endZ, &startZ);
	float *data = new float[(endX-startX) * (endY-startY) * (endZ-startZ)];

	for(int k=startZ; k<endZ; k++)
		for(int j=startY; j<endY; j++)
			for(int i=startX; i<endX; i++)
			{
				fscanf("%f", f);
			}

	fclose(fp);
	return data;
}

int main(int argc, char ** argv)
{
	prirntf("Usage: [plot3d file prefix] [labeled TXT file]\n");

	OSUFlow *osuflow = new OSUFlow();
	printf("read file %s\n", argv[1]);

	VECTOR3 minLen, maxLen, minB, maxB;
	minB[0] = 0; minB[1] = 0; minB[2] = 0;
    maxB[0] = 200; maxB[1] = 200; maxB[2] = 200;

	osuflow->LoadDataCurvilinear((const char*)argv[1], true, minB, maxB); //true: a steady flow field
	osuflow->Boundary(minLen, maxLen);
	printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n",
								minLen[0], maxLen[0], minLen[1], maxLen[1],
								minLen[2], maxLen[2]);


	int dim[3]	;
	CurvilinearGrid *grid = (CurvilinearGrid*)osuflow->GetFlowField()->GetGrid();
	grid->GetDimension(dim[0], dim[1], dim[2]);
	printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);
	int i,j,k;
#if 0
	for (k=0; k<dim[2]; k++)
		for (j=0; j<dim[1]; j++)
			for (i=0; i<dim[0]; i++) {
				VECTOR3 v;
				grid->coordinates_at_vertex(VECTOR3(i,j,k),&v);	//return 1 velocity value
				printf("%f %f %f\n", v[0], v[1], v[2]);
			}
#endif

	int i1, i2, j1, j2, k1, k2;
	if (argc<8) {
		printf("Please input range: (i2 i1 j2 j1 k2 k1)\n");
		scanf("%d %d %d %d %d %d", &i2, &i1, &j2, &j1, &k2, &k1);
	} else {
		i2 = atoi(argv[2]);
		i1 = atoi(argv[3]);
		j2 = atoi(argv[4]);
		j1 = atoi(argv[5]);
		k2 = atoi(argv[6]);
		k1 = atoi(argv[7]);
	}
	// Does not include upbound


	char out_fname[4][256];
	char prefix[4][10]={"lambda2", "q", "delta", "gamma"};
	for (i=0; i<4; i++)
	{
		sprintf(out_fname[i], "%s_%s.raw", argv[8], prefix[i]);
		printf("Output file: %s\n", out_fname[i]);
	}

	float x,y,z;
	VECTOR3 from, to;
	grid->coordinates_at_vertex(VECTOR3(i1, j1, k1), &from);
	grid->coordinates_at_vertex(VECTOR3(i2, j2, k2), &to);
	printf("from: %f %f %f, to: %f %f %f\n", from[0], from[1], from[2], to[0], to[1], to[2]);

	// get min offset unit
	float min_off[3] =  {1e+9,1e+9,1e+9};
	{
		for (i=i1+1; i<i2; i++) {
			VECTOR3 v1, v2;
			grid->coordinates_at_vertex(VECTOR3(i-1, j1, k1), &v1);
			grid->coordinates_at_vertex(VECTOR3(i, j1, k1), &v2);
			min_off[0] = min(v2[0]-v1[0], min_off[0]);
		}
		for (j=j1+1; j<j2; j++) {
			VECTOR3 v1, v2;
			grid->coordinates_at_vertex(VECTOR3(i1, j-1, k1), &v1);
			grid->coordinates_at_vertex(VECTOR3(i1, j, k1), &v2);
			min_off[1] = min(v2[1]-v1[1], min_off[1]);
		}
		for (k=k1+1; k<k2; k++) {
			VECTOR3 v1, v2;
			grid->coordinates_at_vertex(VECTOR3(i1, j1, k-1), &v1);
			grid->coordinates_at_vertex(VECTOR3(i1, j1, k), &v2);
			min_off[2] = min(v2[2]-v1[2], min_off[2]);
		}
		printf("Min grid unit: %f %f %f\n", min_off[0], min_off[1], min_off[2]);
	}

	if (!(from[0]<=to[0] && from[1]<=to[1] && from[2]<=to[2]))
		printf("Input invalid.  Program halts\n");
	int count=0;
	FILE *fp[4];
	for (i=0; i<4; i++)
		fp[i] = fopen(out_fname[i], "wb");

	for (z=from[2]; z< to[2]; z+=unit) {
		for (y=from[1]; y< to[1]; y+=unit)
			for (x=from[0]; x< to[0]; x+=unit)
			{
#if 0
				//VECTOR3 v;
				//osuflow->GetFlowField()->at_phys(VECTOR3(x,y,z), 0, v);
				//printf("%f %f %f\n", v[0], v[1], v[2]);
#endif
				float f[4]; //lambda2, q, delta, gamma;
				osuflow->GetFlowField()->GenerateVortexMetrics(VECTOR3(x,y,z), f[0], f[1], f[2], f[3]);
				for (i=0; i<4; i++)
					fwrite((char *)&f[i], 1, 4, fp[i]);
				count++;
			}
		printf("z=%f\n", z);
	}
	for (i=0; i<4; i++)
		fclose(fp[i]);


	// get dim for given range
	int bdim[3];
	bdim[0] = bdim[1] = bdim[2] = 0;
	{
		int d;
		for (d=0; d<3; d++)
			for (x=from[d]; x<=to[d]; x+=unit)
				bdim[d] ++;
	}

	for (i=0; i<4; i++)
	{
		// get out_fname filename only (no path)
		char *out_fname_no_path = strrchr(out_fname[i], '/');
		if (out_fname_no_path==NULL) out_fname_no_path = out_fname[i]; else out_fname_no_path++;

		char out_desc_fname[256];
		sprintf(out_desc_fname, "%s_%s.nhdr", argv[8], prefix[i]);
		FILE *fp = fopen(out_desc_fname, "wt");
		fprintf(fp,
				"NRRD0001\n"
				"type: float\n"
				"dimension: 3\n"
				"sizes: %d %d %d\n"
				"encoding: raw\n"
				"data file: %s\n"
				"# sampling distance: %f\n"
				"# grid range: %d %d %d - %d %d %d\n"
				"# physical range: %f %f %f - %f %f %f\n"
				"# min grid unit: %f %f %f\n",
				bdim[0], bdim[1], bdim[2], out_fname_no_path,
				unit, i1, j1, k1, i2, j2, k2, from[0], from[1], from[2], to[0], to[1], to[2], min_off[0], min_off[1], min_off[2]);
		fclose(fp);
	}

	printf("Done (%d elems)\n", count);

	return 0;
}
