#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <hdf5.h>
#include "flashhdf5_float.h"
#include <map>
#include <vector>
#include <fstream>

using namespace std; 

void
help(char *name) {

  fprintf(stderr, "%s -f <inFile> -o <outFile> -b <blockFile> --x <xmin> --X <xmax> --y <ymin> --Y <ymax> --z <zmin> --Z <zmax> -L (-L for leaf nodes only) -p (probe file for domain size and other data)\n",name); 
  exit(-1);

}

int main(int argc, char **argv) {

  char inFile[180] = {'\0'};
  char outfile[180] = {'\0'};
  char blocksFileName[180] = {'\0'};
  char dnum[5];
  int nb, i,block_no;
  int cellDims[3];
  char velx[5] = {'\0'};
  char vely[5] = {'\0'};
  char velz[5] = {'\0'};
  float coords[3], bsize[3], bnd_box[6];
  float bounds[6];
  float *fdata1, *fdata2, *fdata3;
  int c;
  float x0,x1,y0,y1,z0,z1;
  std::vector<int> blk_list;
  bool leafsonly = false;
  int nxb,nyb,nzb;
  int nb_out;
  bool probe;
  float min,max;

  strcpy(velx,"velx");
  strcpy(vely,"vely");
  strcpy(velz,"velz");
  probe = false;
  min = 9999999999999999999999.0f;
  max = -9999999999999999999999.0f;
  static struct option opts[] = {
    {"x", 1, 0, 0},
    {"X", 1, 0, 1},
    {"y", 1, 0, 2},
    {"Y", 1, 0, 3},
    {"z", 1, 0, 4},
    {"Z", 1, 0 ,5}
  };

  int option_index = 0;
  while((c = getopt_long(argc, argv, "f:o:b:Lp", opts, &option_index)) != -1) {
      switch(c) {
	case 'f':
	  strcpy(inFile, optarg);
	  break;
	case 'o':
	  strcpy(outfile, optarg);
	  break;
	case 'b':
	  strcpy(blocksFileName, optarg);
	  break;
	case 'L':
	  leafsonly = true;
	  break;
	case 'p':
	  probe = true;
	  break;
	case 0:
	  x0 = atof(optarg); 
	  break;
	case 1:
	  x1 = atof(optarg); 
	  break;
	case 2:
	  y0 = atof(optarg); 
	  break;
	case 3:
	  y1 = atof(optarg); 
	  break;
	case 4:
	  z0 = atof(optarg); 
	  break;
	case 5:
	  z1 = atof(optarg); 
	  break;
	}
    }


  //Error check
  if(inFile[0] == '\0' || outfile[0] == '\0' && !probe)
    help(argv[0]);


  FlashHDFFile fdf(inFile);
  nb = fdf.GetNumberOfBlocks();
  fdf.GetCellDimensions(cellDims);
  nxb = cellDims[0];
  nyb = cellDims[1];
  nzb = cellDims[2];

  // probe
  if(probe) {
    fdf.GetCoordinateRangeEntireDataset(bounds);
    for(i=0; i<nb; i++) {
      fdf.Get3dBlockSize(i,bsize);
      if(bsize[0] < min) {
	min = bsize[0];
      }
      if(bsize[0] > max) {
	max = bsize[0];
      }
    } //end for i
    fprintf(stderr, "X: %f - %f Y: %f - %f Z: %f - %f\n",bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
    fprintf(stderr, "Block Size min: %f Block Size max: %f\n",min,max);
    exit(0);
  } // probe

			
  fprintf(stderr, "Coordinate Ranges: %f %f %f %f %f\n",x0,x1,y0,y1,z0,z1);

  // form a list of blocks in the desired min, max range
  for(i=0; i<nb; i++) {
    if(leafsonly) {
      if(fdf.GetNodeType(i) == 1) {
	fdf.Get3dCoordinate(i,coords);
	if( (coords[0] >= x0 && coords[0] <= x1) && (coords[1] >= y0 && coords[1] <= y1) && (coords[2] >= z0 && coords[2] <= z1) ) { 
	  blk_list.push_back(i);
	}
      }
    } else {
      fdf.Get3dCoordinate(i,coords);
      if( (coords[0] >= x0 && coords[0] <= x1) && (coords[1] >= y0 && coords[1] <= y1) && (coords[2] >= z0 && coords[2] <= z1) ) {
	blk_list.push_back(i);
      }
    } //endif leafsonly
  } //endfor i

  // write blocks list in text file					
  if(blocksFileName[0] != '\0') {
    ofstream ofile;
    ofile.open(blocksFileName, ios::out);
    for(std::vector<int>::iterator itV = blk_list.begin(); itV != blk_list.end(); itV++) {
      ofile << *itV << endl;
    }
    ofile.close();
  }

  //now read each block of data and write 
  fdata1 = new float[nxb*nyb*nzb];
  fdata2 = new float[nxb*nyb*nzb];
  fdata3 = new float[nxb*nyb*nzb];

  FILE *Outfile;
  Outfile = fopen(outfile, "wb");
  nb_out = blk_list.size();
  fwrite(&nb_out, sizeof(int), 1, Outfile);
  fwrite(cellDims, sizeof(int), 3, Outfile);
  for(std::vector<int>::iterator itV = blk_list.begin(); itV != blk_list.end(); itV++) {
    fdf.Get3dCoordinate(*itV,coords);
    fdf.Get3dBlockSize(*itV,bsize);
    fdf.Get3dBoundingBox(*itV,bnd_box);
    fdf.GetScalarVariable(velx,*itV,fdata1);
    fdf.GetScalarVariable(vely,*itV,fdata2);
    fdf.GetScalarVariable(velz,*itV,fdata3);
    //check to see endianess of system, always write little endian
    //this function will always make sure arrays are little endian

    block_no = *itV;
    fwrite(&block_no, sizeof(int), 1, Outfile);
    fwrite(coords, sizeof(float), COORD_SIZE, Outfile);
    fwrite(bsize, sizeof(float), COORD_SIZE, Outfile);
    fwrite(bnd_box, sizeof(float), BND_SIZE, Outfile);
    fwrite(fdata1, sizeof(float), nxb*nyb*nzb, Outfile);
    fwrite(fdata2, sizeof(float), nxb*nyb*nzb, Outfile);
    fwrite(fdata3, sizeof(float), nxb*nyb*nzb, Outfile);
  }

  fprintf(stderr, "Finished writing %d blocks of data to file %s from filr %s.\n",blk_list.size()*3,outfile,inFile);
  fclose(Outfile);
  delete [] fdata1;
  delete [] fdata2;
  delete [] fdata3;
  fdf.Close();
  exit(0);
} 
