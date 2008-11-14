

#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

#include "calc_subvolume.h"

main(int argc, void *argv[]) {

  VECTOR3 minLen, maxLen; 
  VECTOR3 minB, maxB; 
  volume_bounds_type *vb_list; 
  int nproc = 8; 

  printf("hello! entering testmain...\n"); 

  OSUFlow *osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 
  osuflow->LoadData((const char*)argv[1], true); //loading the whole dataset just to get dims. 
                                                 //obviously not very smart. need to change. 
                                                       
  osuflow->Boundary(minLen, maxLen);    // query the dims
  minB[0] = minLen[0]/2; minB[1] = minLen[1]/2;  minB[2] = minLen[2]/2;
  maxB[0] = maxLen[0]/2; maxB[1] = maxLen[1]/2;  maxB[2] = maxLen[2]/2;
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                minLen[0], maxLen[0], minLen[1], maxLen[1], 
                                minLen[2], maxLen[2]); 

  // subdivide the entire domain into nproc subdomains

  int* lattice;
  int lattice_xdim, lattice_ydim, lattice_zdim; 
  vb_list = calc_subvolume(maxLen[0]-minLen[0], maxLen[1]-minLen[1], maxLen[2]-minLen[2], 2, nproc, &lattice, lattice_xdim, lattice_ydim, 
			   lattice_zdim); 

  // create a list of flow field for the subdomains 
  OSUFlow **osuflow_list = new OSUFlow*[nproc];  

  for (int i=0; i<nproc; i++) {
    osuflow_list[i] = new OSUFlow(); 
    printf("PE %d:  %d %d %d : %d %d %d\n", i, vb_list[i].xmin,  vb_list[i].ymin,  vb_list[i].zmin, 
	   vb_list[i].xmax,  vb_list[i].ymax,  vb_list[i].zmax); 
    VECTOR3 minB, maxB; 
    minB[0] = vb_list[i].xmin;  minB[1] = vb_list[i].ymin;     minB[2] = vb_list[i].zmin; 
    maxB[0] = vb_list[i].xmax;  maxB[1] = vb_list[i].ymax;     maxB[2] = vb_list[i].zmax; 
    osuflow_list[i]->LoadData((const char*)argv[1], true, minB, maxB); 

    float from[3], to[3]; 
    from[0] = minB[0];   from[1] = minB[1];   from[2] = minB[2]; 
    to[0] = maxB[0];   to[1] = maxB[1];   to[2] = maxB[2]; 
    osuflow_list[i]->SetRandomSeedPoints(from, to, 20); 

    int nSeeds; 
    VECTOR3* seeds = osuflow_list[i]->GetSeeds(nSeeds); 
    for (int j=0; j<nSeeds; j++) 
      printf(" seed no. %d : [%f %f %f]\n", j, seeds[j][0], 
	     seeds[j][1], seeds[j][2]); 

    list<vtListSeedTrace*> list; 
    osuflow_list[i]->SetIntegrationParams(1, 5); 
    osuflow_list[i]->GenStreamLines(list , FORWARD_DIR, 50, 0); 
    printf(" domain %d done integrations\n", i); 
    printf("list size = %d\n", list.size()); 

    std::list<vtListSeedTrace*>::iterator pIter; 
    pIter = list.begin(); 
    for (; pIter!=list.end(); pIter++) {
      vtListSeedTrace *trace = *pIter; 
      std::list<VECTOR3*>::iterator pnIter; 
      pnIter = trace->begin(); 
      for (; pnIter!=trace->end(); pnIter++) {

      VECTOR3 p = **pnIter; 
      //      printf(" %f %f %f ", p[0], p[1], p[2]); 

      }
    }
  }
}
