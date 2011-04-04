

#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

main(int argc, char** argv) {

  VECTOR3 minLen, maxLen; 

  printf("hello! entering testmainT (time-varying data test)...\n"); 

  OSUFlow *osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 
  osuflow->LoadData((const char*)argv[1], false); //flase: a time-varying flow field 
  osuflow->Boundary(minLen, maxLen); 
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                minLen[0], maxLen[0], minLen[1], maxLen[1], 
                                minLen[2], maxLen[2]); 

  float from[3], to[3]; 
  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 
  osuflow->SetRandomSeedPoints(from, to, 100); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  list<vtListTimeSeedTrace*> list; 
  osuflow->SetIntegrationParams(1, 5); 

  
  osuflow->GenPathLines(list , FORWARD, 50); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)list.size()); 


  std::list<vtListTimeSeedTrace*>::iterator pIter; 


  pIter = list.begin(); 
  for (; pIter!=list.end(); pIter++) {
    vtListTimeSeedTrace *trace = *pIter; 
    std::list<VECTOR4*>::iterator pnIter; 
    pnIter = trace->begin(); 
    for (; pnIter!=trace->end(); pnIter++) {

      VECTOR4 p = **pnIter; 
      printf(" %f %f %f %f", p[0], p[1], p[2], p[3]); 


    }

  }

}
