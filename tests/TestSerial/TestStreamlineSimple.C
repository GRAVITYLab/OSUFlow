

#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

int	// ADD-BY-LEETEN 12/20/2011
main(int argc, char**argv) {

  VECTOR3 minLen, maxLen; 

  printf("Testing streamline geneation...\n"); 

  OSUFlow *osuflow = new OSUFlow(); 
  osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 
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

  list<vtListSeedTrace*> list; 
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreamLines(list , FORWARD_DIR, 50, 0); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)list.size()); 

  std::list<vtListSeedTrace*>::iterator pIter; 


  pIter = list.begin(); 
  for (; pIter!=list.end(); pIter++) {
    vtListSeedTrace *trace = *pIter; 
    std::list<VECTOR3*>::iterator pnIter; 
    pnIter = trace->begin(); 
    for (; pnIter!=trace->end(); pnIter++) {

      VECTOR3 p = **pnIter; 
      printf("%f %f %f, ", p[0], p[1], p[2]); 


    }
    printf("\n");
  }
  return 0;	// ADD-BY-LEETEN 12/20/2011
}
