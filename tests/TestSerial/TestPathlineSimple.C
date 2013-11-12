

#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

int main(int argc, char** argv) {

  VECTOR3 minLen, maxLen; 

  printf("Testing pathline generation\n"); 

  OSUFlow *osuflow = new OSUFlow(); 
  osuflow->LoadData((const char*)argv[1], false); //flase: a time-varying flow field 
  osuflow->Boundary(minLen, maxLen); 
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
                                minLen[0], maxLen[0], minLen[1], maxLen[1], 
                                minLen[2], maxLen[2]); 

  float from[3], to[3]; 
  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 
  osuflow->SetRandomSeedPoints(from, to, 200); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], 
	   seeds[i][1], seeds[i][2]); 

  list<vtListTimeSeedTrace*> list; 
  osuflow->SetIntegrationParams(1, 5);
  osuflow->ScaleField(20.0);
  osuflow->SetMaxError(0.0001);
  osuflow->SetIntegrationParams(1, 0.01, 5);

  
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

#if 0
  // write output to a file
  float zero = 0.0f;
  FILE* f = fopen("field_lines.out", "wb");
  fwrite(from, sizeof(float), 3, f);
  fwrite(&zero, sizeof(float), 1, f);
  fwrite(to, sizeof(float), 3, f);
  fwrite(&zero, sizeof(float), 1, f);

  for (pIter=list.begin(); pIter!=list.end(); pIter++)
  {
    vtListTimeSeedTrace *trace = *pIter; 
    int a = trace->size();
    fwrite(&a, sizeof(int), 1, f);
  }
  int neg = -1;
  fwrite(&neg, sizeof(int), 1, f);

  pIter = list.begin(); 
  for (; pIter!=list.end(); pIter++) {
    vtListTimeSeedTrace *trace = *pIter; 
    std::list<VECTOR4*>::iterator pnIter; 
    pnIter = trace->begin(); 
    for (; pnIter!=trace->end(); pnIter++) {

      VECTOR4 p = **pnIter; 
      //printf(" %f %f %f %f", p[0], p[1], p[2], p[3]); 
      fwrite(&p[0], sizeof(float), 1, f);
      fwrite(&p[1], sizeof(float), 1, f);
      fwrite(&p[2], sizeof(float), 1, f);
      fwrite(&p[3], sizeof(float), 1, f);
    }
  }

  fclose(f);
#endif

  return 0;
}
