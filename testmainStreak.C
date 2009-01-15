

#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

main(int argc, void *argv[]) {

  VECTOR3 minLen, maxLen; 

  printf("hello! entering testmainT (streaklines for time-varying data)...\n"); 

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

  vtStreakTraces list; 
  osuflow->SetIntegrationParams(1, 5); 

  float ctime = 0.0; 
  float max_time = 49.0; 

  osuflow->GenStreakLines(list , FORWARD, ctime); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", list.size()); 

  vtStreakTracesIter pIter; 
  pIter = list.begin(); 
  for (; pIter!=list.end(); pIter++) {
    vtListStreakParticle *trace = *pIter; 
    vtStreakParticleIter  pnIter; 
    pnIter = trace->begin(); 
    for (; pnIter!=trace->end(); pnIter++) {
      vtStreakParticle p = **pnIter; 
      float x = p.itsPoint.phyCoord[0]; 
      float y = p.itsPoint.phyCoord[1]; 
      float z = p.itsPoint.phyCoord[2]; 
      printf(" %f %f %f ", x, y, z); 
    }
  }


  /*
  while (ctime < max_time)  {
    printf(" --------------- Time =  %f ----------------\n", ctime); 
    osuflow->GenStreakLines(list , FORWARD, ctime); 
    printf(" done integrations\n"); 
    printf("list size = %d\n", list.size()); 

    vtStreakTracesIter pIter; 
    pIter = list.begin(); 
    for (; pIter!=list.end(); pIter++) {
      vtListStreakParticle *trace = *pIter; 
      vtStreakParticleIter  pnIter; 
      pnIter = trace->begin(); 
      for (; pnIter!=trace->end(); pnIter++) {
	vtStreakParticle p = **pnIter; 
	float x = p.itsPoint.phyCoord[0]; 
	float y = p.itsPoint.phyCoord[1]; 
	float z = p.itsPoint.phyCoord[2]; 
	//	printf(" %f %f %f ", x, y, z); 
      }
    }

    ctime += 1.0; 
  }
  */
}
