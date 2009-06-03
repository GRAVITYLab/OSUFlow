
//////////////////////////////////////////////////////////////////////
//
//  Sequential Functions for streamline and pathline tracings in 
//  multiple subdomains 
//
//  Han-Wei Shen (hwshen@cse.ohio-state.edu)
//  The Ohio State University 
//  May 30, 2009  at Argonne National Laboratory 
//
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h> 

#include "OSUFlow.h"
#include "Lattice.h"
#include "Lattice4D.h"
#include "calc_subvolume.h"

#include <list>
#include <iterator>

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
//   A sequential function for pathline tracing in all partitions
//
/////////////////////////////////////////////////////////////////////

list<vtListTimeSeedTrace*>* 
ComputePathlines(OSUFlow** osuflow_list, 
		  Lattice4D* lat, 
		  volume_bounds_type *vb_list, 
		  int nparts, 
		  int num_seeds_per_part) 
{

  float from[3], to[3]; 
  VECTOR4 **osuflow_seeds = new VECTOR4*[nparts]; 
  int *osuflow_num_seeds = new int[nparts]; 
  list<vtListTimeSeedTrace*> *sl_list = new list<vtListTimeSeedTrace*>[nparts]; 
  for (int i=0; i<nparts; i++) {

    VECTOR3 * seeds; 
    int num; 

    from[0] = vb_list[i].xmin;  to[0] = vb_list[i].xmax;  
    from[1] = vb_list[i].ymin;  to[1] = vb_list[i].ymax;
    from[2] = vb_list[i].zmin;  to[2] = vb_list[i].zmax; 

    // generate 3D seeds, add time later 
    osuflow_list[i]->SetRandomSeedPoints(from, to, num_seeds_per_part); 
    seeds = osuflow_list[i]->GetSeeds(num); 
    osuflow_num_seeds[i] = num; 

    osuflow_seeds[i] = new VECTOR4[num];    // now copy and augment to 4D seeds 
    for (int j=0; j<num; j++)  {
      osuflow_seeds[i][j][0] = seeds[j][0]; 
      osuflow_seeds[i][j][1] = seeds[j][1]; 
      osuflow_seeds[i][j][2] = seeds[j][2]; 
      // osuflow_seeds[i][j][3] = vb_list[i].tmin; 
       osuflow_seeds[i][j][3] = 0; 
    }
    sl_list[i].clear();   // clear the trace 
  }

  // Perform pathline tracing in all subdomains 
  bool has_seeds = true;      
  int num_seeds_left = nparts*num_seeds_per_part; 

  // loop until all particles stop 
  while(has_seeds == true && num_seeds_left >50) {  

    lat->ResetSeedLists();    // clear up the lattice seed lists

    for (int i=0; i<nparts; i++) {

      if (osuflow_num_seeds[i]==0) {  // domain i is done. 
	printf("skip domain %d \n", i); 
	continue; 
      }
      list<vtListTimeSeedTrace*> list; 
      osuflow_list[i]->SetIntegrationParams(1, 5); 
      osuflow_list[i]->GenPathLines(osuflow_seeds[i],list, FORWARD, 
      				    osuflow_num_seeds[i], 5000); 
      printf(" Generate %d pathlines in partition %d.\n", 
	     list.size(), i); 

      // Insert the pathline traces points
      std::list<vtListTimeSeedTrace*>::iterator pIter; 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
        vtListTimeSeedTrace *trace = *pIter; 
	sl_list[i].push_back(trace); 
      }

      // Redistributing the pathline boundary points
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
	vtListTimeSeedTrace *trace = *pIter; 
	if (trace->size() ==0) continue; 
	std::list<VECTOR4*>::iterator pnIter; 
	pnIter = trace->end(); 
	pnIter--; 
	VECTOR4 p = **pnIter; 
	//check where p should go 
	int ei, ej, ek, et; 
	int neighbor = lat->GetNeighbor(i, p[0], p[1], p[2], p[3], ei, ej, ek, et); 
	if (neighbor!=-1) lat->InsertSeed(ei, ej, ek, et, p); 
	printf(" insert a seed %f %f %f %f to rank %d \n",
	       p[0], p[1], p[2], p[3], neighbor); 
      }
    }
    // Update the seed lists for the next run
    has_seeds = false;  
    num_seeds_left = 0; 
    for (int i=0; i<nparts; i++) {

      osuflow_num_seeds[i] = lat->seedlists[i].size(); 
      num_seeds_left += osuflow_num_seeds[i]; 
      printf("seedlists[%d].size() = %d\n", i, osuflow_num_seeds[i]); 
      if (osuflow_num_seeds[i]!=0) has_seeds = true; 
      else continue; 
      //  don't you want to free the memory? why not?
      //  if (osuflow_seeds[i]!=NULL) delete [] osuflow_seeds[i]; 
      osuflow_seeds[i] = new VECTOR4[osuflow_num_seeds[i]]; 
      std::list<VECTOR4>::iterator seedIter; 
      seedIter = lat->seedlists[i].begin(); 
      int cnt = 0; 
      for (; seedIter!=lat->seedlists[i].end(); seedIter++){
	VECTOR4 p = *seedIter; 
	osuflow_seeds[i][cnt++] = p; 
      }
    }
  }
  return(sl_list); 
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
//   A sequential function for streamline tracing in all partitions
//
/////////////////////////////////////////////////////////////////////

list<vtListSeedTrace*>*
ComputeStreamlines(OSUFlow** osuflow_list, 
		   Lattice* lat, 
		   volume_bounds_type *vb_list, 
		   int nparts, 
		   int num_seeds_per_part) 
{
  int num;   
  float from[3], to[3]; 
  VECTOR3 **osuflow_seeds = new VECTOR3*[nparts]; 
  int *osuflow_num_seeds = new int[nparts]; 
  list<vtListSeedTrace*> *sl_list = new list<vtListSeedTrace*>[nparts]; 

  // initialize the seed sets
  for (int i=0; i<nparts; i++) {

    from[0] = vb_list[i].xmin;  
    from[1] = vb_list[i].ymin;     
    from[2] = vb_list[i].zmin; 

    to[0] = vb_list[i].xmax;  
    to[1] = vb_list[i].ymax;     
    to[2] = vb_list[i].zmax; 

    // set range for seed locations
    osuflow_list[i]->SetRandomSeedPoints(from, to, num_seeds_per_part); 
    osuflow_seeds[i] = osuflow_list[i]->GetSeeds(num); 
    osuflow_num_seeds[i] = num; 
    sl_list[i].clear(); 
  }

  // Now perform particle tracing in all subdomains
  bool has_seeds = true;   
  int num_seeds_left = num_seeds_per_part * nparts; 

  // loop until all particles stop 
  while(has_seeds == true && num_seeds_left >10) {  

    for (int i=0; i<nparts; i++) {

      if (osuflow_num_seeds[i]==0) {  // this partition is done. 
	printf("skip partition %d \n", i); 
	continue; 
      }
      list<vtListSeedTrace*> list; 
      osuflow_list[i]->SetIntegrationParams(1, 5); 
      osuflow_list[i]->GenStreamLines(osuflow_seeds[i], FORWARD_DIR, 
				      osuflow_num_seeds[i], 50, list); 
      printf(" Generate %d streamlines in partition %d. \n", 
	     list.size(), i); 

      //Insert the particle traces 
      std::list<vtListSeedTrace*>::iterator pIter; 
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
        vtListSeedTrace *trace = *pIter; 
	sl_list[i].push_back(trace); 
      }
      // Redistribute the streamline boundary points
      pIter = list.begin(); 
      for (; pIter!=list.end(); pIter++) {
	vtListSeedTrace *trace = *pIter; 
	if (trace->size() ==0) continue; 
	std::list<VECTOR3*>::iterator pnIter; 
	pnIter = trace->end(); // The last point is the boundary point 
	pnIter--; 
	VECTOR3 p = **pnIter; 
	//check which partition p is in 
	int neighbor = lat->CheckNeighbor(i, p[0], p[1], p[2]); 
	int si, sj, sk, ei, ej, ek; 
	lat->GetIndices(i, si, sj, sk); //where am I in the lattice?
	if (neighbor ==0) {ei=si-1; ej = sj; ek = sk;}
	else if (neighbor ==1) {ei=si+1; ej = sj; ek = sk;}
	else if (neighbor ==2) {ei=si; ej = sj-1; ek = sk;}
	else if (neighbor ==3) {ei=si; ej = sj+1; ek = sk;}
	else if (neighbor ==4) {ei=si; ej = sj; ek = sk-1;}
	else if (neighbor ==5) {ei=si; ej = sj; ek = sk+1;}
	if (neighbor!=-1) lat->InsertSeed(ei, ej, ek, p); 
      }
    }
    // now update the seed lists for the next run
    has_seeds = false;  
    num_seeds_left = 0; 
    for (int i=0; i<nparts; i++) {
      osuflow_num_seeds[i] = lat->seedlists[i].size(); 
      num_seeds_left += osuflow_num_seeds[i]; 
      if (osuflow_num_seeds[i]!=0) has_seeds = true; 
      else continue; 
      // ?? you should delete the old seed list, why not? 
      // if (osuflow_seeds[i]!=0) delete [] osuflow_seeds[i]; 
      osuflow_seeds[i] = new VECTOR3[osuflow_num_seeds[i]]; 
      std::list<VECTOR3>::iterator seedIter; 
      seedIter = lat->seedlists[i].begin(); 
      int cnt = 0; 
      for (; seedIter!=lat->seedlists[i].end(); seedIter++){
	VECTOR3 p = *seedIter; 
	osuflow_seeds[i][cnt++] = p; 
      }
      lat->ClearSeedList(i); 
    }
  }
  return(sl_list); 
}
