//////////////////////////////////////////////////////
//
//    A test program to generate particles in 
//    a decomposed domain
//
//    Han-Wei Shen, 11/15/2008 at Argonne National Laboratory
//
#include <stdio.h>
#include <stdlib.h> 
#include "OSUFlow.h"

#include <list>
#include <iterator>

#include "calc_subvolume.h"
#include "Lattice.h" 

#define MAX_ITERATIONS 2

int	// ADD-BY-LEETEN 01/06/2012
main(int argc, void *argv[]) {


  int npart = 2;   // number of subdomains we will create 
  int nproc = 2;   // number of subdomains we will create 
  OSUFlow **osuflow_list = new OSUFlow*[npart];  
	#if	0	// MOD-BY-LEETEN 01/06/2012-FROM:
	  VECTOR3 *osuflow_seeds[npart]; 
	  int osuflow_num_seeds[npart]; 
	#else	// MOD-BY-LEETEN 01/06/2012-TO:
  VECTOR3* *osuflow_seeds = new VECTOR3 *[npart]; 
  int* osuflow_num_seeds = new int[npart]; 
	#endif	// MOD-BY-LEETEN 01/06/2012-END
  VECTOR3 minLen, maxLen; 
  volume_bounds_type *vb_list; 
  int num_seeds_left=0; 
  int **plist; 
  int *num_partitions; 

  OSUFlow *osuflow = new OSUFlow(); 
  printf("read file %s\n", argv[1]); 
  //loading the whole dataset just to get dims. 
  //obviously not very smart. need to change. 
  osuflow->LoadData((const char*)argv[1], true); 
                                                       
  osuflow->Boundary(minLen, maxLen);    // query the dims
  printf(" volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
	 minLen[0], maxLen[0], minLen[1], maxLen[1], 
	 minLen[2], maxLen[2]); 

  // -------------------------------------
  // now subdivide the entire domain into nproc subdomains

  int lattice_xdim, lattice_ydim, lattice_zdim; 
  // partition the domain and create a lattice 
  Lattice* lat = new Lattice(maxLen[0]-minLen[0]+1, maxLen[1]-minLen[1]+1, 
			     maxLen[2]-minLen[2]+1, 0, npart);
  vb_list = lat->GetBoundsList(); 
  // lat->GetLatticeDims(lattice_xdim, lattice_ydim, lattice_zdim); 
  lat->InitSeedLists(); 
  lat->RoundRobin_proc(nproc); 

  plist = new int*[nproc]; 
  num_partitions = new int[nproc]; 
  for (int i=0; i<nproc; i++)
    lat->GetPartitions(i, &(plist[i]), num_partitions[i]);
  


  // -------------------------------------
  // now create a list of flow field for the subdomains 
  for (int i=0; i<npart; i++) {
    osuflow_list[i] = new OSUFlow(); 
    printf("Domain(%d):  %d %d %d : %d %d %d\n", i, vb_list[i].xmin,  
	   vb_list[i].ymin,  vb_list[i].zmin, vb_list[i].xmax,  
	   vb_list[i].ymax,  vb_list[i].zmax); 

    // load subdomain data into OSUFlow
    VECTOR3 minB, maxB; 
    minB[0] = vb_list[i].xmin;  
    minB[1] = vb_list[i].ymin;     
    minB[2] = vb_list[i].zmin; 
    maxB[0] = vb_list[i].xmax;  
    maxB[1] = vb_list[i].ymax;     
    maxB[2] = vb_list[i].zmax; 
    osuflow_list[i]->LoadData((const char*)argv[1], true, minB, maxB); 

    //generating seeds
    float from[3], to[3]; 
    from[0] = minB[0]; to[0] = maxB[0];
    from[1] = minB[1]; to[1] = maxB[1];
    from[2] = minB[2]; to[2] = maxB[2]; 
    osuflow_list[i]->SetRandomSeedPoints(from, to, 10); 
    int num; 
    osuflow_seeds[i] = osuflow_list[i]->GetSeeds(num); 
    osuflow_num_seeds[i] = num; 
  }
  //-------------------------------------------------------
  // Now begin to perform particle tracing in all subdomains
  bool has_seeds = true;      // initially we always have seeds

  int counter = 0; 
  while(has_seeds == true && counter <MAX_ITERATIONS) {  // loop until all particles stop 
    //     while(has_seeds == true) {  // loop until all particles stop 
      lat->ResetSeedLists();    // clear up the lattice seed lists

    for (int proc=0; proc<nproc;  proc++) 
     for (int np=0; np<num_partitions[proc]; np++) {
        int i = plist[proc][np]; 
	if (osuflow_num_seeds[i]==0) {  // npart is already done. 
	  printf("skip domain %d \n", i); 
	  continue; 
	}
	list<vtListSeedTrace*> list; 
	osuflow_list[i]->SetIntegrationParams(1, 5); 
	osuflow_list[i]->GenStreamLines(osuflow_seeds[i], FORWARD_DIR, 
					osuflow_num_seeds[i], 50, list); 

	for (int s=0; s<osuflow_num_seeds[i]; s++) 
	  fprintf(stderr, "Rank %d Seed %d:  %f %f %f \n", i, s, osuflow_seeds[i][s][0], 
		  osuflow_seeds[i][s][1], osuflow_seeds[i][s][2]); 
	printf("domain %d done integrations", i); 
	//	printf(" %d streamlines. \n", list.size()); 
	printf(" %d streamlines. \n", osuflow_num_seeds[i]); 

	std::list<vtListSeedTrace*>::iterator pIter; 
	//------------------------------------------------
	//looping through the trace points. just checking. 
	pIter = list.begin(); 
	for (; pIter!=list.end(); pIter++) {
	  vtListSeedTrace *trace = *pIter; 
	  std::list<VECTOR3*>::iterator pnIter; 
	  pnIter = trace->begin(); 
	  for (; pnIter!=trace->end(); pnIter++) {
	    VECTOR3 p = **pnIter; 
	    // printf(" %f %f %f ", p[0], p[1], p[2]); 
	  }
	}
	//---------------
	// now redistributing the boundary streamline points to its neighbors. 
	pIter = list.begin(); 
	for (; pIter!=list.end(); pIter++) {
	  vtListSeedTrace *trace = *pIter; 
	  if (trace->size() ==0) continue; 
	  std::list<VECTOR3*>::iterator pnIter; 
	  pnIter = trace->end(); 
	  pnIter--; 
	  VECTOR3 p = **pnIter; 
	  //check p is in which neighbor's domain 
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
// 	  printf(" insert a seed to rank %d \n", lat->GetRank(ei,ej, ek)); 
//	  fprintf(stderr,"rank %d inserting point: %.3f\t%.3f\t%.3f\n",i,p[0],p[1],p[2]);

	}
      }

      //-------------

      // now create the seed arrays for the next run
      has_seeds = false;  
      for (int i=0; i<npart; i++) {
	// if (osuflow_seeds[i]!=0) delete [] osuflow_seeds[i]; 
	osuflow_num_seeds[i] = lat->seedlists[i].size(); 
	printf("seedlists[%d].size() = %d\n", i, osuflow_num_seeds[i]); 
	if (osuflow_num_seeds[i]!=0) has_seeds = true; 
	else continue; 
	osuflow_seeds[i] = new VECTOR3[osuflow_num_seeds[i]]; 
	std::list<VECTOR3>::iterator seedIter; 
	seedIter = lat->seedlists[i].begin(); 
	int cnt = 0; 
	for (; seedIter!=lat->seedlists[i].end(); seedIter++){
	  VECTOR3 p = *seedIter; 
	  osuflow_seeds[i][cnt++] = p; 
	}
      }
      counter++; 
      printf(" *** Iteration %d done. \n", counter); 
    }

	// ADD-BY-LEETEN 01/06/2012-BEGIN
  delete [] osuflow_seeds; 
  delete [] osuflow_num_seeds; 
	// ADD-BY-LEETEN 01/06/2012-END


	return 0;	// ADD-BY-LEETEN 01/06/2012
}
