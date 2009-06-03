
#include "Lattice.h"
#include "Lattice4D.h"

list<vtListTimeSeedTrace*>* 
ComputePathlines(OSUFlow** osuflow_list, 
		  Lattice4D* lat, 
		  volume_bounds_type *vb_list, 
		  int nparts, 
		 int num_seeds_per_part) ; 

list<vtListSeedTrace*>*
ComputeStreamlines(OSUFlow** osuflow_list, 
		   Lattice* lat, 
		   volume_bounds_type *vb_list, 
		   int nparts, 
		   int num_seeds_per_part); 
