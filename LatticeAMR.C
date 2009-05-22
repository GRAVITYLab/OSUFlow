
#include "LatticeAMR.h"

/////////////////////////////////////////////////////////
//
// Lattice AMR
//
// constructs and initializes the time-varying regular structured 
// AMR lattice 
//
LatticeAMR::LatticeAMR(float xlen, float ylen, float zlen, float tlen, 
		       int total_level) {

  num_levels = total_level; // total AMR level. 

  xdim = xlen;  // the physical dimension of the whole domain 
  ydim = ylen; 
  zdim = zlen; 
  tdim = tlen;  // number of time steps

  // the bounds of each level 
  xmin = new float[num_levels];   xmax = new float[num_levels]; 
  ymin = new float[num_levels];   ymax = new float[num_levels]; 
  zmin = new float[num_levels];   zmax = new float[num_levels]; 
  tmin = new float[num_levels];   tmax = new float[num_levels]; 

  // the xyz lengths of a block in each level 
  xlength = new float[num_levels]; 
  ylength = new float[num_levels]; 
  zlength = new float[num_levels]; 
  tlength = new float[num_levels]; 

  // the data resolution of blocks in all level (they should be the same) 
  xres = new int[num_levels]; 
  yres = new int[num_levels]; 
  zres = new int[num_levels]; 
  tres = new int[num_levels]; 

  // the dimensions of the lattice in each level // for partitioning purpose 
  idim = new int[num_levels]; 
  jdim = new int[num_levels]; 
  kdim = new int[num_levels]; 
  ldim = new int[num_levels]; 

  // data structures for book-keeping 
  has_data = new bool*[num_levels];     // whether the element has data or not 
  has_data_from_merger = new int*[num_levels]; // used for block merger 
  finest_level = new int*[num_levels];   // what is the finest level of data 
                                         // in this spatial region
  data_ptr = new float**[num_levels]; 
  index_to_rank = new int*[total_level]; // mapping from index to rank 
  nblocks = new int[total_level];        // how many blocks in each region 
                                         // regardless of empty or not 
  npart = 0;                             // number of non-empty blocks total
}

/////////////////////////////////////////////////////////////////
//
//  Destructor
// 
LatticeAMR::~LatticeAMR()
{
  delete []xmin;   delete []ymin;   delete []zmin; delete []tmin;  
  delete []xmax;   delete []ymax;   delete []zmax; delete []tmax; 
  delete []idim;   delete []jdim;   delete []kdim; delete []ldim; 
  delete []xlength;  delete []ylength;  delete []zlength; delete []tlength; 
  delete [] rank_to_index; 
  delete [] nblocks; 

  for (int i=0; i<num_levels; i++) {
    delete [] has_data[i]; 
    delete [] finest_level[i]; 
    delete [] index_to_rank[i]; 
  }
  delete [] has_data; 
  delete [] finest_level; 
  delete [] index_to_rank; 
  delete [] rank_to_index; 
  delete [] parts; 

  if (vb_list!=NULL) delete [] vb_list; 
  if (seedlists!=NULL) delete [] seedlists; 
}

//////////////////////////////////////////////////////////////////////////
//
// Create a new level of lattice. The level has to be smaller than the 
// largest level set in the constructor 
//
// level: the level id [0..max_level-1]
// {x,y,z}_size: the lengths of the block in this level 
// {x,y,z,t}_res: the data resolution in each dimension 
// {x,y,z,t}_{min,max}: the physical bounds of the level 
// 
bool LatticeAMR::CreateLevel(int level, float x_size, float y_size, 
			     float z_size, 
			     int x_res, int y_res, int z_res, 
			     float x_min, float x_max, 
			     float y_min, float y_max, float z_min, 
			     float z_max, float t_min, float t_max) 
{

  if (level <0 || level >=num_levels) return false; 

  xlength[level]=x_size; ylength[level]=y_size; zlength[level]=z_size;
  tlength[level]= 1.0;  // one time step at a time for an AMR block 
  
  // the data resolution 
  xres[level] = x_res; yres[level] = y_res; zres[level] = z_res; 
  tres[level] = 1;      // data always come one time step at a time 
 
  // the physical bounds of this level 
  xmin[level] = x_min; xmax[level] = x_max; 
  ymin[level] = y_min; ymax[level] = y_max; 
  zmin[level] = z_min; zmax[level] = z_max;
  tmin[level] = t_min; tmax[level] = t_max;  // the total time range 
  
  //i,j,k,l dim are the number of blocks in each dimension if the whole 
  //space-time domain is filled with blocks from this level 
  // this is the lattice dimensions for this level 
  // 
  idim[level] = (int) (x_max-x_min)/x_size; 
  jdim[level] = (int) (y_max-y_min)/y_size; 
  kdim[level] = (int) (z_max-z_min)/z_size; 
  ldim[level] = (int) (t_max-t_min)/tlength[level]; 

  printf("**** level %d dims %d %d %d %d \n", level, idim[level], jdim[level], 
  	 kdim[level], ldim[level]); 

  int size = idim[level]*jdim[level]*kdim[level]*ldim[level]; // total number of lattice elements 
  nblocks[level] = size;  // number of blocks at this level, remember not all blocks have data
  printf(" nblock[%d] = %d \n", level, size); 

  has_data[level] = new bool[size]; 
  has_data_from_merger[level] = new int[size]; 
  finest_level[level] = new int[size]; 
  data_ptr[level] = new float*[size]; 
  index_to_rank[level] = new int[size];

  int idx = 0; 
  for (int t=0; t<ldim[level]; t++) 
    for (int k=0; k<kdim[level]; k++) 
      for (int j=0; j<jdim[level]; j++) 
	for (int i=0; i<idim[level]; i++) {
	  has_data[level][idx] = false;    // initial value: no data
	  has_data_from_merger[level][idx] = -1; 
	  finest_level[level][idx] = -1;   // initial value: unknown
	  data_ptr[level][idx] = NULL; 
	  index_to_rank[level][idx] = -1;  // no mapping available 
	  idx++; 
	}
  return true; 
}

//////////////////////////////////////////////////////////////////////////
//
//  Get the index of the block at the given level that contains (x,y,z,t)
//
int LatticeAMR::GetIndexinLevel(int level, float x, float y, float z, float t)
{
  if (level <0 || level >=num_levels) return -1; 
  if (x<xmin[level] || x>xmax[level] || y<ymin[level] || y>ymax[level] || 
      z<zmin[level] || z>zmax[level] || t<tmin[level] || t>ymax[level])
    return -1; 
  
  int i = (int)((x-xmin[level])/xlength[level]); 
  int j = (int)((y-ymin[level])/ylength[level]); 
  int k = (int)((z-zmin[level])/zlength[level]); 
  int l = (int)((t-tmin[level])/tlength[level]); 
  int idx = l*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 
  return idx; 
}

///////////////////////////////////////////////////////////
//
// check in to indicate that the lattice element containing 
// (x,y,z,t) has data at this level 
// It also updates the corresponding lattice element in other level 
// regarding which level has the finest level of data 
//
bool LatticeAMR::CheckIn(int level, float x, float y, float z, float t, 
			 float* data) 
{
  int idx = GetIndexinLevel(level, x, y, z, t); 
  if (idx == -1) { printf(" **** no! out of bound.\n"); return false; } // out of bound
  has_data[level][idx] = true; 
  data_ptr[level][idx] = data; 

  // update the blocks at other levels as well
  for (int i=0; i<num_levels; i++){
    idx = GetIndexinLevel(i,x,y,z,t); 
    if (idx != -1) 
      if (finest_level[i][idx]<level)
	finest_level[i][idx] = level; 
  }
  return true; 
}

///////////////////////////////////////////////////////////////////////
//
//  Call this function after all blocks with data have checked in
//  Go through all levels and collect blocks that have data 
//
void LatticeAMR::CompleteLevels()
{
  npart = 0; 
  // first check how many non-empty blocks
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data[i][j]) npart++; 
  // next allocate the volume bounds type etc. 
  vb_list = new volume_bounds_type_f[npart]; 
  parts = new PartitionAMR4D[npart]; 
  rank_to_index = new int[npart]; 

  int isize, jsize, ksize; 
  float x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max; 
  int offset = 0; 
  int rank = 0; 
  for (int level=0; level<num_levels; level++) {
    isize = idim[level]; 
    jsize = jdim[level]; 
    ksize = kdim[level]; 

    x_min = xmin[level]; x_max = xmax[level]; 
    y_min = ymin[level]; y_max = ymax[level]; 
    z_min = zmin[level]; z_max = zmax[level]; 
    t_min = tmin[level]; t_max = tmax[level]; 

    printf(" --- %f %f %f %f %f %f %f %f \n", x_min, x_max, y_min, y_max, 
	   z_min, z_max, t_min, t_max); 
    // j is an id for all blocks in their own level
    // rank is a global id for all non-empty blocks in all levels 
    for (int j=0; j<nblocks[level]; j++) {
      int iidx, jidx, kidx, tidx; 
      if (has_data[level][j]) {
	tidx = j / (isize*jsize*ksize); 
	int r = j % (isize*jsize*ksize); 
	kidx = r / (isize * jsize) ; 
	jidx = (r % (isize * jsize)) / isize; 
	iidx = r % isize; 
	
	vb_list[rank].xmin= x_min+iidx*xlength[level]; 
	vb_list[rank].xmax= x_min+(iidx+1)*xlength[level]; 
	vb_list[rank].ymin= y_min+jidx*ylength[level]; 
	vb_list[rank].ymax= y_min+(jidx+1)*ylength[level]; 
	vb_list[rank].zmin= z_min+kidx*zlength[level]; 
	vb_list[rank].zmax= z_min+(kidx+1)*zlength[level]; 
	vb_list[rank].tmin= t_min+tidx*tlength[level]; 
	vb_list[rank].tmax= t_min+(tidx+1)*tlength[level]; 

	vb_list[rank].xdim = xres[level]; 
	vb_list[rank].ydim = yres[level]; 
	vb_list[rank].zdim = zres[level]; 
	vb_list[rank].tdim = tres[level]; 
	
	// two-way mapping is established 
	index_to_rank[level][j] = rank; 
	rank_to_index[rank] = offset + j; 
	rank++; 
      }
    }
    offset+= nblocks[level]; 
  }
}



///////////////////////////////////////////////////////

void LatticeAMR::CompleteLevels(int t_interval)
{
  // quickly estimate the max number of partitions. 
  // a rough estimate, which is the worse case in that assuming 
  // one time step per partition for a given spatial region 
  npart = 0; 
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data[i][j]) npart++; 

  printf("*** npart = %d\n", npart); 

  // next allocate the volume bounds type etc. 
  vb_list = new volume_bounds_type_f[npart]; 
  parts = new PartitionAMR4D[npart]; 
  rank_to_index = new int[npart]; 

  // below we want to break the time range into segments. the max length
  // of a segment is 't_interval'. A segment will be break if the length 
  // is larger than t_interval, or there is no data in the corresponding 
  // spatial region 
  // 
  npart = -1; 
  int offset = 0; 
  for (int l=0; l<num_levels; l++) {
    for (int k=0; k<kdim[l]; k++) 
      for (int j=0; j<jdim[l]; j++) 
	for (int i=0; i<idim[l]; i++)  // looping through all spatial regions in all levels
	  {
	    int counter = 0; 
	    for (int t=0; t<ldim[l]; t++) {
	      if (counter>=t_interval) {  
		counter = 0;  // reset the counter. ready to start another time segment 
	      }
	      int idx = t*idim[l]*jdim[l]*kdim[l]+k*idim[l]*jdim[l]+j*idim[l]+i; 
	      if (has_data[l][idx]) {
		if (counter == 0) {    // a new time segment starts
		  npart++; 

		  vb_list[npart].xmin= xmin[l]+i*xlength[l]; 
		  vb_list[npart].xmax= xmin[l]+(i+1)*xlength[l]; 
		  vb_list[npart].ymin= ymin[l]+j*ylength[l]; 
		  vb_list[npart].ymax= ymin[l]+(j+1)*ylength[l]; 
		  vb_list[npart].zmin= zmin[l]+k*zlength[l]; 
		  vb_list[npart].zmax= zmin[l]+(k+1)*zlength[l]; 
		  vb_list[npart].tmin= tmin[l]+t*tlength[l]; 
		  vb_list[npart].tmax= tmin[l]+(t+1)*tlength[l]; 

		  vb_list[npart].xdim = xres[l]; 
		  vb_list[npart].ydim = yres[l]; 
		  vb_list[npart].zdim = zres[l]; 
		  vb_list[npart].tdim = 1; 

		  index_to_rank[l][idx] = npart; 
		  rank_to_index[npart] = offset + idx; 

		  counter++; 
		}
		else  {
		  vb_list[npart].tmax= tmin[l]+(t+1)*tlength[l]; 
		  vb_list[npart].tdim += 1; 
		  index_to_rank[l][idx] = npart; 
		  counter++; 
		}
	      }
	      else {  // do not have data
		counter = 0; // reset the counter
	      }
	    }
	  }
    offset+= nblocks[l]; 
  }
  npart+=1; // because earlier we start npart from -1 

}

///////////////////////////////////////////////////////////////////////
//
//  Call this function after all blocks with data have checked in
//  Go through all levels and collect blocks that have data 
//
void LatticeAMR::MergeAndCompleteLevels()
{
  npart = 0; 
  // first check how many non-empty blocks
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data_from_merger[i][j]!=-1) npart++; 
  // next allocate the volume bounds type etc. 
  vb_list = new volume_bounds_type_f[npart]; 
  parts = new PartitionAMR4D[npart]; 
  rank_to_index = new int[npart]; 

  int isize, jsize, ksize; 
  float x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max; 
  int offset = 0; 
  int rank = 0; 
  for (int level=0; level<num_levels; level++) {
    isize = idim[level]; 
    jsize = jdim[level]; 
    ksize = kdim[level]; 

    x_min = xmin[level]; x_max = xmax[level]; 
    y_min = ymin[level]; y_max = ymax[level]; 
    z_min = zmin[level]; z_max = zmax[level]; 
    t_min = tmin[level]; t_max = tmax[level]; 
    // j is to index all blocks in this level
    // rank is to index non-empty blocks in all levels 
    for (int j=0; j<nblocks[level]; j++) {
      int iidx, jidx, kidx, tidx; 
      if (has_data_from_merger[level][j]!=-1) {
	tidx = j / (isize*jsize*ksize); 
	int r = j % (isize*jsize*ksize); 
	kidx = r / (isize * jsize) ; 
	jidx = (r % (isize * jsize)) / isize; 
	iidx = r % isize;
	
	vb_list[rank].xmin= x_min+iidx*xlength[level]; 
	vb_list[rank].xmax= x_min+(iidx+1)*xlength[level]; 
	vb_list[rank].ymin= y_min+jidx*ylength[level]; 
	vb_list[rank].ymax= y_min+(jidx+1)*ylength[level]; 
	vb_list[rank].zmin= z_min+kidx*zlength[level]; 
	vb_list[rank].zmax= z_min+(kidx+1)*zlength[level]; 
	vb_list[rank].tmin= t_min+tidx*tlength[level]; 
	vb_list[rank].tmax= t_min+(tidx+1)*tlength[level]; 

	// need to calculate the resolutions and also
	// update the data pointer  

	int DLevel = has_data_from_merger[level][j]; 
	if (DLevel==level){
	  vb_list[rank].xdim = xres[level]; 
	  vb_list[rank].ydim = yres[level]; 
	  vb_list[rank].zdim = zres[level]; 
	  vb_list[rank].tdim = tres[level]; 
	}
	else {

	  int Imin, Imax, Jmin, Jmax, Kmin, Kmax, Tmin, Tmax; 
	  MapCells(iidx, jidx, kidx, tidx, level, DLevel, Imin, Imax, 
		   Jmin, Jmax, Kmin, Kmax, Tmin, Tmax); 
	  int x_res = xres[DLevel]*(Imax-Imin); 
	  int y_res = yres[DLevel]*(Jmax-Jmin); 
	  int z_res = zres[DLevel]*(Kmax-Kmin); 
	  int t_res = tres[DLevel]*(Tmax-Tmin); 
	  
	  float *data = new float[x_res*y_res*z_res*t_res*3]; 

	  if (data == NULL) {
	    printf(" Panic. cannot allocate memory. \n"); 
	    exit(1); 
	  }
	  printf(" *** %d %d %d %d %d %d %d %d\n", Imin, Imax, Jmin, Jmax, 
		 Kmin, Kmax, Tmin, Tmax); 
	  for (int tFor = Tmin; tFor<Tmax; tFor++) 
	    for (int iFor=Imin; iFor<Imax; iFor++)
	      for (int jFor=Jmin; jFor<Jmax; jFor++)
		for (int kFor=Kmin; kFor<Kmax; kFor++) {
		  int idx = tFor*idim[DLevel]*jdim[DLevel]*kdim[DLevel]+kFor*idim[DLevel]*jdim[DLevel]+jFor*idim[DLevel]+iFor; 
		  // now copy data over 
		  if (data_ptr[DLevel][idx] == NULL) {
		    printf(" Panic. level %d [%d %d %d %d] null data ptr\n", DLevel, 
			   iFor, jFor, kFor, tFor); 
		    exit(1); 
		  }
		  if (has_data[DLevel][idx] == false) {
		    printf(" Panic. level %d has No data!\n", DLevel); 
		    exit(1); 
		  }
		  float* from_data = data_ptr[DLevel][idx]; 
		  for(int dtFor=0; dtFor<tres[DLevel]; dtFor++)
		    for(int dkFor=0; dkFor<zres[DLevel]; dkFor++)
		      for(int djFor=0; djFor<yres[DLevel]; djFor++)
			for(int diFor=0; diFor<xres[DLevel];diFor++) {
			  int source = dtFor*zres[DLevel]*yres[DLevel]*xres[DLevel]+dkFor*yres[DLevel]*xres[DLevel]+djFor*xres[DLevel]+diFor; 
			  int target_i = (iFor-Imin)*xres[DLevel]+diFor; 
			  int target_j = (jFor-Jmin)*yres[DLevel]+djFor; 
			  int target_k = (kFor-Kmin)*zres[DLevel]+dkFor; 
			  int target_t = (tFor-Tmin)*tres[DLevel]+dtFor; 
			  int target = target_t*x_res*y_res*z_res+target_k*x_res*y_res+target_j*x_res+target_i; 
			  data[target*3] = from_data[source*3]; 
			  data[target*3+1] = from_data[source*3+1]; 
			  data[target*3+2] = from_data[source*3+2]; 
			}
		  if (from_data!=NULL) 
		    delete [] from_data; // why cann't you free the data? 
		}

	  data_ptr[level][j] = data; 
	  vb_list[rank].xdim = x_res; 
	  vb_list[rank].ydim = y_res; 
	  vb_list[rank].zdim = z_res; 
	  vb_list[rank].tdim = t_res; 
	  }
	index_to_rank[level][j] = rank; 
	rank_to_index[rank] = offset + j; 
	rank++; 
	}
    }
    offset+= nblocks[level]; 
  }
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++) {
      if (has_data_from_merger[i][j]!=-1) has_data[i][j] = true; 
      else has_data[i][j]=false; 
    }

}

/////////////////////////////////////////////////////////////////
// 
//   Get the data pointer for the block of the given rank 
//

float** LatticeAMR::GetDataPtr(int rank) 
{
  if (rank < 0 || rank >=npart) return (NULL); 

  int n_steps = vb_list[rank].tdim; 
  float ** ppvector = new float*[n_steps]; 

  int remainder = rank_to_index[rank]; 
  int l; 
  for (l=0; l<num_levels; l++) {
    if (remainder-nblocks[l] < 0) break; 
    else 
      remainder -=nblocks[l]; 
  }
  if (l == num_levels) return(NULL); //rank is too big 

  int t_jump = idim[l]*jdim[l]*kdim[l]; 
  int index_offset = 0; 
  for (int i=0; i<n_steps; i++) {
    ppvector[i] = data_ptr[l][remainder+i*t_jump]; 
  }
  return (ppvector); 
}


/////////////////////////////////////////////////////////////////
//
// find the block coords (i,j,k,l) at the given level that contains (x,y,z,t)
// This function does not cheeck whether the element has data or not  
//
int LatticeAMR::GetCoordsinLevel(int level, float x, float y, float z, float t, 
				  int &i, int&j, int&k, int&l)
{
  if (level <0 || level >=num_levels) return -1; 
  if (x<xmin[level] || x>xmax[level] || y<ymin[level] || y>ymax[level] || 
      z<zmin[level] || z>zmax[level] || t<tmin[level] || t>ymax[level])
    return -1; 
  
  i = (int)((x-xmin[level])/xlength[level]); 
  j = (int)((y-ymin[level])/ylength[level]); 
  k = (int)((z-zmin[level])/zlength[level]); 
  l = (int)((t-tmin[level])/tlength[level]); 
  int idx = l*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 
  return idx; 
}

//////////////////////////////////////////////////////////////////
//
//    Return the finest level with data that contains (x,y,z,t) 
//
int LatticeAMR::GetFinestLevel(float x, float y, float z, float t) 
{
  int idx; 
  for (int i=num_levels-1; i>=0; i--) {
    idx= GetIndexinLevel(i, x,y,z,t); 
    if (idx != -1) 
      if (has_data[i][idx] == true) return(i); 
  }
  return(-1); // no data at that point 
}

////////////////////////////////////////////////////////////////////
//
// Find the subdomain rank that contains the physical location (x,y,z,t) 
// with the fineset level of data 
//
int LatticeAMR::GetRank(float x, float y, float z, float t) {

  int idx; 
  for (int i=num_levels-1; i>=0; i--) {
    idx= GetIndexinLevel(i, x,y,z,t); 
    if (idx!=-1)
      if (has_data[i][idx] == true) {
	int rank = index_to_rank[i][idx]; 
	printf("*** rank = %d \n", rank); 
	return rank; 
      }
  }
  return(-1); 
}

//////////////////////////////////////////////////////////////////////
//
//  Get the subdomain rank for the lattice[i,j,k] element at 
//  the given level.
//  note: the lattice element may not contain data. In that case, 
//  the return value will be -1. 
//
int LatticeAMR::GetRank(int i, int j, int k, int t, int level) {

  if (level <0 || level >=num_levels) return -1; 
  if (i<0 || j<0 || k<0 || t<0 ||
      i>=idim[level] || j>=jdim[level] || k>=kdim[level] || t>=ldim[level]) 
    return -1; 
  
  int idx = t*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 

  return index_to_rank[level][idx]; 
}
////////////////////////////////////////////////////////////////////////
//
// find the indices and level of the lattice element that has its 
// subdomain number = rank 
//
int LatticeAMR::GetIndices(int rank, int &iidx, int &jidx, int &kidx, int& tidx, int& level) 
{
  if (rank < 0 || rank >=npart) return (-1); 

  int remainder = rank_to_index[rank]; 
  int l; 
  for (l=0; l<num_levels; l++) {
    if (remainder-nblocks[l] < 0) break; 
    else 
      remainder -=nblocks[l]; 
  }
  if (l == num_levels) return(-1); //rank is too big 

  int isize = idim[l]; 
  int jsize = jdim[l]; 
  int ksize = kdim[l]; 
  int tsize = ldim[l]; 

  tidx = remainder / (isize*jsize*ksize); 
  int r = remainder % (isize*jsize*ksize); 
  kidx = r / (isize * jsize) ; 
  jidx = (r % (isize * jsize)) / isize; 
  iidx = r % isize; 
  level = l;

  return(1); 
}

///////////////////////////////////////////////////////////////
//
// find the indices and level of the lattice element that has 
// the finest resolution of data at (x,y,z,t)
// it also returns the rank of the lattice element 
//
int LatticeAMR::GetIndices(float x, float y, float z, float t, 
			  int &iidx, int &jidx, int &kidx, int&tidx, int&level) 
{
  int idx, rank; 
  for (level=num_levels-1; level>=0; level--) {
    idx= GetCoordsinLevel(level, x,y,z,t, iidx, jidx, kidx,tidx); 
    if (idx !=-1) 
      if (has_data[level][idx] == true)  {
	rank = index_to_rank[level][idx]; 
	return rank; 
      }
  }
  return(-1); // no data at that point 
}

//////////////////////////////////////////////////////////////
//
//check if the point (x,y,z,t) is in the lattice element [i,j,k,l] 
// at the given level 
//
bool LatticeAMR::isIn(float x, float y, float z, float t, int i, int j, int k, 
		      int l, int level) {

  if (level <0 || level >= num_levels) return (-1); 
  if (i < 0 || i > idim[level] - 1 || j < 0 || j > jdim[level] - 1
      || k < 0 || k > kdim[level] - 1 || l <0 || l> ldim[level]-1) 
    return(false); 

  float min_x, max_x, min_y, max_y, min_z, max_z, min_t, max_t; 

  min_x = xmin[level] + i*xlength[level]; 
  max_x = xmin[level] + (i+1)*xlength[level]; 
  min_y = ymin[level] + i*ylength[level]; 
  max_y = ymin[level] + (i+1)*ylength[level]; 
  min_z = zmin[level] + i*zlength[level]; 
  max_z = zmin[level] + (i+1)*zlength[level]; 
  min_t = tmin[level] + i*tlength[level]; 
  max_t = tmin[level] + (i+1)*tlength[level]; 

  if (min_x > x || max_x < x) 
    return (false); 
  else if (min_y > y || max_y < y) 
    return (false); 
  else if (min_z > z || max_z < z) 
    return (false); 
  else if (min_t > t || max_t < t) 
    return(false); 

  return(true); 

}

///////////////////////////////////////////////////////////
//
// look up the volume bounds of lattice[i,j,k,t] at the given level 
//
int LatticeAMR::GetBounds(int i, int j, int k, int t, int level, volume_bounds_type_f& vb)  {

  if (level <0 || level >=num_levels) return (-1); 

  if (i < 0 || i >= idim[level]) 
    return(-1); 
  else if (j < 0 || j >= jdim[level]) 
    return(-1); 
  else if (k < 0 || k >= kdim[level]) 
    return(-1); 
  else if (t < 0 || t >=ldim[level]) 
    return(-1); 

  vb.xmin = xmin[level] + i*xlength[level]; 
  vb.xmax = xmin[level] + (i+1)*xlength[level]; 
  vb.ymin = ymin[level] + i*ylength[level]; 
  vb.ymax = ymin[level] + (i+1)*ylength[level]; 
  vb.zmin = zmin[level] + i*zlength[level]; 
  vb.zmax = zmin[level] + (i+1)*zlength[level]; 
  vb.tmin = tmin[level] + i*tlength[level]; 
  vb.tmax = tmin[level] + (i+1)*tlength[level]; 

  return(1); 
}

//////////////////////////////////////////////////////////////
//
// look up the volume bounds of the subdomain 'rank' 
//
int LatticeAMR::GetBounds(int rank, volume_bounds_type_f &vb)
{
  if (rank < 0 || rank >=npart) return(-1); 

  vb = vb_list[rank];
  return(1); 

}

///////////////////////////////////////////////////////////
//
// returns neighbor rank where x,y,z,t is in 
//
int LatticeAMR::CheckNeighbor(int myrank, float x, float y, float z, float t) {

  return GetRank(x,y,z,t); 

}
////////////////////////////////////////////////////////
//
// GetNeighbor
//
// Returns the indices and level of the neighbor that has data and contains the point (x,y,z,t)
// Also returns the rank of the element (-1 means no data found at the point in all levels)
// myrank: global subdomain id (not used any more) 
//
//
int LatticeAMR::GetNeighbor(int myrank, float x, float y, float z, float t, 
			   int &ei, int &ej, int &ek, int &et, int &level) {


  int idx = GetIndices(x,y,z,t,ei, ej, ek, et, level); 

  return idx; 

}

/////////////////////////////////////////////////////////////////
//
//   Map a cell from one level to the cell(s) in the over level 
// 
bool LatticeAMR::MapCells(int fromI, int fromJ, int fromK, int fromT, 
		 int fromLevel, int toLevel, int& toImin, int& toImax, 
		 int& toJmin, int& toJmax, int& toKmin, int& toKmax, 
		 int& toTmin, int& toTmax)
{
  float minX = xmin[fromLevel]+fromI*xlength[fromLevel]; 
  float maxX = xmin[fromLevel]+(fromI+1)*xlength[fromLevel]; 
  float minY = ymin[fromLevel]+fromJ*ylength[fromLevel]; 
  float maxY = ymin[fromLevel]+(fromJ+1)*ylength[fromLevel]; 
  float minZ = zmin[fromLevel]+fromK*zlength[fromLevel]; 
  float maxZ = zmin[fromLevel]+(fromK+1)*zlength[fromLevel]; 
  float minT = tmin[fromLevel]+fromT*tlength[fromLevel]; 
  float maxT = tmin[fromLevel]+(fromT+1)*tlength[fromLevel]; 

  toImin = (int)((minX - xmin[toLevel])/(float)xlength[toLevel]); 
  if (toImin<0 || toImin>idim[toLevel]) return false; 
  toImax = (int)((maxX - xmin[toLevel])/(float)xlength[toLevel]); 
  if (toImax<0 || toImax>idim[toLevel]) return false; 
  toJmin = (int)((minY - ymin[toLevel])/(float)ylength[toLevel]); 
  if (toJmin<0 || toJmin>jdim[toLevel]) return false; 
  toJmax = (int)((maxY - ymin[toLevel])/(float)ylength[toLevel]); 
  if (toJmax<0 || toJmax>jdim[toLevel]) return false; 
  toKmin = (int)((minZ - zmin[toLevel])/(float)zlength[toLevel]); 
  if (toKmin<0 || toKmin>kdim[toLevel]) return false; 
  toKmax = (int)((maxZ - zmin[toLevel])/(float)zlength[toLevel]); 
  if (toKmax<0 || toKmax>kdim[toLevel]) return false; 
  toTmin = (int)((minT - tmin[toLevel])/(float)tlength[toLevel]); 
  if (toTmin<0 || toTmin>ldim[toLevel]) return false; 
  toTmax = (int)((maxT - tmin[toLevel])/(float)tlength[toLevel]); 
  if (toTmax<0 || toTmax>ldim[toLevel]) return false; 
  return true; 
}

//////////////////////////////////////////////////////////////////////
//
//    To see if data from blocks of a higher LOD can be 
//    merged into the block (i,j,k,t) at this level
//
//    The purpose of merger is to reduce the number of blocks in the 
//    data set for efficiency 
//    
//
bool LatticeAMR::Mergeable(int i, int j, int k, int t, int level, 
			   int& mergeLevel)
{
  int toImin, toImax, toJmin, toJmax, toKmin, toKmax, toTmin, toTmax; 

  int idx = t*idim[level]*jdim[level]*kdim[level]+k*idim[level]*jdim[level]+j*idim[level]+i; 

  // Blocks from which level (as deep as possible) can be merged to this level? 
  mergeLevel = level; 
  if (has_data[level][idx]==true) { //the block has data. no merging needed 
    has_data_from_merger[level][idx] = level; //default: its own level 
    return true; 
  }
  else // check finer levels
    if (level+1 >= num_levels) // this is the finest level 
      return false; // this means in this region no data is available 
    else {
      // Can the space of this block be filled by blocks in the next level? 
      // if yes, what are the indices of those next level blocks 
      bool mappable = MapCells(i,j,k,t,level,level+1, toImin, toImax, toJmin, toJmax, 
		toKmin, toKmax, toTmin, toTmax); 
      if (mappable == false) {
	return false; // out of next level's bound 
      }

      bool first = true; 
      for (int tFor = toTmin; tFor<toTmax; tFor++) 
	for (int kFor = toKmin; kFor<toKmax; kFor++)
	  for (int jFor = toJmin; jFor<toJmax; jFor++)
	    for (int iFor = toImin; iFor<toImax; iFor++) {
	      // recursively call the next level 
	      int mLevel; 
	      if (Mergeable(iFor, jFor, kFor, tFor, level+1, 
			    mLevel)== false) {
		return(false); 
	      }
	      else {
		if (first == true)  {
		  mergeLevel = mLevel; 
		  first = false; 
		} 
		else  // all the merge levels from children have to be the same 
		  if (mLevel != mergeLevel) {
		    return(false); 
		  }
	      }
	    }
      // mergeable 
      // reset the children's record 
      for (int tFor = toTmin; tFor<toTmax; tFor++) 
	for (int iFor = toImin; iFor<toImax; iFor++)
	  for (int jFor = toJmin; jFor<toJmax; jFor++)
	    for (int kFor = toKmin; kFor<toKmax; kFor++){
	      int idx2 = tFor*idim[level+1]*jdim[level+1]*kdim[level+1]+kFor*idim[level+1]*jdim[level+1]+jFor*idim[level+1]+iFor; 
	      has_data_from_merger[level+1][idx2] = -1; 
	    }
      has_data_from_merger[level][idx]=mergeLevel; 
      return(true); 
    }
}
////////////////////////////////////////////////////////////
//
//  Merge smaller blocks of the same resolutions into a larger 
//  block whenever possible 
//

void LatticeAMR::MergeBlocks()
{
  int l = 0; // 0-th level;  coarsest resolution 
             // other levels will be done through recursive calls 
  for (int tFor=0; tFor<ldim[l]; tFor++)
    for (int kFor=0; kFor<kdim[l]; kFor++)
      for (int jFor=0; jFor<jdim[l]; jFor++)
	for (int iFor=0; iFor<idim[l]; iFor++) {
	  int mergeLevel; 
	  bool mg = Mergeable(iFor, jFor, kFor, tFor, l, mergeLevel); 
	  if (mg == true) 
	    printf("%d %d %d %d got level %d merger\n", iFor, jFor, kFor, tFor, 
		   mergeLevel); 
	}

  int cnt = 0; 
  for (int i=0; i<num_levels; i++)
    for (int j=0; j<nblocks[i]; j++)
      if (has_data_from_merger[i][j]!=-1) cnt++; 
  printf(" cnt = %d \n", cnt); 

}

//----------------------------------------------------------------------------

void LatticeAMR::InitSeedLists() {

  seedlists = new list<VECTOR4>[npart]; 

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

void LatticeAMR::ResetSeedLists() {

  for (int i = 0; i < npart; i++)
    seedlists[i].clear(); 

}

void LatticeAMR::ResetSeedLists(int i) {

  seedlists[i].clear(); 

}
//--------------------------------------------------------------------------

bool LatticeAMR::InsertSeed(int i, int j, int k, int t, int level, VECTOR4 p) {

  int rank = GetRank(i,j,k, t, level); 

  //  printf(" insert to rank %d \n", rank); 
  if (rank ==-1) return(false); 
  else {
    seedlists[rank].push_back(p); 
    return(true); 
  }
}


bool LatticeAMR::InsertSeed(int i, VECTOR4 p) {

  if (i<0 || i>=npart) return(false); 
  else {
    seedlists[i].push_back(p); 
    return(true); 
  }
}


//---------------------------------------------------------------------------
//
// assign the partitions to the processors in a round-robin manner 
//
void LatticeAMR::RoundRobin_proc(int nproc) {

  for (int i = 0; i < npart; i++) 
    parts[i].Proc = i % nproc; 

}

//---------------------------------------------------------------------------
//
int LatticeAMR::GetProc(int rank) {

  return(parts[rank].Proc); 

}

//---------------------------------------------------------------------------
//
// look up the lattice[i,j,k] element at the given level 
//  for the processor of the subdomain
//
int LatticeAMR::GetProc(int i, int j, int k, int t, int level) 
{
  if (level <0 || level >=num_levels) return (-1); 

  if (i < 0 || i >= idim[level])
    return(-1); 
  else if (j < 0 || j >= jdim[level])
    return(-1); 
  else if (k < 0 || k >= kdim[level])
    return(-1); 
  else if (t < 0 || t >= ldim[level])
    return(-1); 

  int rank = GetRank(i,j,k,t,level); 
  if (rank == -1) return (-1); 
  else 
    return (parts[rank].Proc); 
}

//----------------------------------------------------------------------------
//
// query the paritions that are assigned to processor 'proc' 
//
void LatticeAMR::GetPartitions(int proc, int**p_list, int& num) {

  int cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      cnt++;
  }

  num = cnt; 
  (*p_list) = new int[cnt]; 
  cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      (*p_list)[cnt++] = i; 
  }

}
//----------------------------------------------------------------------------
//
// GetPartitions
//
// query the paritions that are assigned to processor 'proc' 
//
// p_list must be allocated large enough prior to calling
//
// Tom Peterka, 12/5/08
//
void LatticeAMR::GetPartitions(int proc, int*p_list, int& num) {

  int cnt = 0; 

  for (int i = 0; i < npart; i++)
    if (parts[i].Proc == proc)
      cnt ++; 

  num = cnt; 
  cnt = 0; 

  for (int i = 0; i < npart; i++) {
    if (parts[i].Proc == proc)
      p_list[cnt++] = i;
  }

}
