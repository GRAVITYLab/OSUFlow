
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

  // the physical range of each level 
  xmin = new float[num_levels];   xmax = new float[num_levels]; 
  ymin = new float[num_levels];   ymax = new float[num_levels]; 
  zmin = new float[num_levels];   zmax = new float[num_levels]; 
  tmin = new float[num_levels];   tmax = new float[num_levels]; 

  // the physical size of the blocks in each level 
  xlength = new float[num_levels]; 
  ylength = new float[num_levels]; 
  zlength = new float[num_levels]; 
  tlength = new float[num_levels]; 

  // the dimensions of the lattice in each level 
  idim = new int[num_levels]; 
  jdim = new int[num_levels]; 
  kdim = new int[num_levels]; 
  ldim = new int[num_levels]; 

  // data structures for book-keeping 
  has_data = new bool*[num_levels];     // whether the element has data or not 
  finest_level = new int*[num_levels];  // what is the finest level of data 
                                         // in this region
  data_ptr = new float**[num_levels]; 
  index_to_rank = new int*[total_level]; // mapping from index to rank 
  nblocks = new int[total_level];        // how many blocks in each region 
                                         // regardless of empty or not 
  npart = 0;                            // number of number empty blocks total
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
bool LatticeAMR::CreateLevel(int level, float x_size, float y_size, 
			     float z_size, float t_size, 
			     float x_min, float x_max, 
			     float y_min, float y_max, float z_min, 
			     float z_max, float t_min, float t_max) 
{

  if (level <0 || level >=num_levels) return false; 

  xlength[level]=x_size; ylength[level]=y_size; zlength[level]=z_size;
 
  xmin[level] = x_min; xmax[level] = x_max; 
  ymin[level] = y_min; ymax[level] = y_max; 
  zmin[level] = z_min; zmax[level] = z_max;
  tmin[level] = t_min; tmax[level] = t_max; 
  
  idim[level] = (int) (x_max-x_min)/x_size; 
  jdim[level] = (int) (y_max-y_min)/y_size; 
  kdim[level] = (int) (z_max-z_min)/z_size; 
  ldim[level] = (int) (t_max-t_min)/t_size; 

  printf("**** level %d dims %d %d %d %d \n", level, idim[level], jdim[level], 
	 kdim[level], ldim[level]); 

  int size = idim[level]*jdim[level]*kdim[level]*ldim[level]; 
  nblocks[level] = size;  // number of blocks in this level 

  has_data[level] = new bool[size]; 
  finest_level[level] = new int[size]; 
  data_ptr[level] = new float*[size]; 
  index_to_rank[level] = new int[size];

  int idx = 0; 
  for (int t=0; t<ldim[level]; t++) 
    for (int k=0; k<kdim[level]; k++) 
      for (int j=0; j<jdim[level]; j++) 
	for (int i=0; i<idim[level]; i++) {
	  has_data[level][idx] = false;    // initial value: no data
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
// It also updates the corresponding element in other level 
//
bool LatticeAMR::CheckIn(int level, float x, float y, float z, float t, 
			 float* data) 
{
  int idx = GetIndexinLevel(level, x, y, z, t); 
  if (idx == -1) return false; 
  has_data[level][idx] = true; 
  data_ptr[level][idx] = data; 

  // update other blocks at other levels as well
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
  vb_list = new volume_bounds_type[npart]; 
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
	
	index_to_rank[level][j] = rank; 
	rank_to_index[rank] = offset + j; 

	rank++; 
      }
    }
    offset+= nblocks[level]; 
  }
}

/////////////////////////////////////////////////////////////////
// 
//   Get the data pointer for the block of the given rank 
//
float* LatticeAMR::GetDataPtr(int rank) 
{
  if (rank < 0 || rank > npart) return (NULL); 

  int remainder = rank_to_index[rank]; 
  int l; 
  for (l=0; l<num_levels; l++) {
    if (remainder-nblocks[l] < 0) break; 
    else 
      remainder -=nblocks[l]; 
  }
  if (l == num_levels) return(NULL); //rank is too big 
  return (data_ptr[l][remainder]); 
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
  if (rank < 0 || rank > npart) return (-1); 

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
// the finest data at (x,y,z,t)
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
  max_x = xmax[level] + (i+1)*xlength[level]; 
  min_y = ymin[level] + i*ylength[level]; 
  max_y = ymax[level] + (i+1)*ylength[level]; 
  min_z = zmin[level] + i*zlength[level]; 
  max_z = zmax[level] + (i+1)*zlength[level]; 
  min_t = tmin[level] + i*tlength[level]; 
  max_t = tmax[level] + (i+1)*tlength[level]; 

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
int LatticeAMR::GetBounds(int i, int j, int k, int t, int level, volume_bounds_type& vb)  {

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
  vb.xmax = xmax[level] + (i+1)*xlength[level]; 
  vb.ymin = ymin[level] + i*ylength[level]; 
  vb.ymax = ymax[level] + (i+1)*ylength[level]; 
  vb.zmin = zmin[level] + i*zlength[level]; 
  vb.zmax = zmax[level] + (i+1)*zlength[level]; 
  vb.tmin = tmin[level] + i*tlength[level]; 
  vb.tmax = tmax[level] + (i+1)*tlength[level]; 

  return(1); 
}

//////////////////////////////////////////////////////////////
//
// look up the volume bounds of the subdomain 'rank' 
//
int LatticeAMR::GetBounds(int rank, volume_bounds_type &vb)
{
  if (rank < 0 || rank >npart) return(-1); 

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
// returns neighbor that contains x,y,z,t point with data 
// and sets i,j,k of neighbor
// it also returns the rank of the element (-1 means no data there)
// myrank: global partition number
//
//
int LatticeAMR::GetNeighbor(int myrank, float x, float y, float z, float t, 
			   int &ei, int &ej, int &ek, int &et, int &level) {


  int idx = GetIndices(x,y,z,t,ei, ej, ek, et, level); 

  return idx; 

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
