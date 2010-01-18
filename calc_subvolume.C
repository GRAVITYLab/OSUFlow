/*
 *  cm5 volume renderer: calc_subvolume
 *
 * Purpose:
 *
 * Calculate the portion of the volume that each processor will be
 * dealing with and send it a message describing the bounds.
 *
 * James S. Painter (painter@cs.utah.edu)     August 15 1992
 * Univesity of Utah, Computer Science Department 
 *
 * Copyright (c) 1992 University of Utah
 * Copying, use and development for non-commercial purposes permitted.
 *                  All rights for commercial use reserved.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#include "calc_subvolume.h"

static volume_bounds_type* vb_list;  
static volume_bounds_type* vb_list_adj;  

int x_list[1000], y_list[1000], z_list[1000]; 
int x_cnt =0, y_cnt=0, z_cnt =0; 

// a stupid list sort. not the fastest algorithm. but get the job done.
void sort_list(int* list, int& cnt) 
{
  int* tmp_list = new int[cnt]; 
  for (int i=0; i<cnt; i++) tmp_list[i] = list[i]; 
  for (int i=0; i<cnt; i++) {
    int largest = -1; 
    int idx = i; 
    for (int j = i; j<cnt; j++) {
      if (tmp_list[j]>largest) {idx = j; largest = tmp_list[j]; }
    }
    int tmp = tmp_list[i]; 
    tmp_list[i] = tmp_list[idx]; 
    tmp_list[idx] =tmp; 
  }
  list[0] = tmp_list[cnt-1]; 
  int count = 1; 
  for (int i=cnt-2; i>=0; i--) {
    if (tmp_list[i]>list[count-1]) {
      list[count] = tmp_list[i]; 
      count++; 
    }
  }
  cnt = count; 

  delete[] tmp_list;

}

//////////////////////////////////////////////////////////////////////////

int split_volume ( 
		  int level,
		  volume_bounds_type *vb, 
		  volume_bounds_type *vb1, 
		  volume_bounds_type *vb2) {

  int split = 0;

  *vb1 = *vb2 = *vb;
  switch (level % 3)
    {
    case 0:  /* split in X */
      if (vb->xmin != vb->xmax-1) {
	vb2->xmin = vb1->xmax = (vb->xmin+vb->xmax+(level%2)) / 2;
	x_list[x_cnt++] = vb2->xmin; 
	split = 'X';
	break;
      }
    case 1: /* split in Y */
      if (vb->ymin != vb->ymax-1) {
	vb2->ymin = vb1->ymax = (vb->ymin+vb->ymax+(level%2)) / 2;
	y_list[y_cnt++] = vb2->ymin; 
	split = 'Y';
	break;
      }
    case 2: /* split in Z */
      if (vb->zmin != vb->zmax-1) {
	vb2->zmin = vb1->zmax = (vb->zmin+vb->zmax+(level%2)) / 2;
	z_list[z_cnt++] = vb2->zmin; 
	split = 'Z';
	break;
      }
    }

  if (!split && 
      (vb->xmin != vb->xmax+1 || 
       vb->ymin != vb->ymax+1 || 
       vb->zmin != vb->zmax+1)) {
    /* We failed to split for some reason.  Increment level and try again 
     * to force a split on the next dimension.  This will happen only
     * if we "bottom out" in one dimension before the others -- unlikely.
     */
    split = split_volume( level-1, vb, vb1, vb2 );
  }

  return split;

}
///////////////////////////////////////////////////////////////////////////////

/*
 *  Recursive helper function for calc_subvolume.
 */
static void calc_subvolume_h (
			      volume_bounds_type vb, /* the subvolume to divide */
			      int pmin,	/* The processors available */
			      int pmax,		
			      int level) { /* recursion level */

  volume_bounds_type vb1,vb2;	/* divided subvolumes sub */
  int pmid;			/* dividing point of processors */
  int split;			/* boolean: true if split was possible */
  int split_volume(int,volume_bounds_type*,volume_bounds_type*,
		   volume_bounds_type*);

  /* Sanity check arguments */
  assert(pmin <= pmax && 
	 vb.xmin <= vb.xmax && 
	 vb.ymin <= vb.ymax && 
	 vb.zmin <= vb.zmax );

  /* If we are down to one processor, it is the base case.  */
  if (pmin == pmax) {
    vb_list[pmin].xmin = vb.xmin;  vb_list[pmin].xmax = vb.xmax; 
    vb_list[pmin].ymin = vb.ymin;  vb_list[pmin].ymax = vb.ymax; 
    vb_list[pmin].zmin = vb.zmin;  vb_list[pmin].zmax = vb.zmax; 
    return;
  }

  split = split_volume( level, &vb, &vb1, &vb2 );

  if (split) {      /* Recurse on subproblems */
    pmid = (pmax + pmin)/2;
    calc_subvolume_h( vb1, pmin, pmid, level+1 );
    calc_subvolume_h( vb2, pmid+1, pmax, level+1 );
  }
  else       /* Couldn't split the problem.  */
    printf( "Too many processors for this problem size!\n" );

}
//////////////////////////////////////////////////
volume_bounds_type* calc_subvolume(int vxdim, int vydim, int vzdim, 
				   int ghost, int npart, 
				   int &lattice_xdim,
				   int& lattice_ydim, int& lattice_zdim) {

  volume_bounds_type vb;
  vb_list = new volume_bounds_type[npart]; 
  int* lattice; 

  /* Call recursive procedure that does the real work */
  vb.xmin = vb.ymin = vb.zmin = 0;
  vb.xmax = vxdim - 1;
  vb.ymax = vydim - 1;
  vb.zmax = vzdim - 1;
  x_list[x_cnt++] = 0; y_list[y_cnt++] = 0; z_list[z_cnt++] = 0; 
  calc_subvolume_h(vb, 0, npart - 1, 0);

  //   for (int i=0; i<x_cnt; i++) printf("x[%d] = %d ", i, x_list[i]); 
  //   printf("\n"); 
  //   for (int i=0; i<y_cnt; i++) printf("y[%d] = %d ", i, y_list[i]); 
  //   printf("\n"); 
  //   for (int i=0; i<z_cnt; i++) printf("z[%d] = %d ", i, z_list[i]); 
  //   printf("\n"); 

  sort_list(x_list, x_cnt); 
  sort_list(y_list, y_cnt); 
  sort_list(z_list, z_cnt); 

    for (int i=0; i<x_cnt; i++) printf("x[%d] = %d ", i, x_list[i]); 
    printf("\n"); 
    for (int i=0; i<y_cnt; i++) printf("y[%d] = %d ", i, y_list[i]); 
    printf("\n"); 
    for (int i=0; i<z_cnt; i++) printf("z[%d] = %d ", i, z_list[i]); 
    printf("\n"); 

  lattice = new int[x_cnt*y_cnt*z_cnt]; 
  int xidx, yidx, zidx; 
  for (int i=0; i<npart; i++) {

    // create a lattice for the subdomains

    for (int j=0; j<x_cnt; j++) 
      if (vb_list[i].xmin == x_list[j]) xidx = j; 
    for (int j=0; j<y_cnt; j++) 
      if (vb_list[i].ymin == y_list[j]) yidx = j; 
    for (int j=0; j<z_cnt; j++) 
      if (vb_list[i].zmin == z_list[j]) zidx = j; 
    //     printf("vblist %d [%d %d %d] at lattice [%d %d %d]\n", i, vb_list[i].xmin, 
    // 	   vb_list[i].ymin, vb_list[i].zmin, zidx, yidx, xidx); 

    lattice[zidx*x_cnt*y_cnt+yidx*x_cnt+xidx] = i; 

    // add ghost dimensions
    if (vb_list[i].xmin-ghost >=0) vb_list[i].xmin-=ghost; 
    if (vb_list[i].xmax+ghost <vxdim) vb_list[i].xmax+=ghost; 
    if (vb_list[i].ymin-ghost >=0) vb_list[i].ymin-=ghost; 
    if (vb_list[i].ymax+ghost <vydim) vb_list[i].ymax+=ghost; 
    if (vb_list[i].zmin-ghost >=0) vb_list[i].zmin-=ghost; 
    if (vb_list[i].zmax+ghost <vzdim) vb_list[i].zmax+=ghost; 

  }
  lattice_xdim = x_cnt; 
  lattice_ydim = y_cnt; 
  lattice_zdim = z_cnt; 

  // readjusting the index. so that lattice[i][j][k] stores the k*xdim*ydim+j*xdim+i -th partition 
  // without the following code, it is not necessary the case 
  vb_list_adj = new volume_bounds_type[npart]; 
  for (int i=0; i<lattice_xdim; i++) 
    for (int j=0; j<lattice_ydim; j++) 
      for (int k=0; k<lattice_zdim; k++) {
	int idx = k * lattice_ydim*lattice_xdim + j * lattice_xdim + i; 
	vb_list_adj[idx].xmin = vb_list[lattice[idx]].xmin; 
	vb_list_adj[idx].xmax = vb_list[lattice[idx]].xmax; 
	vb_list_adj[idx].ymin = vb_list[lattice[idx]].ymin; 
	vb_list_adj[idx].ymax = vb_list[lattice[idx]].ymax; 
	vb_list_adj[idx].zmin = vb_list[lattice[idx]].zmin; 
	vb_list_adj[idx].zmax = vb_list[lattice[idx]].zmax; 
      }
  
  delete [] vb_list; 
  delete[] lattice;
  return (vb_list_adj); 

}



