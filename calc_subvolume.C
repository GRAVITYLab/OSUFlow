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

volume_bounds_type* vb_list; 

//////////////////////////////////////////////////////////////////////////

int split_volume
  ( 
   int level,
   volume_bounds_type *vb, 
   volume_bounds_type *vb1, 
   volume_bounds_type *vb2 
  )
{
  int split = 0;

  *vb1 = *vb2 = *vb;
  switch (level % 3)
    {
    case 0:  
      /* split in X */
      if (vb->xmin != vb->xmax-1)
	{
	  vb2->xmin = vb1->xmax = (vb->xmin+vb->xmax+(level%2)) / 2;
	  split = 'X';
	  break;
	}
    case 1:
      /* split in Y */
      if (vb->ymin != vb->ymax-1)
	{
	  vb2->ymin = vb1->ymax = (vb->ymin+vb->ymax+(level%2)) / 2;
	  split = 'Y';
	  break;
	}
    case 2:
      /* split in Z */
      if (vb->zmin != vb->zmax-1)
	{
	  vb2->zmin = vb1->zmax = (vb->zmin+vb->zmax+(level%2)) / 2;
	  split = 'Z';
	  break;
	}
    }
  if (!split && 
      (vb->xmin != vb->xmax+1 || 
       vb->ymin != vb->ymax+1 || 
       vb->zmin != vb->zmax+1))
    {
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
static void calc_subvolume_h
  (
   volume_bounds_type vb,	/* the subvolume to divide */

   int pmin,			/* The processors available */
   int pmax,		

   int level			/* recursion level */
  )
{
  volume_bounds_type vb1,vb2;	/* divided subvolumes sub */
  int pmid;			/* dividing point of processors */
  int split;			/* boolean: true if split was possible */
  int split_volume(int,volume_bounds_type*,volume_bounds_type*,
		    volume_bounds_type*);

  /* 
   * Sanity check arguments 
   */
  assert(pmin <= pmax && 
	 vb.xmin <= vb.xmax && 
	 vb.ymin <= vb.ymax && 
	 vb.zmin <= vb.zmax );

  /*
   * If we are down to one processor, it is the base case.  
   */
  if (pmin == pmax)
    {
      /* 
       *  Send off the subvolume processor pmin is to work on.
       */ 
      
      //      printf(" split volume [%d] .... %d %d %d %d %d %d\n", 
      //	     pmin, vb.xmin, vb.xmax, 
      //	     vb.ymin, vb.ymax, vb.zmin, vb.zmax); 

      vb_list[pmin].xmin = vb.xmin;  vb_list[pmin].xmax = vb.xmax; 
      vb_list[pmin].ymin = vb.ymin;  vb_list[pmin].ymax = vb.ymax; 
      vb_list[pmin].zmin = vb.zmin;  vb_list[pmin].zmax = vb.zmax; 

      
      //      MPI_Send(&vb, sizeof(vb), MPI_BYTE, pmin, 0, MPI_COMM_WORLD); 

      return;
    }

  split = split_volume( level, &vb, &vb1, &vb2 );

  if (split)
    {
      /* Recurse on subproblems */
      pmid = (pmax + pmin)/2;
      calc_subvolume_h( vb1, pmin, pmid, level+1 );
      calc_subvolume_h( vb2, pmid+1, pmax, level+1 );
    }
  else
    {
      /* Couldn't split the problem.  */
      printf( "Too many processors for this problem size!\n" );
    }
}

//////////////////////////////////////////////////
volume_bounds_type*
calc_subvolume(int vxdim,int vydim,int vzdim, int nproc)
{
  volume_bounds_type vb;
  vb_list = new volume_bounds_type[nproc]; 

  /* Call recursive procedure that does the real work */
  vb.xmin = vb.ymin = vb.zmin = 0;
  vb.xmax = vxdim-1;
  vb.ymax = vydim-1;
  vb.zmax = vzdim-1;
  calc_subvolume_h(vb, 0, nproc-1, 0);
  return (vb_list); 
}



