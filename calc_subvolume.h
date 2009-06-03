/*
 *  cm5 volume renderer: calc_subvolume
 *
 * Purpose:

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

#ifndef _CALC_SUBVOLUME_H
#define _CALC_SUBVOLUME_H 

/* Data Types */
typedef
  struct volume_bounds_tag	/* Bounds within a volume */
{
  int xmin, xmax; 
  int ymin, ymax;
  int zmin, zmax;
  int tmin, tmax;   // computational time steps
  int xdim, ydim, zdim, tdim; // data resolution 
} volume_bounds_type;

typedef
  struct volume_bounds_tag_f	/* Bounds within a volume */
{
  float xmin, xmax; 
  float ymin, ymax;
  float zmin, zmax;
  int tmin, tmax;   // computational time steps
  int xdim, ydim, zdim, tdim; // data resolution 
} volume_bounds_type_f;

/* Function Prototypes */

/*
 * Given the total data volume size, split the volume up into pieces
 * and send a message to each processing node informing it what part of
 * the volume it will be dealing with.
 */ 
volume_bounds_type* calc_subvolume(int vxdim, int vydim, int vzdim, 
				   int ghost, int nproc, 
				   int& l_xdim, int& l_ydim, int& l_zdim);

#endif
