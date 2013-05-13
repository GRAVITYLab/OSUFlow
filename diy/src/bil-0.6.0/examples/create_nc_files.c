/***************************************************************************
 *   Copyright (C) 2010  Wes Kendall                                       *
 *   kendall@eecs.utk.edu                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 3 of the    *
 *   License, or (at your option) any later version.                       *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Lesser General Public License for more details.                   *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this program; if not, write to the                 *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// create_nc_files.c
// Wes Kendall
// 05/29/2010 
// creates files of 3 dimensions and fills it with 
// increasing integers for testing purposes. the user supplies
// 3 arguments to the program: the dimension size to use for (the same
// size will be used for each dimension), the number of files to generate,
// and the directory to store the files in (this argument is optional and
// defaults to the user's local directory)
////

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <netcdf.h>
#include <string.h>
#include "bil.h"

int main(int argc, char **argv)
{
	if (argc < 3 || argc > 4)
	{
		fprintf(stderr, "usage %s dim_size num_files <file_dir>\n", argv[0]);
		exit(1);
	}
	const int dim_size = atoi(argv[1]);
	const int num_files = atoi(argv[2]);
	const char *file_dir = (argc == 4) ? argv[3] : "./";

	// do some checking of arguments
	assert(strlen(file_dir) < 512);
	assert(dim_size <= 512 && dim_size > 0);
	assert(num_files <= 16 && num_files > 0);
	
	// allocate enough room to store one variable
	int *data = malloc(sizeof(int) * dim_size * dim_size * dim_size);
	assert(data != NULL);
	
	int i;
	for (i = 0; i < num_files; i++)
	{
		// create a file with a default type of file name so the 
		// testing program can read it in
		char file_name[1024];
		sprintf(file_name, "%s/%d_%d_%d.%d.nc", file_dir,
				dim_size, dim_size, dim_size, i);

		int nc_file;
		int dim_ids[3], var_ids[2];
		// create the file with 3 dimensions and 2 variables
		assert(nc_create(file_name, 0, &nc_file) == NC_NOERR);
		assert(nc_def_dim(nc_file, "dim0", dim_size, &(dim_ids[0])) == NC_NOERR);
		assert(nc_def_dim(nc_file, "dim1", dim_size, &(dim_ids[1])) == NC_NOERR);
		assert(nc_def_dim(nc_file, "dim2", dim_size, &(dim_ids[2])) == NC_NOERR);
		assert(nc_def_var(nc_file, "var0", NC_INT, 3, dim_ids, &(var_ids[0]))
				== NC_NOERR);
		assert(nc_def_var(nc_file, "var1", NC_INT, 3, dim_ids, &(var_ids[1]))
				== NC_NOERR);
		nc_enddef(nc_file);

		// write the variable data. make the variable data offset by which file
		// it is in and increase the count by one by the highest varying dimension
		// (the x dimension)
		size_t start[3] = {0, 0, 0};
		size_t size[3] = {dim_size, dim_size, dim_size};
		int j;
		for (j = 0; j < dim_size * dim_size * dim_size; j++)
		{
			data[j] = j + (dim_size * dim_size * dim_size * i);
		}
		assert(nc_put_vara_int(nc_file, var_ids[0], start, size, data) == NC_NOERR);

		// write a second variable and simply add 1 to it to test reading
		// of multiple variables
		for (j = 0; j < dim_size * dim_size * dim_size; j++)
		{
			data[j] = j + (dim_size * dim_size * dim_size * i) + 1;
		}
		assert(nc_put_vara_int(nc_file, var_ids[1], start, size, data) == NC_NOERR);
		nc_close(nc_file);
	}

	free(data);
	return 0;
}
