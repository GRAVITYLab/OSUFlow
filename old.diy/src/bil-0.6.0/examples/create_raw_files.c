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

// create_raw_files.c
// Wes Kendall
// 05/29/2010
// creates raw files of 3 dimensions and fills them with 
// increasing integers for testing purposes. the program takes
// 3 inputs: the dimension size (this size will be used for all
// dimensions), the number of files to generate, and an optional 
// directory to store the files in (defaults to the user's current
// directory)
////

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

	// allocate enough room for a variable
	int *data = malloc(sizeof(int) * dim_size * dim_size * dim_size);
	assert(data != NULL);

	int i;
	for (i = 0; i < num_files; i++)
	{
		// create a file named accordingly so that the test program can
		// read it in
		char file_name[1024];
		sprintf(file_name, "%s/%d_%d_%d.%d.raw", file_dir, 
				dim_size, dim_size, dim_size, i);

		// open the file
		FILE *fp = fopen(file_name, "wb");
		assert(fp != NULL);
		// write the variable to the file. generate monotonically increasing
		// numbers that start at the ending of the previous file's numbers. This
		// is for testing purposes
		int j;
		for (j = 0; j < dim_size * dim_size * dim_size; j++)
		{
			data[j] = j + (dim_size * dim_size * dim_size * i);
		}
		fwrite(data, dim_size * dim_size * dim_size, sizeof(int), fp);
		fclose(fp);
	}

	// clean up
	free(data);
	return 0;
}
