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

// bil_pio.c
// Wes Kendall
// 07/18/2010
// Functions for aggregating and issuing I/O in BIL.
////

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#ifndef DISABLE_PNETCDF
#include <pnetcdf.h>
#endif // DISABLE_PNETCDF

#include "bil.h"
#include "bil_pio.h"
#include "bil_misc.h"
#include "bil_block.h"
#include "bil_sched.h"
#include "bil_timing.h"

int BIL_Pio_compare_block_id(const void* a, const void* b) {
  BIL_Block* block_a = (BIL_Block*)a;
  BIL_Block* block_b = (BIL_Block*)b;
  if (block_a->block_id < block_b->block_id) {
    return -1;
  } else if (block_a->block_id > block_b->block_id) {
    return 1;
  } else {
    return 0;
  }
}

int BIL_Pio_compare_block_read_rank(const void* a, const void* b) {
  BIL_Block* block_a = (BIL_Block*)a;
  BIL_Block* block_b = (BIL_Block*)b;
  if (block_a->read_rank < block_b->read_rank) {
    return -1;
  } else if (block_a->read_rank > block_b->read_rank) {
    return 1;
  } else {
    if (block_a->block_id < block_b->block_id) {
      return -1;
    } else if (block_a->block_id > block_b->block_id) {
      return 1;
    } else {
      return 0;
    }
  }
}

int BIL_Pio_compare_block_request_rank(const void* a, const void* b) {
  BIL_Block* block_a = (BIL_Block*)a;
  BIL_Block* block_b = (BIL_Block*)b;
  if (block_a->request_rank < block_b->request_rank) {
    return -1;
  } else if (block_a->request_rank > block_b->request_rank) {
    return 1;
  } else {
    if (block_a->block_id < block_b->block_id) {
      return -1;
    } else if (block_a->block_id > block_b->block_id) {
      return 1;
    } else {
      return 0;
    }
  }
}

int BIL_Pio_compare_block_group_name(const void* a, const void* b) {
  BIL_Block* block_a = (BIL_Block*)a;
  BIL_Block* block_b = (BIL_Block*)b;
  int file_name_cmp = memcmp(block_a->file_name, block_b->file_name,
                             BIL_MAX_FILE_NAME_SIZE);
  if (file_name_cmp == 0) {
    return memcmp(block_a->var_name, block_b->var_name, BIL_MAX_VAR_NAME_SIZE);
  } else {
    return file_name_cmp;
  }
}

#ifndef DISABLE_PNETCDF
void BIL_Pio_nc_to_mpi_type(nc_type nc_var_type, MPI_Datatype* mpi_var_type_ret,
                            int* var_size_ret) {
  MPI_Datatype mpi_var_type = MPI_BYTE;
  switch(nc_var_type) {
    case NC_BYTE:
      mpi_var_type = MPI_BYTE;
      break;
    case NC_CHAR:
      mpi_var_type = MPI_CHAR;
      break;
    case NC_SHORT:
      mpi_var_type = MPI_SHORT;
      break;
    case NC_INT:
      mpi_var_type = MPI_INT;
      break;
    case NC_FLOAT:
      mpi_var_type = MPI_FLOAT;
      break;
    case NC_DOUBLE:
      mpi_var_type = MPI_DOUBLE;
      break;
    default:
      BIL_Misc_error(BIL_UNSUPPORTED_NC_TYPE_ERR);
  }

  // Set return values.
  MPI_Type_size(mpi_var_type, var_size_ret);
  *mpi_var_type_ret = mpi_var_type;
}
#endif

void BIL_Pio_read_raw_blocks(MPI_Comm all_readers_comm, MPI_Comm io_comm,
                             int num_blocks, BIL_Block* blocks) {
  int i;
  for (i = 0; i < num_blocks; i++) {
    MPI_File fp;
    BIL_Timing_fopen_start(all_readers_comm);
    assert(MPI_File_open(io_comm, blocks[i].file_name, MPI_MODE_RDONLY,
                         BIL->io_hints, &fp) == MPI_SUCCESS);
    BIL_Timing_fopen_stop(all_readers_comm);

    // Get variable and subarray datatype for I/O.
    MPI_Datatype var_type;
    assert(MPI_Type_contiguous(blocks[i].var_size, MPI_BYTE, &var_type)
           == MPI_SUCCESS);
    assert(MPI_Type_commit(&var_type) == MPI_SUCCESS); // added by TP
    MPI_Datatype file_type;
    assert(MPI_Type_create_subarray(blocks[i].num_dims,
                                    blocks[i].file_dim_sizes,
                                    blocks[i].sizes, blocks[i].starts,
                                    MPI_ORDER_C, var_type, 
                                    &file_type) == MPI_SUCCESS);
    assert(MPI_Type_commit(&file_type) == MPI_SUCCESS);

    assert(MPI_File_set_view(fp, BIL->io_header_size, var_type, file_type,
                             (char *)"native",
                             MPI_INFO_NULL) == MPI_SUCCESS);
    // Allocate data and read it collectively.
    blocks[i].data = BIL_Misc_malloc(blocks[i].total_size * blocks[i].var_size);
    BIL_Timing_io_start(all_readers_comm);
    assert(MPI_File_read(fp, blocks[i].data, blocks[i].total_size,
                         var_type, MPI_STATUS_IGNORE) == MPI_SUCCESS);
    BIL_Timing_io_stop(all_readers_comm,
                       blocks[i].total_size * blocks[i].var_size);

    // Clean up.
    MPI_File_close(&fp);
    MPI_Type_free(&var_type);
    MPI_Type_free(&file_type);
  }
}

#ifndef DISABLE_PNETCDF
void BIL_Pio_read_nc_blocks(MPI_Comm all_readers_comm, MPI_Comm io_comm,
                            int num_blocks, BIL_Block* blocks) {
  int i;
  for (i = 0; i < num_blocks; i++) {
    int fp;
    BIL_Timing_fopen_start(all_readers_comm);
    assert(ncmpi_open(io_comm, blocks[i].file_name, NC_NOWRITE,
           BIL->io_hints, &fp) == NC_NOERR);
    BIL_Timing_fopen_stop(all_readers_comm);
  
    ncmpi_begin_indep_data(fp);
  
    // Find the id, type, and size of the variable.
    int var_id;
    assert(ncmpi_inq_varid(fp, blocks[i].var_name, &var_id) == NC_NOERR);
    nc_type var_type;
    assert(ncmpi_inq_vartype(fp, var_id, &var_type) == NC_NOERR);
  
    // Create extra variables specifically for the netCDF API.
    MPI_Offset nc_dim_starts[BIL_MAX_NUM_DIMS];
    MPI_Offset nc_dim_sizes[BIL_MAX_NUM_DIMS];
    int j;
    for (j = 0; j < blocks[i].num_dims; j++) {
      nc_dim_starts[j] = blocks[i].starts[j];
      nc_dim_sizes[j] = blocks[i].sizes[j];
    }
    MPI_Datatype nc_var_type;
    BIL_Pio_nc_to_mpi_type(var_type, &nc_var_type, &(blocks[i].var_size));
    
    // Allocate room for data and read it independently.
    blocks[i].data = BIL_Misc_malloc(blocks[i].total_size * blocks[i].var_size);
    BIL_Timing_io_start(all_readers_comm);
    assert(ncmpi_get_vara(fp, var_id, nc_dim_starts, nc_dim_sizes,
                          blocks[i].data, blocks[i].total_size,
                          nc_var_type) == NC_NOERR);
    BIL_Timing_io_stop(all_readers_comm,
                       blocks[i].total_size * blocks[i].var_size);
    // Clean up.
    ncmpi_end_indep_data(fp);
    ncmpi_close(fp);
  }
}
#endif

void BIL_Pio_extract_block(BIL_Block* sub_block, BIL_Block* block) {
  assert(sub_block->var_size == block->var_size);
  int data_offset = 0, i, j;
  for (i = sub_block->starts[0]; 
       i < sub_block->starts[0] + sub_block->sizes[0]; i++) {
    for (j = sub_block->starts[1];
         j < sub_block->starts[1] + sub_block->sizes[1]; j++) {
      memcpy(sub_block->data + data_offset,
             block->data + 
               ((((i - block->starts[0]) * block->sizes[1] * block->sizes[2])
               + ((j - block->starts[1]) * block->sizes[2])
               + (sub_block->starts[2] - block->starts[2])) * block->var_size),
             sub_block->sizes[2]  * block->var_size);
      data_offset += sub_block->sizes[2] * block->var_size;
    }
  }
  assert(data_offset == sub_block->total_size * sub_block->var_size);
}

void BIL_Pio_insert_block(BIL_Block* sub_block, BIL_Block* block) {
  assert(sub_block->var_size == block->var_size);
  int data_offset = 0, i, j;
  for (i = sub_block->starts[0]; 
       i < sub_block->starts[0] + sub_block->sizes[0]; i++) {
    for (j = sub_block->starts[1];
         j < sub_block->starts[1] + sub_block->sizes[1]; j++) {
      memcpy(block->data + 
               ((((i - block->starts[0]) * block->sizes[1] * block->sizes[2])
               + ((j - block->starts[1]) * block->sizes[2])
               + (sub_block->starts[2] - block->starts[2])) * block->var_size),
             sub_block->data + data_offset,
             sub_block->sizes[2]  * block->var_size);
      data_offset += sub_block->sizes[2] * block->var_size;
    }
  }
  assert(data_offset == sub_block->total_size * sub_block->var_size);
}

void BIL_Pio_exchange_blocks(BIL_Sched_IO* io_sched,
                             BIL_Sched_IO* inv_io_sched) {
  BIL_Timing_comm_start();
  int i;
  // Convenience variables
  int num_recv_blocks = inv_io_sched->num_io_blocks;
  BIL_Block* recv_blocks = inv_io_sched->io_blocks;
  int num_io_blocks = io_sched->num_io_blocks;
  BIL_Block* io_blocks = io_sched->io_blocks;
  
  // Sort your recv blocks by their read rank.
  qsort(recv_blocks, num_recv_blocks, sizeof(BIL_Block),
        BIL_Pio_compare_block_read_rank);
  // Sort the I/O blocks by their group name to perform binary searching.
  qsort(io_blocks, num_io_blocks, sizeof(BIL_Block),
        BIL_Pio_compare_block_group_name);
  qsort(BIL->blocks, BIL->num_blocks, sizeof(BIL_Block), 
        BIL_Pio_compare_block_group_name);

  // In netcdf mode, you will not know the variable size of the blocks you
  // are receiving until now. The variable sizes are stored in the
  // BIL global variable after reading. Fill in these sizes into the inverse
  // I/O schedule.
  for (i = 0; i < num_recv_blocks; i++) {
    BIL_Block* io_block =
      bsearch(&(recv_blocks[i]), BIL->blocks, BIL->num_blocks, 
              sizeof(BIL_Block), BIL_Pio_compare_block_group_name);
    recv_blocks[i].var_size = io_block->var_size;
  }

  int* all_num_send_blocks = 
    BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  memset(all_num_send_blocks, 0, sizeof(int) * BIL->world_size);

  for (i = 0; i < num_recv_blocks; i++) {
    all_num_send_blocks[recv_blocks[i].read_rank]++;
  }
  int* recv_counts = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  for (i = 0; i < BIL->world_size; i++) {
    recv_counts[i] = 1;
  }
  int num_send_blocks = 0;
  MPI_Reduce_scatter(all_num_send_blocks, &num_send_blocks, 
                     recv_counts, MPI_INT, MPI_SUM, BIL->world_comm);

  BIL_Misc_free(all_num_send_blocks);
  BIL_Misc_free(recv_counts);

  // Receive the blocks that you will send.
  BIL_Block* send_blocks =
    BIL_Misc_malloc(sizeof(BIL_Block) * num_send_blocks);
  MPI_Request* recv_requests =
    BIL_Misc_malloc(sizeof(MPI_Request) * num_send_blocks);
  
  for (i = 0; i < num_send_blocks; i++) {
    MPI_Irecv(send_blocks + i, sizeof(BIL_Block), MPI_BYTE, MPI_ANY_SOURCE,
              0, BIL->world_comm, recv_requests + i);
  }
  // Send the blocks that you will receive.
  MPI_Request* send_requests =
    BIL_Misc_malloc(sizeof(MPI_Request) * inv_io_sched->num_io_blocks);
  for (i = 0; i < num_recv_blocks; i++) {
    MPI_Isend(recv_blocks + i, sizeof(BIL_Block), MPI_BYTE,
              recv_blocks[i].read_rank, 0, BIL->world_comm,
              send_requests + i);
  }

  assert(MPI_Waitall(num_send_blocks, recv_requests, MPI_STATUSES_IGNORE)
         == MPI_SUCCESS);
  assert(MPI_Waitall(num_recv_blocks, send_requests,
         MPI_STATUSES_IGNORE) == MPI_SUCCESS);
  // This barrier is need on my mac to get this to work. Not sure why, but
  // it seems that Waitall doesn't block until the isends/irecvs are complete.
  MPI_Barrier(BIL->world_comm);
  BIL_Misc_free(recv_requests);
  BIL_Misc_free(send_requests);

  // Sort your send blocks by their request rank.
  qsort(send_blocks, num_send_blocks, sizeof(BIL_Block),
        BIL_Pio_compare_block_request_rank);

  // Find the total amount of send data.
  int64_t total_send_data = 0;
  for (i = 0; i < num_send_blocks; i++) {
    total_send_data += send_blocks[i].var_size * send_blocks[i].total_size;
  }
  // Allocate a big array for the sending data. Once this array is allocated,
  // the I/O data will be copied to it in such a manner that MPI_Alltoallv
  // can be performed.
  void* send_data = BIL_Misc_malloc(total_send_data);
  // For convenience, set the blocks send data to the appropriate offset
  // in the send_data array.
  if (num_send_blocks > 0) {
    send_blocks[0].data = send_data;
    for (i = 1; i < num_send_blocks; i++) {
      send_blocks[i].data = send_blocks[i - 1].data +
        (send_blocks[i - 1].var_size * send_blocks[i - 1].total_size);
    }
  }

  // Pack each block into send_data.
  for (i = 0; i < num_send_blocks; i++) {
    // Find the corresponding I/O block.
    BIL_Block* io_block =
      bsearch(&(send_blocks[i]), io_blocks, num_io_blocks, sizeof(BIL_Block),
              BIL_Pio_compare_block_group_name);
    assert(io_block != NULL);
    BIL_Pio_extract_block(&(send_blocks[i]), io_block);
  }
  // Free the I/O data now.
  for (i = 0; i < num_io_blocks; i++) {
    BIL_Misc_free(io_blocks[i].data);
  }

  // Allocate a big array for receiving data. This array will be populated
  // during the MPI_Alltoallv call.
  int64_t total_recv_data = 0;
  for (i = 0; i < num_recv_blocks; i++) {
    total_recv_data += recv_blocks[i].var_size * recv_blocks[i].total_size;
  }
  void* recv_data = NULL;
  // Special case - single block I/O where the user has already allocated a
  // buffer.
  if (num_recv_blocks == 1 && BIL->blocks[0].data != NULL) {
    recv_data = BIL->blocks[0].data;
  } else if (num_recv_blocks > 0) {
    recv_data = BIL_Misc_malloc(total_recv_data);
    recv_blocks[0].data = recv_data;
    for (i = 1; i < num_recv_blocks; i++) {
      recv_blocks[i].data = recv_blocks[i - 1].data +
        (recv_blocks[i - 1].var_size * recv_blocks[i - 1].total_size);
    }
  }

  // Find out the recv counts and offsets. Keep in mind that you may receive
  // multiple blocks from the same process.
  recv_counts = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  memset(recv_counts, 0, sizeof(int) * BIL->world_size);
  int* recv_offsets = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  memset(recv_offsets, 0, sizeof(int) * BIL->world_size);
  int cur_recv_offset = 0;
  for (i = 0; i < num_recv_blocks; i++) {
    if (recv_counts[recv_blocks[i].read_rank] == 0) {
      recv_offsets[recv_blocks[i].read_rank] = cur_recv_offset;
    }
    recv_counts[recv_blocks[i].read_rank] +=
      recv_blocks[i].total_size * recv_blocks[i].var_size;
    cur_recv_offset += recv_blocks[i].total_size * recv_blocks[i].var_size;
  }

  // Find out the send counts and offsets. Keep in mind you might send
  // multiple blocks to the same process.
  int* send_counts = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  memset(send_counts, 0, sizeof(int) * BIL->world_size);
  int* send_offsets = BIL_Misc_malloc(sizeof(int) * BIL->world_size);
  memset(send_offsets, 0, sizeof(int) * BIL->world_size);
  int cur_send_offset = 0;
  for (i = 0; i < num_send_blocks; i++) {
    if (send_counts[send_blocks[i].request_rank] == 0) {
      send_offsets[send_blocks[i].request_rank] = cur_send_offset;
    }
    send_counts[send_blocks[i].request_rank] +=
      send_blocks[i].total_size * send_blocks[i].var_size;
    cur_send_offset += send_blocks[i].total_size * send_blocks[i].var_size;
  }
  
  // Exchange blocks.
  MPI_Alltoallv(send_data, send_counts, send_offsets, MPI_BYTE,
                recv_data, recv_counts, recv_offsets, MPI_BYTE,
                BIL->world_comm);

  // Free your send data.
  BIL_Misc_free(send_data);
  // Copy your recv blocks into BIL->blocks. Remember that multiple recv blocks
  // might account for only one block in BIL->blocks.
  qsort(BIL->blocks, BIL->num_blocks, sizeof(BIL_Block), 
        BIL_Pio_compare_block_id);
  if (num_recv_blocks > 1) {
    for (i = 0; i < num_recv_blocks; i++) {
      // Find the corresponding I/O block.
      BIL_Block* io_block = bsearch(&(recv_blocks[i]), BIL->blocks, 
                                    BIL->num_blocks, sizeof(BIL_Block),
                                    BIL_Pio_compare_block_id);
      assert(io_block != NULL);
      if (io_block->data == NULL) {
        io_block->data =
          BIL_Misc_malloc(io_block->var_size * io_block->total_size);
      }
      BIL_Pio_insert_block(&(recv_blocks[i]), io_block);
    }
    // Free the recv_data.
    BIL_Misc_free(recv_data);
  } else if (num_recv_blocks == 1 && BIL->blocks[0].data == NULL) {
    assert(BIL->num_blocks == 1);
    // The reason we are deferencing this void** is because the user passed
    // a NULL buffer that was allocated for them.
    BIL->blocks[0].data = recv_blocks[0].data;
  }

  BIL_Misc_free(send_counts);
  BIL_Misc_free(send_offsets);
  BIL_Misc_free(recv_counts);
  BIL_Misc_free(recv_offsets);
  BIL_Misc_free(send_blocks);

  BIL_Timing_comm_stop();
}

void BIL_Pio_issue(BIL_Sched_IO* io_sched, BIL_Sched_IO* inv_io_sched,
                   int num_groups, int* blocks_to_group_map) {
  int* var_sizes = BIL_Misc_malloc(sizeof(int) * num_groups);
  memset(var_sizes, 0, sizeof(int) * num_groups);
  int s, i;
  int which_block = 0;
  // The I/O may have to be performed in multiple stages. This occurs mostly
  // when the number of files is greater than the number of processes. Depending
  // on how multiple variables are aggregated, this could also occur if many
  // variables are being read in.
  for (s = 0; s < io_sched->num_io_stages; s++) {
    if (BIL->io_type == BIL_RAW) {
      BIL_Pio_read_raw_blocks(BIL->world_comm, MPI_COMM_SELF,
                              1, &(io_sched->io_blocks[which_block]));
      which_block++;
    }
#ifndef DISABLE_PNETCDF
    else if (BIL->io_type == BIL_NC) {
      BIL_Pio_read_nc_blocks(BIL->world_comm, MPI_COMM_SELF,
                             1, &(io_sched->io_blocks[which_block]));
      which_block++;
    }
#endif
    else {
      assert("BIL I/O type not supported" == 0);
    }
  }

#ifndef DISABLE_PNETCDF
  // If we are in netcdf mode, find the size of all the variables read in.
  if (BIL->io_type == BIL_NC) {
    for (i = 0; i < io_sched->num_io_blocks; i++) {
      var_sizes[io_sched->io_blocks[i].io_group] =
        io_sched->io_blocks[i].var_size;
    }
    MPI_Allreduce(MPI_IN_PLACE, var_sizes, num_groups, MPI_INT, MPI_MAX,
                  BIL->world_comm);
    for (i = 0; i < BIL->num_blocks; i++) {
      BIL->blocks[i].var_size = var_sizes[blocks_to_group_map[i]];
    }
  }
#endif
  BIL_Misc_free(var_sizes);
  MPI_Barrier(BIL->world_comm);
  BIL_Pio_exchange_blocks(io_sched, inv_io_sched);
}
