//------------------------------------------------------------------------------
//
// parallel io class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2011 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------

#ifndef _IO
#define _IO

#include <stddef.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include "diy.h"
#include "blocking.hpp"
#include "mpi.h"

using namespace std;

class IO {

 public:

  IO(int did, int dim, int tb, int mb, MPI_Comm comm);
  ~IO(){};
  void WriteAnaInit(const char *filename, bool compress = false);
  void ReadAnaInit(const char *filename, bool swap_bytes, bool compress,
		   int& glo_num_blocks, int& loc_num_blocks);
  void WriteAnaFinalize();
  void ReadAnaFinalize();
  void ReadDataInit(const char *filename);
  void ReadDataFinalize();
  template<typename T> 
  void ReadAllData(T*** data, const int *extents, int nb, 
		   Blocking *blocking);
  void WriteAllAna(void **ana, int nb, int max_nb, int **hdrs, 
		   void (*type_func)(void*, int, int, MPI_Datatype*));
  void ReadAllAna(void** &ana, int **hdrs, 
		  void* (*create_type_func)(int, int, int *, 
					    MPI_Datatype *));
  static void ReadInfo(MPI_File fd, bool swap_bytes, MPI_Comm comm,
		       int& glo_num_blocks, int& loc_num_blocks);

 private:

  int WriteFooter(MPI_File fd, const int64_t *blk_sizes, int nb_out);
  int WriteDatatype(MPI_File fd, void* addr,
		    MPI_Datatype *d, int num_d, MPI_Offset& ofst);
  int ReadDatatype(MPI_File fd, void* addr, 
		   MPI_Datatype *d, int num_d, MPI_Offset& ofst);
  template<typename T> 
  void ReadData(T*& data, const int *starts, 
		const int *sizes, const int *extents,
		bool swap_bytes = false);
  int ReadHeader(MPI_File fd, int null, int *hdr, 
		 MPI_Offset& ofst);
  void ReorderFooter(int64_t *in_blks, int64_t *out_blks, 
		     int *num_blks, int tot_blks);
  static int ReadFooter(MPI_File fd, MPI_Comm comm, int64_t* &ftr, int *tb,
			bool swap_bytes = false);
  static void handle_error(int errcode, MPI_Comm comm, char *str);

  int did; // domain id
  int dim; // number of dimensions in the dataset
  int tot_b; // total number of blocks
  int max_b; // maximum number of blocks per process
  int read_ana_b; // local number of analysis blocks read in
  vector <MPI_Offset> block_starts; // block starting offsets for all my blocks
  vector <int64_t> block_sizes; // datatype sizes for all my blocks
  MPI_File fd_in, fd_out; // input, output file descriptor
  MPI_Offset ofst; // current file pointer
  MPI_Comm comm; // communicator
  int rank, groupsize; // MPI usual
  bool swap_bytes; // whether to swap bytes for endian conversion
  bool compress; // whether to compress output

};

#endif
