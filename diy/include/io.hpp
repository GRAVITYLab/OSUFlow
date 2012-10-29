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

  IO(int dim, int tb, int mb, MPI_Comm comm);
  ~IO(){};
  void WriteAnaInit(const char *filename, bool compress = false);
  void ReadAnaInit(const char *filename, bool swap_bytes, 
		   bool compress = false);
  void WriteAnaFinalize();
  void ReadAnaFinalize();
  void ReadDataInit(const char *filename);
  void ReadDataFinalize();
  template<typename T> 
  void ReadAllData(T*** data, const int64_t *extents, int nb, 
		   Blocking *blocking);
  void WriteAllAna(void **ana, int nb, int max_nb, int **hdrs, 
		   int num_hdr_elems,
		   void* (*type_func)(void*, int, MPI_Datatype*));
  int ReadAllAna(void** &ana, int **hdrs, 
		 void* (*create_type_func)(int, int *, MPI_Datatype *));

 private:

  void handle_error(int errcode, char *str);
  int WriteFooter(MPI_File fd, const int64_t *blk_sizes, int nb_out);
  int WriteDatatype(MPI_File fd, void* addr,
		    MPI_Datatype *d, int num_d, MPI_Offset& ofst);
  int ReadDatatype(MPI_File fd, void* addr, 
		   MPI_Datatype *d, int num_d, MPI_Offset& ofst);
  template<typename T> 
  void ReadData(T*& data, const int64_t *starts, 
		const int64_t *sizes, const int64_t *extents,
		bool swap_bytes = false);
  int ReadHeader(MPI_File fd, int null, int *hdr, 
		 MPI_Offset& ofst);
  int ReadFooter(MPI_File fd, int64_t* &ftr, int *tb,
		 bool swap_bytes = false);
  void ReorderFooter(int64_t *in_blks, int64_t *out_blks, 
		     int *num_blks, int tot_blks);

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
