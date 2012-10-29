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
#include "diy.h"
#include "blocking.hpp"
#include "mpi.h"

struct snd_item { // send-item message (includes both header and item)
  int *hdr; // header
  MPI_Request hreq; // header send request
  MPI_Request ireq; // item send request
};

struct rcv_hdr { // receive-header message
  int *hdr; // header
  MPI_Request req; // header receive request
  int tag; // message tag
  int proc; // source process rank
};

class IO {

 public:

  IO(int dim, int tb, int mb, MPI_Comm comm);
  ~IO(){};
  void WriteAnaInit(const char *filename, bool write_footer);
  void WriteAnaFinalize();
  void ReadDataInit(const char *filename);
  void ReadDataFinalize();
  template<typename T> 
  void ReadAllData(T*** data, const int64_t *extents, int nb, 
		   Blocking *blocking);
  void WriteAllAna(void **ana, int nb, int max_nb, MPI_Datatype *dtype);

  void SendItem(char *item, int *hdr, int dest_rank, int id,
		MPI_Datatype *dtype);
  void StartRecvItem(int src_rank, int id);
  void FinishRecvItems(char **items, int *ids,
		       char * (*CreateItem)(int *),
		       MPI_Datatype* (*RecvItemDtype)(char *));
  unsigned int ReadFooter(MPI_File fd, unsigned int **ftr, uint32_t *tb);

 private:

  void handle_error(int errcode, char *str);
  unsigned WriteFooter(MPI_File fd, const unsigned int *ftr, uint32_t nb_out);
  unsigned int WriteDatatype(MPI_File fd, int null, MPI_Aint addr,
			     MPI_Datatype *d, MPI_Offset& ofst);
  unsigned int ReadDatatype(MPI_File fd, int null, MPI_Aint addr, 
			    MPI_Datatype *d, MPI_Offset& ofst);
  template<typename T> 
  void ReadData(T*& data, const int64_t *starts, 
		const int64_t *sizes, const int64_t *extents);

  int dim; // number of dimensions in the dataset
  uint32_t tot_b; // total number of blocks
  int max_b; // maximum number of blocks per process
  MPI_File fd_in, fd_out; // input, output file descriptor
  MPI_Offset ofst; // current file pointer
  vector <size_t> sizes; // datatype sizes for all my blocks
  bool footer; // whether to write a footer
  MPI_Comm comm; // communicator

  vector<snd_item> snd_items; // posted send-MSC messages
  vector<rcv_hdr> rcv_hdrs; // posted receive-header messages

};

#endif
