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

#include "io.hpp"

//--------------------------------------------------------------------------
//
// the following declarations are needed in order to link template functions
// across compilation units
//
template 
void IO::ReadAllData<char>(char*** data, const int64_t *extents, int nb, 
			   Blocking *blocking);
template 
void IO::ReadAllData<unsigned char>(unsigned char*** data, 
				    const int64_t *extents, int nb, 
				    Blocking *blocking);
template 
void IO::ReadAllData<int16_t>(int16_t*** data, const int64_t *extents, int nb, 
			      Blocking *blocking);
template 
void IO::ReadAllData<uint16_t>(uint16_t*** data, const int64_t *extents, 
			       int nb, Blocking *blocking);
template 
void IO::ReadAllData<int32_t>(int32_t*** data, const int64_t *extents, int nb, 
			      Blocking *blocking);
template 
void IO::ReadAllData<uint32_t>(uint32_t*** data, const int64_t *extents, 
			       int nb, Blocking *blocking);
template 
void IO::ReadAllData<float>(float*** data, const int64_t *extents, int nb, 
			    Blocking *blocking);
template 
void IO::ReadAllData<double>(double*** data, const int64_t *extents, int nb, 
			     Blocking *blocking);
template 
void IO::ReadAllData<long double>(long double*** data, const int64_t *extents, 
				  int nb, Blocking *blocking);
//--------------------------------------------------------------------------
//
// the following functions decode the template argument and map it to
// an MPI_Datatype of the appropriate type
//
// Courtesy of Dries Kimpe
// ref: http://www.springerlink.com/content/3txyrj15lad5hucl
//

template <typename T> 
void GetType(MPI_Datatype *d) { // catch-all
  fprintf(stderr, "Unknown typename for MPI_Datatype\n");  
  assert(false);
}

#define MAKEFUNC(a, b) \
  template <> void GetType<a> (MPI_Datatype *d) { *d = b; }

// add supported types here
MAKEFUNC(char, MPI_CHAR);
MAKEFUNC(unsigned char, MPI_UNSIGNED_CHAR);
MAKEFUNC(int16_t, MPI_SHORT);
MAKEFUNC(uint16_t, MPI_UNSIGNED_SHORT);
MAKEFUNC(int32_t, MPI_INT);
MAKEFUNC(uint32_t, MPI_UNSIGNED);
MAKEFUNC(float, MPI_FLOAT);
MAKEFUNC(double, MPI_DOUBLE);
MAKEFUNC(long double, MPI_LONG_DOUBLE);

template <typename T> 
void GetDatatype(MPI_Datatype *p) {
  GetType<T> (p);
}

//----------------------------------------------------------------------------
//
// constructor
//
// dim: number of dimensions
// tb: total number of blocks
// mb: maximum number of blocks per process
// comm: MPI communicator
//
IO::IO(int dim, int tb, int mb, MPI_Comm comm) {

  this->dim = dim;
  this->tot_b = tb;
  this->max_b = mb;
  this->comm = comm;

}
//----------------------------------------------------------------------------
//
// initializes parallel writing of analysis block
// filename: output filename
// write_footer: whether to write a footer with starting locations of blocks
//
void IO::WriteAnaInit(const char *filename, bool write_footer) {

  assert(MPI_File_open(comm, (char *)filename, 
		       MPI_MODE_WRONLY | MPI_MODE_CREATE,
		       MPI_INFO_NULL, &fd_out) == MPI_SUCCESS);
  MPI_File_set_size(fd_out, 0); // start with an empty file every time
  ofst = 0;
  footer = write_footer;

}
//----------------------------------------------------------------------------
//
// finalizes parallel writing of analysis block
//
void IO::WriteAnaFinalize() {

  int rank, groupsize; // MPI usual
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // gather number of output blocks from each process to root
  int *num_blks = NULL;
  unsigned int nb = sizes.size();
  if (rank == 0)
    num_blks = new int[groupsize];
  MPI_Gather(&nb, 1, MPI_UNSIGNED, num_blks, 1, MPI_UNSIGNED, 0, comm);

  int  *displs = new int[groupsize];
  unsigned int *all_sizes = NULL; // all datatype sizes for everyone's blocks
  int tot_out_blks = 0;
  if (rank == 0) {
    for (int i = 0; i < groupsize; i++) {
      tot_out_blks += num_blks[i];
      displs[i] = (i == 0 ? 0 : displs[i - 1] + num_blks[i - 1]);
    }
    all_sizes = new unsigned int[tot_out_blks];
  }
  // gather sizes of output blocks from each process to root
  MPI_Gatherv(&sizes[0], sizes.size(), MPI_UNSIGNED, all_sizes, num_blks, 
	      displs, MPI_UNSIGNED, 0, comm);

  // root writes the footer
  if (rank == 0) {
    WriteFooter(fd_out, all_sizes, tot_out_blks);
    delete[] num_blks;
    delete[] all_sizes;
  }

  MPI_File_close(&fd_out);

}
//----------------------------------------------------------------------------
//
// initializes parallel reading of analysis or of data
// filename: output filename
//
void IO::ReadDataInit(const char *filename) {

  assert(MPI_File_open(comm, (char *)filename, MPI_MODE_RDONLY,
		       MPI_INFO_NULL, &fd_in) == MPI_SUCCESS);

}
//----------------------------------------------------------------------------
//
// finalizes parallel reading of analysis or of data
//
void IO::ReadDataFinalize() {

  MPI_File_close(&fd_in);

}
//----------------------------------------------------------------------------
//
// reads all blocks of data from from a file in parallel with all other procs
//
// data: array of pointers to data blocks (output)
// extents: sizes of the entire data domain
// nb: number of blocks for this process
//
// side effects: memory for data will be allocated
//
template<typename T> 
void IO::ReadAllData(T*** data, const int64_t *extents, int nb, 
		     Blocking *blocking) {

  int64_t starts[MAX_DIM], sizes[MAX_DIM]; // starts and sizes of a block
  int i, j;

  *data = new T*[max_b];

  // a fairly naive implementation
  for (j = 0; j < max_b; j++) { // for max blocks per process

    if (j >= nb) { // null block
      for (i = 0; i < dim; i++) {
	starts[i] = 0; // 0's used to trigger a null read
	sizes[i] = 0;
      }
      ReadData((*data)[j], starts, sizes, extents);
    }
    else { // non-null block
      blocking->BlockStartsSizes(j, starts, sizes);
      ReadData((*data)[j], starts, sizes, extents);
    }

  }

}
//----------------------------------------------------------------------------
//
// reads a block of data from from a file in parallel with all other procs
//
// data: pointer to data (output)
// starts: starting indices of the block
// sizes: sizes of the block
// extents: sizes of the entire data domain
//
// side effects: memory for data will be allocated by this function
//
template<typename T> 
void IO::ReadData(T* &data, const int64_t *starts, 
		  const int64_t *sizes, const int64_t *extents) {

  int si[MAX_DIM]; // size of entire dataset
  int su[MAX_DIM]; // size of the block
  int st[MAX_DIM]; // start of the block
  int errcode;
  int tot_size = 1; // number of data elements in the block
  MPI_Datatype filetype;
  MPI_Datatype dtype;
  MPI_Status status;
  int i;

  for(i = 0; i < dim; i++)
    tot_size *= sizes[i];

  // short-circuit for size 0 reads
  if (!tot_size) {
    for (i = 0; i < dim; i++) { // these need to be set even though not used
      si[dim - 1 - i] = extents[i];
      st[dim - 1 - i] = 0; // any value will do
      su[dim - 1 - i] = 1; // any nonzero value will do
    }
    MPI_Type_create_subarray(dim, si, su, st, MPI_ORDER_C, MPI_BYTE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(fd_in, 0, MPI_BYTE, filetype, (char *)"native",
		      MPI_INFO_NULL);
    errcode = MPI_File_read_all(fd_in, data, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_read_all empty data");
    return;
  }

  // set subarray params
  // reversed orders are intentional: starts and subsizes are [z][y][x]
  // starts, sizes are [x][y][z]
  for (i = 0; i < dim; i++) {
    si[dim - 1 - i] = extents[i];
    st[dim - 1 - i] = starts[i];
    su[dim - 1 - i] = sizes[i];
  }

  data = new T[tot_size];
  GetDatatype<T>(&dtype); // get matching MPI datatype for T

  // do the read
  MPI_Type_create_subarray(dim, si, su, st, MPI_ORDER_C, dtype, &filetype);
  MPI_Type_commit(&filetype);
  MPI_File_set_view(fd_in, 0, MPI_BYTE, filetype, (char *)"native",
		    MPI_INFO_NULL);
  errcode = MPI_File_read_all(fd_in, data, tot_size, dtype, &status);
  if (errcode != MPI_SUCCESS)
    handle_error(errcode, (char *)"MPI_File_read_all nonempty data");
  assert(status.count == (int)(sizeof(T)) * tot_size);

  MPI_Type_free(&filetype);
  MPI_Type_free(&dtype);

}
//----------------------------------------------------------------------------
//
// writes all analysis blocks in parallel with all other procs
// ana: array of pointers to analysis blocks
// nb: number of blocks
// max_nb: maximum number of blocks in any process
// dtype: pointer to MPI datatype of blocks
//
void IO::WriteAllAna(void **ana, int nb, int max_nb, MPI_Datatype *dtype) {

  // a naive implementation for now; worst case: one process has all objects
  for (int i = 0; i < max_nb; i++) {
    if (i < nb) // non-null block
      WriteDatatype(fd_out, 0, (MPI_Aint)ana[i], dtype, ofst);
    else // write null block merge handle
      WriteDatatype(fd_out, 1, 0, dtype, ofst);
  }

}
//----------------------------------------------------------------------------
//
// sends one item to a destination rank
// item: item to be sent
//  specifically, the base address of the displacements in the item's datatype
//  usually the address of the item but could be MPI_BOTTOM if absolute
//  addresses are used in the datatype
// hdr: header containing quantity information
// dest_rank: destination rank
// id: sending identification (used to tag communication)
// dtype: MPI datatype of block
//
void IO::SendItem(char *item, int *hdr, int dest_rank, int id, 
		  MPI_Datatype *dtype) {

  MPI_Request req;
  snd_item si; // sent item

  si.hdr = hdr;
  MPI_Isend(hdr, MAX_HDR_ELEMENTS, MPI_INT, dest_rank, 2 * id, comm, &req);
  si.hreq = req;
  MPI_Isend(item, 1, *dtype, dest_rank, 2 * id + 1, comm, &req);
  si.ireq = req;
  snd_items.push_back(si);

}
//----------------------------------------------------------------------------
//
// initiates receiving one item from a source rank
// item is not available until FinishRecvItems is called
//
// src_rank: source rank
// id: receiving identification (used to tag communication)
//
// side effects: allocates space for the new block
//
// returns: pointer to the item
//
void IO::StartRecvItem(int src_rank, int id) {

  int *hdr = new int[MAX_HDR_ELEMENTS];
  MPI_Request req;
  rcv_hdr rh; // received header

  MPI_Irecv(hdr, MAX_HDR_ELEMENTS, MPI_INT, src_rank, 2 * id, comm, &req);
  rh.hdr = hdr;
  rh.req = req;
  rh.proc = src_rank;
  rh.tag = 2 * id;
  rcv_hdrs.push_back(rh);

}
//----------------------------------------------------------------------------
//
// completes item communication
//
// items: array of pointers to items (output)
// ids: array of receiving identifications (communication tags)
// CreateItem: pointer to user-supplied function that takes a header
//   and creates an item, returning a pointer to it
// RecvItemDtype: pointer to user-supplied function
//   that takes an item and creates an MPI datatype for it
//
// char *'s are used as generic pointers to any bytes
//
// side effects: allocates space for the new items in mh
//
// returns: pointer to the array of item pointers
//
void IO::FinishRecvItems(char **items, int *ids,
			 char * (*CreateItem)(int *),
			 MPI_Datatype* (*RecvItemDtype)(char *)) {

  MPI_Request req;
  vector<MPI_Request> reqs, reqs1;

  // flush completion of header sends
  for (unsigned int i = 0; i < snd_items.size(); i++)
    reqs.push_back(snd_items[i].hreq);
  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);

  // flush completion of header receives
  reqs.clear();
  for (unsigned int i = 0; i < rcv_hdrs.size(); i++)
    reqs.push_back(rcv_hdrs[i].req);
  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);

  // post item receives
  reqs.clear();
  for (unsigned int i = 0; i < rcv_hdrs.size(); i++) {
    items[i] = CreateItem(rcv_hdrs[i].hdr);
    MPI_Datatype *dm = RecvItemDtype(items[i]);
    MPI_Irecv(MPI_BOTTOM, 1, *dm, rcv_hdrs[i].proc, rcv_hdrs[i].tag + 1, 
	      comm, &req);
    MPI_Type_free(dm);
    delete dm;
    reqs1.push_back(req);
    ids[i] = rcv_hdrs[i].tag / 2;
  }
  rcv_hdrs.clear();

  // flush completion of item sends
  reqs.clear();
  for (unsigned int i = 0; i < snd_items.size(); i++)
    reqs.push_back(snd_items[i].ireq);
  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUS_IGNORE);
  snd_items.clear();

  // flush completion of item receives
  MPI_Waitall(reqs1.size(), &reqs1[0], MPI_STATUS_IGNORE);

}
//----------------------------------------------------------------------------
//
// independently writes the file footer
// footer in file is always ordered by global block id
//
// fd: open MPI file handle
// ftr: footer data
// tb: total number of blocks
//
// returns: number of bytes written
//
unsigned int IO::WriteFooter(MPI_File fd, const unsigned int *ftr, 
			     uint32_t nb_out) {

  MPI_Status status;
  int errcode;
  int blk_start = 0; // start of block in the file (value written in footer)
  int tot_bytes = 0; // total bytes written

  // write the block starts
  for (int i = 0; i < (int)nb_out; i++) {
    errcode = MPI_File_write(fd, &blk_start, 1, MPI_UNSIGNED, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_write footer block start");
    assert(status.count == sizeof(unsigned int));
    tot_bytes += status.count;
    blk_start += ftr[i];
  }

  // write the total number of blocks
  errcode = MPI_File_write(fd, &nb_out, 1, MPI_UNSIGNED, &status);
  assert(status.count == sizeof(uint32_t));
  if (errcode != MPI_SUCCESS)
    handle_error(errcode, (char *)"MPI_File_write footer number of blocks");
  tot_bytes += status.count;

  return tot_bytes;

}
//----------------------------------------------------------------------------
//
// independently reads the file footer
// footer in file is always ordered by global block id
// output footer is in the same order
//
// fd: open MPI file handle
// ftr: footer data (output)
// tb: total number of blocks (output)
//
// side effects: allocates ftr
//
// returns: number of bytes read
//
unsigned int IO::ReadFooter(MPI_File fd, unsigned int **ftr, 
			       uint32_t *tb) {

  MPI_Status status;
  int errcode;
  int tot_bytes = 0; // total bytes read
  int i;

  MPI_File_seek(fd, -sizeof(uint32_t), MPI_SEEK_END);

  // read the total number of blocks
  errcode = MPI_File_read(fd, tb, 1, MPI_UNSIGNED, &status);
  assert(status.count == sizeof(uint32_t));
  if (errcode != MPI_SUCCESS)
    handle_error(errcode, (char *)"MPI_File_read footer number of blocks");
  tot_bytes += status.count;

  // allocate footer
  if (*tb > 0) {
    assert((*ftr = new unsigned int[*tb]) != NULL);
    // assumes that sizeof int and unsigned int does not change from 
    // writing to reading machines
    MPI_File_seek(fd, -*tb * sizeof(unsigned int) - sizeof(int), MPI_SEEK_END);

    // read the block starts
    for (i = 0; i < (int)*tb; i++) {
      errcode = MPI_File_read(fd, &(*ftr)[i], 1, MPI_UNSIGNED, &status);
      if (errcode != MPI_SUCCESS)
	handle_error(errcode, (char *)"MPI_File_read footer block start");
      assert(status.count == sizeof(unsigned int));
      tot_bytes += status.count;
    }
  }

  return tot_bytes;

}
//----------------------------------------------------------------------------
//
// collectively writes an MPI datatype
//
// fd: open MPI file handle
// null: 0 = actual write, 1 = empty write
// addr: starting address for source data
// d: pointer to MPI datatype
// ofst: current file pointer (bytes), updated by the write
//
// returns: number of bytes written
//
unsigned int IO::WriteDatatype(MPI_File fd, int null, MPI_Aint addr, 
			       MPI_Datatype *d, MPI_Offset& ofst) {

  MPI_Status status;
  int errcode;
  char b; // unused

  MPI_File_set_view(fd, ofst, MPI_BYTE, MPI_BYTE, (char *)"native", 
		    MPI_INFO_NULL);
  if (null) { // empty write
    errcode = MPI_File_write_all(fd, &b, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_write_all empty datatype");
  }
  else {
    errcode = MPI_File_write_all(fd, (void *)addr, 1, *d, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_write_all nonempty datatype");
  }

  ofst += status.count;
  return status.count;

}
//----------------------------------------------------------------------------
//
// collectively reads an MPI datatype
//
// fd: open MPI file handle
// null: 0 = actual read, 1 = empty read
// addr: starting address for destination data, allocated by caller
// d: pointer to MPI datatype
// ofst: current file pointer (bytes), updated by the read
//
// returns: number of bytes read
//
unsigned int IO::ReadDatatype(MPI_File fd, int null, MPI_Aint addr, 
			      MPI_Datatype *d, MPI_Offset& ofst) {

  MPI_Status status;
  int errcode;
  char b; // unused

  MPI_File_set_view(fd, ofst, MPI_BYTE, MPI_BYTE, (char *)"native", 
		    MPI_INFO_NULL);
  if (null) { // empty read
    errcode = MPI_File_read_all(fd, &b, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_read_all empty datatype");
  }
  else {
    errcode = MPI_File_read_all(fd, (void *)addr, 1, *d, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, (char *)"MPI_File_read_all nonempty datatype");
  }

  ofst += status.count;
  return status.count;

}
//----------------------------------------------------------------------------
//
// MPI error handler
// decodes and prints MPI error messages
//
void IO::handle_error(int errcode, char *str) {

  char msg[MPI_MAX_ERROR_STRING];
  int resultlen;
  MPI_Error_string(errcode, msg, &resultlen);
  fprintf(stderr, "%s: %s\n", str, msg);
  MPI_Abort(comm, 1);

}
//-----------------------------------------------------------------------
