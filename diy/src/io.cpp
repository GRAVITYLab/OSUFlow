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
#include "util.hpp"

//--------------------------------------------------------------------------

extern bool dtype_absolute_address; // addresses in current datatype
                                     // are absolute w.r.t. MPI_BOTTOM
                                     // or relative w.r.t. base address

//--------------------------------------------------------------------------
//
// the following declarations are needed in order to link template functions
// across compilation units
//
template 
void IO::ReadAllData<char>(char*** data, const int *extents, int nb, 
			   Blocking *blocking);
template 
void IO::ReadAllData<unsigned char>(unsigned char*** data, 
				    const int *extents, int nb, 
				    Blocking *blocking);
template 
void IO::ReadAllData<int16_t>(int16_t*** data, const int *extents, int nb, 
			      Blocking *blocking);
template 
void IO::ReadAllData<uint16_t>(uint16_t*** data, const int *extents, 
			       int nb, Blocking *blocking);
template 
void IO::ReadAllData<int32_t>(int32_t*** data, const int *extents, int nb, 
			      Blocking *blocking);
template 
void IO::ReadAllData<uint32_t>(uint32_t*** data, const int *extents, 
			       int nb, Blocking *blocking);
template 
void IO::ReadAllData<float>(float*** data, const int *extents, int nb, 
			    Blocking *blocking);
template 
void IO::ReadAllData<double>(double*** data, const int *extents, int nb, 
			     Blocking *blocking);
template 
void IO::ReadAllData<long double>(long double*** data, const int *extents, 
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
// did: domain id
// dim: number of dimensions
// tb: total number of blocks
// mb: maximum number of blocks per process
// comm: MPI communicator
//
IO::IO(int did, int dim, int tb, int mb, MPI_Comm comm) {

  this->did = did;
  this->dim = dim;
  this->tot_b = tb;
  this->max_b = mb;
  this->comm = comm;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

}
//----------------------------------------------------------------------------
//
// initializes parallel writing of analysis block
// filename: output filename
// compress: whether to apply compression to nonempty writes (false by default)
//
void IO::WriteAnaInit(const char *filename, bool compress) {

  int retval = MPI_File_open(comm, (char *)filename, 
			     MPI_MODE_WRONLY | MPI_MODE_CREATE,
			     MPI_INFO_NULL, &fd_out);
  assert(retval == MPI_SUCCESS);
  MPI_File_set_size(fd_out, 0); // start with an empty file every time
  ofst = 0;
  this->compress = compress;
  block_sizes.clear();

}
//----------------------------------------------------------------------------
//
// initializes parallel reading of analysis block
// filename: input filename
// swap_bytes: whether to swap bytes for endian conversion
//  only applies to reading the headers and footer
//  user must swap bytes manually for datatypes because they are custom
// compress: whether to apply decompression to nonempty reads
// glo_num__blocks: total number of blocks in the global domain (output)
// loc_num_blocks: local number of blocks on this process (output)
//
void IO::ReadAnaInit(const char *filename, bool swap_bytes, bool compress,
		     int& glo_num_blocks, int& loc_num_blocks) {

  int64_t *ftr; // footer, allocated by ReadFooter
  int64_t *all_sizes = NULL; // sizes for all blocks in footer

  int retval = MPI_File_open(comm, (char *)filename, MPI_MODE_RDONLY,
			     MPI_INFO_NULL, &fd_in);
  assert(retval == MPI_SUCCESS);
  ofst = 0;
  this->compress = compress;
  block_sizes.clear();

  // root reads footer and distributes offset to each process
  if (rank == 0) {
    ReadFooter(fd_in, comm, ftr, &tot_b, swap_bytes);
    all_sizes = new int64_t[tot_b];
    for (int i = 0; i < tot_b - 1; i++)
      all_sizes[i] = ftr[i + 1] - ftr[i];
    MPI_Offset file_size;
    MPI_File_get_size(fd_in, &file_size);
    all_sizes[tot_b - 1] = file_size - (tot_b + 1) * sizeof(int64_t);
  }
  MPI_Bcast(&tot_b, 1, MPI_INT, 0, comm);
  // maximum number of blocks per process
  max_b = (tot_b >= groupsize ? tot_b / groupsize : 1);

  // get block starting offsets and block sizes
  // cotiguous block distribution to processes w.r.t order blocks appear
  // in the file (as opposed to round robin order)

  block_starts.resize(max_b);
  block_sizes.resize(max_b);
  MPI_Scatter(ftr, max_b, MPI_LONG_LONG, &block_starts[0], max_b, 
	      MPI_LONG_LONG, 0, comm);
  MPI_Scatter(all_sizes, max_b, MPI_LONG_LONG, &block_sizes[0], max_b, 
	      MPI_LONG_LONG, 0, comm);

  if (rank == 0) {
    delete[] ftr;
    delete[] all_sizes;
  }

  // count and save my local number of blocks read in
  int i;
  for (i = 0; i < max_b && block_starts[i] >= 0; i++)
	 ;
  read_ana_b = i;

  glo_num_blocks = tot_b;
  loc_num_blocks = read_ana_b;

}
//----------------------------------------------------------------------------
//
// reads info from file
//
// fd: file descriptor of open file
// swap_bytes: whether to swap bytes for endian conversion
// comm: MPI communicator
// glo_num__blocks: total number of blocks in the global domain (output)
// loc_num_blocks: local number of blocks on this process (output)
//
void IO::ReadInfo(MPI_File fd, bool swap_bytes, MPI_Comm comm,
		  int& glo_num_blocks, int& loc_num_blocks) {

  int64_t *ftr; // footer, allocated by ReadFooter
  int tot_b; // total number of global blocks
  int max_b; // maximum number of blocks in any process
  int rank, groupsize; // MPI usual

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // root reads footer and distributes offset to each process
  if (rank == 0)
    ReadFooter(fd, comm, ftr, &tot_b, swap_bytes);
  MPI_Bcast(&tot_b, 1, MPI_INT, 0, comm);
  // maximum number of blocks per process
  max_b = (tot_b >= groupsize ? tot_b / groupsize : 1);

  // get block starting offsets
  // cotiguous block distribution to processes w.r.t order blocks appear
  // in the file (as opposed to round robin order)
  MPI_Offset block_starts[max_b]; // block starting offsets for all my blocks
  MPI_Scatter(ftr, max_b, MPI_LONG_LONG, &block_starts[0], max_b, 
	      MPI_LONG_LONG, 0, comm);

  // count and save my local number of blocks read in
  int i;
  for (i = 0; i < max_b && block_starts[i] >= 0; i++)
	 ;
  glo_num_blocks = tot_b;
  loc_num_blocks = i;

  // cleanup
  if (rank == 0)
    delete[] ftr;

}
//----------------------------------------------------------------------------
//
// finalizes parallel writing of analysis block
//
void IO::WriteAnaFinalize() {

  // gather number of output blocks from each process to root
  int *num_blks = NULL;
  int nb = (int)block_sizes.size();
  if (rank == 0)
    num_blks = new int[groupsize];
  MPI_Gather(&nb, 1, MPI_INT, num_blks, 1, MPI_INT, 0, comm);

  int *displs = new int[groupsize];
  int64_t *all_sizes = NULL; // all datatype sizes for everyone's blocks
  int64_t *reorder_sizes = NULL; // reordered sizes
  int tot_out_blks = 0;
  if (rank == 0) {
    for (int i = 0; i < groupsize; i++) {
      tot_out_blks += num_blks[i];
      displs[i] = (i == 0 ? 0 : displs[i - 1] + num_blks[i - 1]);
    }
    all_sizes = new int64_t[tot_out_blks];
    reorder_sizes = new int64_t[tot_out_blks];
  }
  // gather sizes of output blocks from each process to root
  MPI_Gatherv(&block_sizes[0], (int)block_sizes.size(), MPI_LONG_LONG, 
	      all_sizes, num_blks, displs, MPI_LONG_LONG, 0, comm);

  // root writes the footer
  if (rank == 0) {
    ReorderFooter(all_sizes, reorder_sizes, num_blks, tot_out_blks);
    WriteFooter(fd_out, reorder_sizes, tot_out_blks);
    delete[] num_blks;
    delete[] all_sizes;
    delete[] reorder_sizes;
  }

  MPI_File_close(&fd_out);

}
//----------------------------------------------------------------------------
//
// finalizes parallel reading of analysis block
//
void IO::ReadAnaFinalize() {

  MPI_File_close(&fd_in);
  block_starts.clear();
  block_sizes.clear();
  
}
//----------------------------------------------------------------------------
//
// reorders block sizes from independent order to collective order
// indeoendent order is all blocks of a process grouped together
// collective order is the order that processes would write collectively, ie,
//   loop over processes contributing blocks, taking one block from
//   each process at a time
//
// in_blks: sizes of input blocks in independent order
// out_blks: sizes of output blocks reordered in collective order
// num_blks: number of blocks that each process has
// tot_blks: total number of output blocks
//
// caller needs to allocate out_blks
//
void IO::ReorderFooter(int64_t *in_blks, int64_t *out_blks, int *num_blks, 
		       int tot_blks) {

  bool *in_used = new bool[tot_blks]; // this input element was used already
  for (int i = 0; i < tot_blks; i++)
    in_used[i] = false;

  int in_start = 0; // in_blks starting index
  int in_cur = 0; // in_blks current index
  int num_cur = 0; // num_blks current index
  int out_cur = 0; // out_blks current index

  while (out_cur < tot_blks) {

    out_blks[out_cur++] = in_blks[in_cur]; // copy the current element
    in_used[in_cur] = true;
    if (out_cur == tot_blks) // all done
      break;

    while (num_blks[num_cur] == 0) // skip over 0 values in num_blks
      num_cur++;

    if (in_cur + num_blks[num_cur] < tot_blks) // continue marching along
      in_cur += num_blks[num_cur++]; // jump to next value in in_blks
    else { // reset back to start
      in_start++;
      while (in_used[in_start]) // next unused input
	in_start++;
      in_cur = in_start;
      num_cur = 0;
    }

  }

  delete[] in_used;

}
//----------------------------------------------------------------------------
//
// initializes parallel reading of analysis or of data
// filename: output filename
//
void IO::ReadDataInit(const char *filename) {

  int retval = MPI_File_open(comm, (char *)filename, MPI_MODE_RDONLY,
			     MPI_INFO_NULL, &fd_in);
  assert(retval == MPI_SUCCESS);

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
void IO::ReadAllData(T*** data, const int *extents, int nb, 
		     Blocking *blocking) {

  int starts[DIY_MAX_DIM], sizes[DIY_MAX_DIM]; // starts and sizes of a block
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
// swap_bytes: whether to swap data bytes after reading (default = false)
//
// side effects: memory for data will be allocated by this function
//
template<typename T> 
void IO::ReadData(T* &data, const int *starts, 
		  const int *sizes, const int *extents,
		  bool swap_bytes) {

  int si[DIY_MAX_DIM]; // size of entire dataset
  int su[DIY_MAX_DIM]; // size of the block
  int st[DIY_MAX_DIM]; // start of the block
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
    errcode = MPI_File_read_all(fd_in, data, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, (char *)"MPI_File_read_all empty data");
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
    handle_error(errcode, comm, (char *)"MPI_File_read_all nonempty data");
  int count;
  MPI_Get_count(&status, dtype, &count);
  assert(count == tot_size);

  if (swap_bytes && sizeof(data[0] > 1))
    swap((char *)&data[0], tot_size, sizeof(data[0]));

  MPI_Type_free(&filetype);
  MPI_Type_free(&dtype);

}
//----------------------------------------------------------------------------
//
// writes all analysis blocks in parallel with all other procs
// ana: array of pointers to analysis blocks
// nb: number of blocks
// max_nb: maximum number of blocks in any process
// hdrs: headers, one per analysis block (NULL if not used)
// type_func: pointer to function that creates MPI datatype for item 
//
void IO::WriteAllAna(void **ana, int nb, int max_nb, int **hdrs, 
		     void (*type_func)(void*, int, int, MPI_Datatype*)) {

  int s; // temp
  MPI_Offset big_size; // MPI_Offset version of datatype size
  MPI_Offset scan_size; // exclusive scan of sizes before me
  MPI_Offset tot_size; // total size of the iteration
  MPI_Offset unused = 0;
  MPI_Offset unused1 = 0;
  MPI_Datatype dtype; // analysis block datatype
  MPI_Datatype ctype; // combined header + analysis block datatype

  // a naive implementation for now
  for (int i = 0; i < max_nb; i++) {

    MPI_Offset ofst0 = ofst;

    if (i < nb) { // non-null block

      // combine header
      type_func(ana[i], did, i, &dtype); 
      if (hdrs) {
	void *addr; // base address of original datatype
	if (dtype_absolute_address)
	  addr = DIY_BOTTOM;
	else
	  addr = ana[i];
	struct map_block_t map[] = {
	  { MPI_INT, ADDR, DIY_MAX_HDR_ELEMENTS, DIY_Addr(hdrs[i]) },
	  { dtype,   ADDR, 1,                    DIY_Addr(addr)    },
	};
	DIY_Create_struct_datatype(0, 2, map, &ctype);
      }
      else 
	ctype = dtype;

      if (compress) { // compressed
	vector<unsigned char> comp_buf;
	int comp_size;
	CompressBlock(ana[i], ctype, comm, &comp_buf, &comp_size);

	// use MPI_Offset to not overflow the scan
	// in MPI 2.2, should use MPI_OFFSET datatype, but BG/P is MPI 2.0
	// using MPI_LONG_LONG for backwards portability
	big_size = comp_size;
	MPI_Exscan(&big_size, &scan_size, 1, MPI_LONG_LONG, MPI_SUM, comm);

	if (rank > 0)
	  ofst += scan_size;
	MPI_Datatype d = MPI_BYTE;
	WriteDatatype(fd_out, &comp_buf[0], &d, 
		      comp_size, ofst);
      } 
      else { // uncompressed

	MPI_Type_size(ctype, &s);
	big_size = s; 	// use MPI_Offset to not overflow the scan
	// in MPI 2.2, should use MPI_OFFSET datatype, but BG/P is MPI 2.0
	// using MPI_LONG_LONG for backwards portability
	MPI_Exscan(&big_size, &scan_size, 1, MPI_LONG_LONG, MPI_SUM, comm);
	ofst0 = ofst;
	if (rank > 0)
	  ofst += scan_size;
	if (dtype_absolute_address)
	  WriteDatatype(fd_out, 0, &ctype, 1, ofst);
	else
	  WriteDatatype(fd_out, ana[i], &ctype, 1, ofst);

      }
      MPI_Allreduce(&big_size, &tot_size, 1, MPI_LONG_LONG, MPI_SUM, comm);
      ofst = ofst0 + tot_size;
      DIY_Destroy_datatype(&dtype);
      if (hdrs)
	DIY_Destroy_datatype(&ctype);

    }

    else  { // null block
      MPI_Exscan(&unused, &unused, 1, MPI_LONG_LONG, MPI_SUM, comm);
      MPI_Datatype unused_type = MPI_BYTE;
      WriteDatatype(fd_out, 0, &unused_type, 0, ofst);
      MPI_Allreduce(&unused, &unused1, 1, MPI_LONG_LONG, MPI_SUM, comm);
    }

  }

}
//----------------------------------------------------------------------------
//
// reads all analysis blocks in parallel with all other procs
//
// ana: array of pointers to analysis blocks (output)
// hdrs: headers, one per analysis block
//   (allocated by caller, pass NULL if not used)
// create_type_func: pointer to function that takes a block local id, 
//   block header, and creates (allocates) a block and creates a DIY datatype 
//
void IO::ReadAllAna(void** &ana, int **hdrs, 
		    void* (*create_type_func)(int, int, int *, 
					      MPI_Datatype *)) {

  MPI_Datatype dtype;

  ana = new void*[read_ana_b];

  // a naive implementation for now
  for (int i = 0; i < max_b; i++) {

    if (i < read_ana_b) { // non-null block

      if (compress) { // compressed

	// read entire block
	dtype = MPI_BYTE;
	unsigned char *comp_buf = new unsigned char[block_sizes[i]];
	int size = (int)block_sizes[i]; // number of bytes in block
	ReadDatatype(fd_in, comp_buf, &dtype, size, block_starts[i]);

	// decompress block
	vector<unsigned char> decomp_buf;
	int decomp_size; // size of decompressed block in bytes
	DecompressBlockToBuffer(comp_buf, size, &decomp_buf, 
			      &decomp_size);

	// copy header
	memcpy(hdrs[i], &decomp_buf[0], DIY_MAX_HDR_ELEMENTS * sizeof(int));
	if (swap_bytes)
	  swap((char *)hdrs[i], DIY_MAX_HDR_ELEMENTS, sizeof(int));

	// allocate memory and create datatype
	if (hdrs)
	  ana[i] = create_type_func(did, i, hdrs[i], &dtype);
	else
	  ana[i] = create_type_func(did, i, NULL, &dtype);

	// decompress block and fill datatype
	int hdr_size = DIY_MAX_HDR_ELEMENTS * sizeof(int);
	DecompressBlockToDatatype(comp_buf + hdr_size, size - hdr_size,
				  ana[i], dtype, comm);

	delete[] comp_buf;

      }

      else { // uncompressed

	// read header
	if (hdrs)
	  ReadHeader(fd_in, 0, hdrs[i], block_starts[i]);

	// allocate memory and create datatype
	if (hdrs)
	  ana[i] = create_type_func(did, i, hdrs[i], &dtype);
	else
	  ana[i] = create_type_func(did, i, NULL, &dtype);

	// read analysis block
	if (dtype_absolute_address)
	  ReadDatatype(fd_in, MPI_BOTTOM, &dtype, 1, block_starts[i]);
	else
	  ReadDatatype(fd_in, ana[i], &dtype, 1, block_starts[i]);

      }

      DIY_Destroy_datatype(&dtype);

    }

    else { // null block

      if (compress) {

	ReadDatatype(fd_in, NULL, NULL, 0, block_starts[i]);

      }

      else {

	// read null header and null analysis block
	if (hdrs)
	  ReadHeader(fd_in, 1, NULL, block_starts[i]);
	ReadDatatype(fd_in, NULL, NULL, 0, block_starts[i]);

      }

    }

  }

}
//----------------------------------------------------------------------------
//
// independently writes the file footer
//
// fd: open MPI file handle
// blk_sizes: block sizes
// nb_out: total number of output blocks
//
// returns: number of bytes written
//
int IO::WriteFooter(MPI_File fd, const int64_t *blk_sizes, int nb_out) {

  MPI_Status status;
  int errcode;
  int tot_bytes = 0; // total bytes written
  int64_t *footer = new int64_t[nb_out + 1]; // footer ready for writing

  // populate the footer
  footer[0] = 0;
  for (int i = 1; i < nb_out; i++)
    footer[i] = footer[i - 1] + blk_sizes[i - 1];
  footer[nb_out] = nb_out;

  // write the footer to disk
  MPI_File_seek(fd, 0, MPI_SEEK_END);
  errcode = MPI_File_write(fd, footer, nb_out + 1, MPI_LONG_LONG, &status);
  int count;
  MPI_Get_count(&status, MPI_LONG_LONG, &count);
  assert(count == nb_out + 1);
  tot_bytes += (count  * (int)sizeof(int64_t));

  // print the file size
  MPI_Offset ofst;
  MPI_Offset disp;
  MPI_File_get_position(fd, &ofst);
  MPI_File_get_byte_offset(fd, ofst, &disp);
  int size = (int)(disp / 1048576);
  fprintf(stderr, "Output file size is %d MB\n", size);

  delete[] footer;
  return tot_bytes;

}
//----------------------------------------------------------------------------
//
// independently reads the file footer
// footer in file is always ordered by global block id
// output footer is in the same order
// allocates ftr to be the correct size (caller should not allocate it)
//
// fd: open MPI file handle
// comm: MPI communicator
// ftr: footer data (output)
// tb: total number of blocks (output)
// swap_bytes: whether to swap bytes for endian conversion
//
// returns: number of bytes read
//
int IO::ReadFooter(MPI_File fd, MPI_Comm comm, int64_t* &ftr, int *tb, 
		   bool swap_bytes) {

  MPI_Status status;
  int errcode;
  int tot_bytes = 0; // total bytes read
  int count;
  MPI_Offset size;
  int groupsize; // MPI comm size

  MPI_File_get_size(fd, &size);
  MPI_Comm_size(comm, &groupsize);

  // read the total number of blocks
  int64_t tb64; // 64 bit version of tb
  errcode = MPI_File_read_at(fd, size - sizeof(int64_t), &tb64, 1, 
			     MPI_LONG_LONG, &status);
  *tb = (int)tb64;
  MPI_Get_count(&status, MPI_LONG_LONG, &count);
  assert(count == 1);
  if (errcode != MPI_SUCCESS)
    handle_error(errcode, comm,
		 (char *)"MPI_File_read footer number of blocks");
  tot_bytes += (count * (int)sizeof(int64_t));
  if (swap_bytes)
      swap((char *)tb, 1, sizeof(int64_t));

  if (*tb <= 0)
    return 0;

  // allocate footer to min. of # blocks and #procs
  int bf = (*tb >= groupsize ? *tb / groupsize : 1); // nominal blocks per proc
  int pad_b = bf * groupsize; // padded number of footer blocks
  ftr = new int64_t[pad_b];

  // read footer
  errcode = MPI_File_read_at(fd, size - (*tb + 1) * sizeof(int64_t),
			     ftr, *tb, MPI_LONG_LONG, &status);
  if (errcode != MPI_SUCCESS)
    handle_error(errcode, comm, (char *)"MPI_File_read footer block start");
  MPI_Get_count(&status, MPI_LONG_LONG, &count);
  assert(count == *tb);
  tot_bytes += (count * (int)sizeof(int64_t));
  if (swap_bytes)
    swap((char *)ftr, *tb, sizeof(int64_t));

  // pad footer with -1 values for any procs that don't have enough blocks
  for (int i = *tb; i < pad_b; i++)
    ftr[i] = -1;

  return tot_bytes;

}
//----------------------------------------------------------------------------
//
// collectively reads block header
//
// fd: open MPI file handle
// null: 0 = actual read, 1 = empty read
// hdr: header data, allocated by caller
// ofst: current file pointer (bytes), updated by the read
//
// returns: number of bytes read
//
int IO::ReadHeader(MPI_File fd, int null, int *hdr, 
		   MPI_Offset& ofst) {

  MPI_Status status;
  int errcode;
  char b; // unused

  int count = 0;

  if (null) { // empty write
    errcode = MPI_File_read_at_all(fd, ofst, &b, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, (char *)"MPI_File_read_at_all empty header");
  }

  else {
    errcode = MPI_File_read_at_all(fd, ofst, (void *)hdr, DIY_MAX_HDR_ELEMENTS, 
				    MPI_INT, &status);
    MPI_Get_count(&status,MPI_INT,&count);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, (char *)"MPI_File_read_at_all header");
    assert(count == DIY_MAX_HDR_ELEMENTS);
    ofst += count * sizeof(int);
  }

  return(count * (int)sizeof(int));

}
//----------------------------------------------------------------------------
//
// collectively writes an MPI datatype
//
// fd: open MPI file handle
// addr: starting address for source data
// d: pointer to MPI datatype
// num_d: number of instances of d, the datatype (0 = null write)
// ofst: current file pointer (bytes), updated by the write
//
// returns: number of bytes written
//
int IO::WriteDatatype(MPI_File fd, void* addr, 
		      MPI_Datatype *d, int num_d, MPI_Offset& ofst) {

  MPI_Status status;
  int errcode;
  char b; // unused
  int bytes = 0; // number of bytes written

  if (!num_d) { // empty write
    errcode = MPI_File_write_at_all(fd, ofst, &b, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, (char *)"MPI_File_write_all empty datatype");
  }
  else {
    errcode = MPI_File_write_at_all(fd, ofst, addr, num_d, *d, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, 
		   (char *)"MPI_File_write_all nonempty datatype");
    MPI_Get_count(&status, MPI_BYTE, &bytes);
    block_sizes.push_back(bytes);
  }

  ofst += bytes;

  return bytes;

}
//----------------------------------------------------------------------------
//
// collectively reads an MPI datatype
//
// byte swapping not handled here because diy doesn't know anything about
// the user's datatype--up to user to deal with endian byte swapping
//
// fd: open MPI file handle
// addr: starting address for destination data, allocated by caller
// d: pointer to MPI datatype
// num_d: number of instances of d, the datatype (0 = null write)
// ofst: current file pointer (bytes), updated by the read
//
// returns: number of bytes read
//
int IO::ReadDatatype(MPI_File fd, void* addr, 
		     MPI_Datatype *d, int num_d, MPI_Offset& ofst) {

  MPI_Status status;
  int errcode;
  char b; // unused
  int bytes = 0; // number of bytes read

  if (!num_d) { // empty read
    errcode = MPI_File_read_at_all(fd, ofst, &b, 0, MPI_BYTE, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, (char *)"MPI_File_read_all empty datatype");
  }
  else {
    errcode = MPI_File_read_at_all(fd, ofst, addr, num_d, *d, &status);
    if (errcode != MPI_SUCCESS)
      handle_error(errcode, comm, 
		   (char *)"MPI_File_read_all nonempty datatype");
    MPI_Get_count(&status, MPI_BYTE, &bytes);
  }

  ofst += bytes;

  return bytes;

}
//----------------------------------------------------------------------------
//
// MPI error handler
// decodes and prints MPI error messages
//
void IO::handle_error(int errcode, MPI_Comm comm, char *str) {

  char msg[MPI_MAX_ERROR_STRING];
  int resultlen;
  MPI_Error_string(errcode, msg, &resultlen);
  fprintf(stderr, "%s: %s\n", str, msg);
  MPI_Abort(comm, 1);

}
//-----------------------------------------------------------------------
