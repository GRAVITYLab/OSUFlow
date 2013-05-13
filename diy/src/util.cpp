//---------------------------------------------------------------------------
//
// utilities
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
//----------------------------------------------------------------------------

#include "diy.h"
#include "util.hpp"
#include <stdio.h>
#include <assert.h>

//-----------------------------------------------------------------------
//
// swaps bytes
//
// n: address of items
// nitems: number of items
// item_size: either 2, 4, or 8 bytes
// returns quietly if item_size is 1
//
void swap(char *n, int nitems, int item_size) {

  int i;

  switch(item_size) {
  case 1:
    return;
    break;
  case 2:
    for (i = 0; i < nitems; i++) {
      swap2(n);
      n += 2;
    }
    break;
  case 4:
    for (i = 0; i < nitems; i++) {
      swap4(n);
      n += 4;
    }
    break;
  case 8:
    for (i = 0; i < nitems; i++) {
      swap8(n);
      n += 8;
    }
    break;
  default:
    fprintf(stderr, "Error: size of data must be either 1, 2, 4, or 8 bytes per item\n");

  }

}
//-----------------------------------------------------------------------
//
// Swaps 8  bytes from 1-2-3-4-5-6-7-8 to 8-7-6-5-4-3-2-1 order.
// cast the input as a char and use on any 8 byte variable
//
void swap8(char *n) {

  char *n1;
  char c;

  n1 = n + 7;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//-----------------------------------------------------------------------------
//
// Swaps 4 bytes from 1-2-3-4 to 4-3-2-1 order.
// cast the input as a char and use on any 4 byte variable
//
void swap4(char *n) {

  char *n1;
  char c;

  n1 = n + 3;
  c = *n;
  *n = *n1;
  *n1 = c;

  n++;
  n1--;
  c = *n;
  *n = *n1;
  *n1 = c;

}
//----------------------------------------------------------------------------
//
// Swaps 2 bytes from 1-2 to 2-1 order.
// cast the input as a char and use on any 2 byte variable
//
void swap2(char *n){

  char c;

  c = *n;
  *n = n[1];
  n[1] = c;

}
//----------------------------------------------------------------------------
//
// creates an MPI struct datatype from a typemap
// use the datatype with the address MPI_BOTTOM
//
// addr: base address of the start of the typemap
// map: typemap
// type: MPI datatype (output)
//
void CreateDtype(MPI_Aint addr, vector<map_block_t> *map, 
		 MPI_Datatype *type) {

  // skip map blocks with count 0
  for (int i = 0; i < (int)map->size(); i++) {
    if ((*map)[i].count <= 0) {
      map->erase(map->begin() + i);
      i--; // next time, test the new element that got moved up to this slot
    }
  }

  // form the required vectors
  vector <MPI_Aint> addrs(map->size(), addr);
  vector <int> counts(map->size(), 0);
  vector <MPI_Datatype> base_types(map->size(), 0);
  for(int i = 0; i < (int)map->size(); i++) {
    if ((*map)[i].disp_type == OFST)
      addrs[i] = addrs[i] + (*map)[i].disp;
    else
      addrs[i] = (*map)[i].disp;
    counts[i] = (*map)[i].count;
    base_types[i] = (*map)[i].base_type;
  }

  // create the datatype
  MPI_Type_create_struct((int)map->size(), &counts[0], 
			 (MPI_Aint *)&addrs[0], &base_types[0], type);

}
//----------------------------------------------------------------------------
//
// block compression
//
// addr: address of start of datatype
// dtype: MPI datatype
// comm: MPI communicator
// comp_buf: pointer to compressed buffer, MPI datatype MPI_BYTE (output)
// comp_size: size of compressed buffer in bytes
//
// side effects: grows comp_buf, user's responsibility to free it
//
void CompressBlock(void* addr, MPI_Datatype dtype, MPI_Comm comm,
		   vector<unsigned char> *comp_buf, int *comp_size) {

#ifdef ZLIB

  z_stream strm; // structure used to pass info t/from zlib
  unsigned char *inbuf; // input buffer for deflate (packed by MPI_Pack)
  unsigned char temp_buf[CHUNK]; // temporary
  int num_out; // number of compressed bytes in one deflation
  int ret, flush;

  // pack datatype
  int size; // packed size in bytes
  MPI_Pack_size(1, dtype, comm, &size);
  inbuf = new unsigned char[size];
  memset(inbuf, 0, size); // just to be sure all bytes are initialized
  int pos = 0; // position in inbuf
  MPI_Pack((void *)addr, 1, dtype, inbuf, size, &pos, comm);

  // init deflate state
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);
  assert(ret == Z_OK);

  *comp_size = 0;
  int num_in_chunks = (int)(ceil((float)size / CHUNK)); // number of input chunks

  // read the character buffer in chunks
  ret = Z_STREAM_ERROR; // ret will be reset in loop below, asserted aftertwards
  for (int i = 0; i < num_in_chunks; i++) { 
    strm.next_in = inbuf + i * CHUNK;
    if (i < num_in_chunks - 1) {
      strm.avail_in = CHUNK;
      flush = Z_NO_FLUSH;
    }
    else {
      strm.avail_in = (size % CHUNK ? size - i * CHUNK : CHUNK);
      flush = Z_FINISH;
    }
    do { // compresss each chunk
      strm.avail_out = CHUNK;
      strm.next_out = temp_buf;
      ret = deflate(&strm, flush);
      assert(ret != Z_STREAM_ERROR);
      num_out = CHUNK - strm.avail_out;
      *comp_size += num_out;
      comp_buf->insert(comp_buf->end(), temp_buf, temp_buf + num_out);
    } while (strm.avail_out == 0);
    assert(strm.avail_in == 0); //confirm all input used
  }
  assert(ret == Z_STREAM_END); // final check that everything completed ok
  deflateEnd(&strm);

  delete[] inbuf;

#else

  addr = addr; // quiet oompiler warnings
  dtype = dtype;
  comm = comm;
  comp_buf = comp_buf; 
  comp_size = comp_size;

  fprintf(stderr, "Error: CompressBlock() attempting to do compression "
	  "without zlib installed\n");

#endif

}
//-----------------------------------------------------------------------
//
// block decompression with the result in a datatype
//
// in_buf: pointer to compressed buffer
// in_size: size of compressed buffer (bytes)
// addr: address of start of datatype
// dtype: MPI datatype
// comm: MPI communicator
//
// side effects: fills datatype with decompressed result
//
void DecompressBlockToDatatype(unsigned char *in_buf, int in_size,
			       void* addr, MPI_Datatype dtype, 
			       MPI_Comm comm) {

#ifdef ZLIB

  z_stream strm; // structure used to pass info t/from zlib
  unsigned char temp_buf[CHUNK]; // temporary
  int num_out; // number of compressed bytes in one deflation
  int ret; // return value
  vector<unsigned char> *decomp_buf;
  int decomp_size;

  // allocate inflate state
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit(&strm);
  assert(ret == Z_OK);

  decomp_size = 0;

  // number of input chunks
  int num_in_chunks = (int)(ceil((float)in_size / CHUNK)); 

  // read the character buffer in chunks
  for (int i = 0; i < num_in_chunks; i++) { 
    strm.next_in = in_buf + i * CHUNK;
    if (i < num_in_chunks - 1) 
      strm.avail_in = CHUNK;
    else
      strm.avail_in = (in_size % CHUNK ? in_size - i * CHUNK : CHUNK);
    do { // decompress each chunk
      strm.avail_out = CHUNK;
      strm.next_out = temp_buf;
      ret = inflate(&strm, Z_NO_FLUSH);
      assert(ret != Z_STREAM_ERROR);
      num_out = CHUNK - strm.avail_out;
      decomp_size += num_out;
      decomp_buf->insert(decomp_buf->end(), temp_buf, temp_buf + num_out);
    } while (strm.avail_out == 0);
    if (ret == Z_STREAM_END) // in case decompression finishes early
      break;
  }
  inflateEnd(&strm);

  // fill datatype
  MPI_Unpack(&decomp_buf[0], decomp_size, 0, addr, 1, dtype, comm);

#else

  in_buf = in_buf; // quiet compiler warnings
  in_size = in_size;
  addr = addr;
  dtype = dtype;
  comm = comm;

  fprintf(stderr, "Error: DeCompressBlock() attempting to de decompression "
	  "without zlib installed\n");

#endif

}
//-----------------------------------------------------------------------
//
// block decompression with the result in a contiguous buffer
//
// in_buf: input block buffer (MPI_BYTE datatype)
// in_size: input size in bytes
// decomp_buf: decompressed buffer
// decomp_size: decompressed size in bytes (output)
//
// side effects: grows decomp_buf, user's responsibility to free it
//
void DecompressBlockToBuffer(unsigned char* in_buf, int in_size, 
			     vector<unsigned char> *decomp_buf, 
			     int *decomp_size) {

#ifdef ZLIB

  z_stream strm; // structure used to pass info t/from zlib
  unsigned char temp_buf[CHUNK]; // temporary
  int num_out; // number of compressed bytes in one deflation
  int ret; // return value

  // allocate inflate state
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit(&strm);
  assert(ret == Z_OK);

  *decomp_size = 0;
  int num_in_chunks = (int)(ceil((float)in_size / CHUNK)); // number of input chunks

  // read the character buffer in chunks
  for (int i = 0; i < num_in_chunks; i++) { 
    strm.next_in = in_buf + i * CHUNK;
    if (i < num_in_chunks - 1) 
      strm.avail_in = CHUNK;
    else
      strm.avail_in = (in_size % CHUNK ? in_size - i * CHUNK : CHUNK);
    do { // decompress each chunk
      strm.avail_out = CHUNK;
      strm.next_out = temp_buf;
      ret = inflate(&strm, Z_NO_FLUSH);
      assert(ret != Z_STREAM_ERROR);
      num_out = CHUNK - strm.avail_out;
      *decomp_size += num_out;
      decomp_buf->insert(decomp_buf->end(), temp_buf, temp_buf + num_out);
    } while (strm.avail_out == 0);
    if (ret == Z_STREAM_END) // in case decompression finishes early
      break;
  }
  inflateEnd(&strm);

#else

  in_buf = in_buf; // quiet compiler warnings
  in_size = in_size;
  decomp_buf = decomp_buf;
  decomp_size = decomp_size;

  fprintf(stderr, "Error: DeCompressBlock() attempting to de decompression "
	  "without zlib installed\n");

#endif

}
//-----------------------------------------------------------------------
