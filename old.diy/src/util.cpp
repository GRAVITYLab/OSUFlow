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

#include "util.hpp"
#include <stdio.h>

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
// swap8(n)
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
// swap4(n)
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
// swap2(n)
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

