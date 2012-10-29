/*
  Misc utilities used by parallel merge and sample sort functions.
*/
#include <assert.h>
#include <stdlib.h>

/*
  A check for malloc calls.
*/
void* cmalloc(int size) {
  if (size == 0) {
    return NULL;
  }
  void* buf = malloc(size);
  assert(buf != NULL);
  return buf;
}

/*
  A check for free calls.
*/
void cfree(void* buf) {
  if (buf != NULL) {
    free(buf);
  }
}

/*
  Ceil function to avoid including math lib.
*/
int ceil_pos(float val) {
  return (val - (int)val > 0) ? (int)(val + 1) : (int)val;
}

/*
  An extended binary searching algorithm that also gives information regarding
  if an element lands in between a key.
*/
int bsearchx(void *key, void *item_data, int num_items, int item_size,
             int (*compare)(const void *, const void *)) {
  int cmp;
  char* base = (char*)item_data;
  size_t lim;
  char* item;

  /* Check the left edge first. This is beneficial for tree querying */
  if (num_items > 1 && compare(key, item_data + item_size) < 0) {
    int cmp = compare(key, item_data);
    if (cmp < 0) {
      return -1;
    } else if (cmp > 0) {
      return 1;
    } else {
      return 0;
    }
  }

  for (lim = num_items - 1; lim > 0; lim >>= 1) {
    item = base + (lim >> 1) * item_size;
    cmp = compare(key, item);
    if (cmp == 0) {
      return ((item - (char*)item_data) / item_size) * 2;
    }
    if (cmp > 0) {
      if (compare(key, item + item_size) < 0) {
        return (((item - (char*)item_data) / item_size) * 2) + 1;
      }
      /* If key > item then move right */
      base = item + item_size;
      lim--;
    }
  }
  /* Check other edge cases last */
  cmp = compare(key, item_data + (item_size * (num_items - 1)));
  if (cmp > 0) {
    return (num_items * 2) - 1;
  } else if (cmp == 0) {
    return (num_items * 2) - 2;
  } else {
    return -1;
  }
}
