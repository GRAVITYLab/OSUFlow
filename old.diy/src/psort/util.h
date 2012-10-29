#ifndef UTIL_H
#define UTIL_H 1

#ifdef __cplusplus
extern "C"
{
#endif

void* cmalloc(int size);

void cfree(void* buf);

int ceil_pos(float val);

int bsearchx(void *key, void *item_data, int num_items, int item_size,
             int (*compare)(const void *, const void *));

#ifdef __cplusplus
}
#endif

#endif
