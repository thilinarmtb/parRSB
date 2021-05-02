#ifndef _GENMAP_SORT_IMPL_H_
#define _GENMAP_SORT_IMPL_H_

#include "genmap-sort.h"
#ifndef _GENMAP_SORT_H_
#error "genmap-sort.h is not included."
#endif

double get_scalar(struct array *a, uint i, uint off, uint usize, gs_dom type);
void get_extrema(void *extrema, struct sort *s, uint field, struct comm *c);

int set_dest(uint *proc, uint np, ulong start, uint size, ulong nelem);

int load_balance(struct array *a, size_t size, struct comm *c,
                 struct crystal *cr);

int sort_local(struct sort *s);

#endif
