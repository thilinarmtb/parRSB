#ifndef _GENMAP_ILU_H_
#define _GENMAP_ILU_H_

#include <genmap-csr-mat.h>
#include <genmap.h>

struct ilu_data {
  csr_mat M;
  genmap_handle h;
};

int ilu_setup(genmap_handle h, struct comm *c, struct ilu_data *d);

int ilu_apply(GenmapScalar *u, GenmapScalar *rhs, struct ilu_data *d,
              buffer *buf);

int ilu_free(struct ilu_data *d);

#endif
