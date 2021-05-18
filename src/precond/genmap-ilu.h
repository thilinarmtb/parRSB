#ifndef _GENMAP_ILU_H_
#define _GENMAP_ILU_H_

#include <genmap-csr-mat.h>
#include <genmap.h>

struct serial_ilu_data {
  csr_mat M;
  genmap_handle h;
};

int serial_ilu_setup(genmap_handle h, struct comm *c, struct serial_ilu_data *d);

int serial_ilu_apply(GenmapScalar *u, GenmapScalar *rhs, struct serial_ilu_data *d,
              buffer *buf);

int serial_ilu_free(struct serial_ilu_data *d);

#endif
