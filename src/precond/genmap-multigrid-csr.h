#ifndef _GENMAP_MULTIGRID_CSR_H_
#define _GENMAP_MULTIGRID_CSR_H_

#include <genmap-csr-mat.h>
#include <genmap-multigrid.h>

struct mg_data_csr {
  struct comm *c;

  int nlevels;
  uint *level_off;

  int *nsmooth;
  GenmapScalar *sigma;

  struct gs_data **J;

  GenmapScalar *buf;

  csr_mat *M;
};

void mg_setup_aux_csr(struct mg_data_csr *d, slong *wrk, buffer *buf);

#endif
