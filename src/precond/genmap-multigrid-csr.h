#ifndef _GENMAP_MULTIGRID_CSR_H_
#define _GENMAP_MULTIGRID_CSR_H_

#include <genmap-multigrid-csr-mat.h>
#include <genmap-multigrid.h>
#include <genmap.h>

typedef struct mg_level_csr *mgLevel;
struct mg_data_csr;

struct mg_level_csr {
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // Interpolation/Restriction operator
  csr_mat M;         // Laplacian for the level
};

void mg_level_setup(struct mg_data_csr *d, slong *wrk, buffer *buf);

int mg_level_get_nsmooth_csr(struct mg_level *lvl);

GenmapScalar mg_level_get_sigma_csr(struct mg_level *lvl);

struct mg_data_csr {
  struct comm c;
  int nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y, *x, *b, *u, *rhs, *buf;
};

void mg_setup_csr(genmap_handle h, struct comm *c, struct mg_data *d);

int mg_get_nlevels_csr(struct mg_data *d);

uint *mg_get_level_off_csr(struct mg_data *d);

void mg_free_csr(struct mg_data *d);

#endif
