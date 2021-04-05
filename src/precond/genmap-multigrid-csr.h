#ifndef _GENMAP_MULTIGRID_CSR_H_
#define _GENMAP_MULTIGRID_CSR_H_

#include <genmap-multigrid-csr-mat.h>
#include <genmap.h>

typedef struct mg_level_csr *mgLevel;

struct mg_level_csr {
  struct mg_data_csr *data;
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *J; // Interpolation/Restriction operator
  csr_mat M;         // Laplacian for the level
};

struct mg_data_csr {
  struct comm c;
  genmap_handle h;
  struct gs_data *top;
  buffer bfr;
  int nlevels;
  mgLevel *levels;
  uint *level_off;
  GenmapScalar *y, *x, *b, *u, *rhs, *buf;
};

typedef struct {
  ulong r, c;
  uint proc;
} csr_entry;

typedef struct {
  ulong r, c, rn, cn;
  uint p;
  GenmapScalar v;
} entry;

void mg_setup_csr(genmap_handle h, struct comm *c, struct mg_data_csr *d);

void mg_free_csr(struct mg_data_csr *d);

#endif
