#ifndef _GENMAP_CSR_MAT_H_
#define _GENMAP_CSR_MAT_H_

#include <genmap.h>

typedef struct {
  ulong r, c;
  uint proc;
  GenmapScalar v;
} csr_entry;

typedef struct {
  ulong r, c, rn, cn;
  uint p;
  GenmapScalar v;
} entry;

struct csr_mat_ {
  uint rn;
  ulong row_start;
  uint *row_off;

  ulong *col;
  GenmapScalar *v, *diag;

  struct gs_data *gsh;
};

typedef struct csr_mat_ *csr_mat;

void csr_mat_setup(csr_mat *M, struct array *entries, struct comm *c,
                   buffer *bfr);

void csr_mat_apply(GenmapScalar *y, csr_mat M, GenmapScalar *x, buffer *buf);

void csr_mat_print(csr_mat M, struct comm *c);

int csr_mat_get(double *val, csr_mat M, uint i, uint j);

int csr_mat_copy(csr_mat D, csr_mat S);

int csr_mat_free(csr_mat M);

#endif
