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
  ulong *row_id;

  GenmapScalar *diag;

  ulong *col;
  GenmapScalar *v;

  struct gs_data *gsh;
};

typedef struct csr_mat_ *csr_mat;

void csr_mat_setup(struct csr_mat_ **M, struct array *entries, struct comm *c,
                   buffer *bfr);

void csr_mat_apply(GenmapScalar *y, csr_mat M, GenmapScalar *x, buffer *buf);

void csr_mat_print(csr_mat M, struct comm *c);

int csr_mat_get_local(double *val, uint *off, csr_mat M, uint i, ulong j);

int csr_mat_get_global(double *val, uint *off, csr_mat M, ulong i, ulong j);

int csr_mat_copy(csr_mat D, csr_mat S);

int csr_mat_free(csr_mat M);

#endif
