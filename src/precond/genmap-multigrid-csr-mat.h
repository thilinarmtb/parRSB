#ifndef _GENMAP_MULTIGRID_CSR_MAT_H_
#define _GENMAP_MULTIGRID_CSR_MAT_H_

#include <genmap.h>

typedef struct csr_mat_ *csr_mat;

struct csr_mat_ {
  uint rn;
  ulong row_start;

  uint *row_off;
  ulong *col;
  GenmapScalar *v, *diag;

  struct gs_data *gsh;
};

void csr_mat_setup(struct array *entries, struct comm *c, csr_mat *M);
void csr_mat_apply(GenmapScalar *y, csr_mat M, GenmapScalar *x);
void csr_mat_print(csr_mat M, struct comm *c);
int csr_mat_free(csr_mat M);
void csr_mat_gather(csr_mat M, struct gs_data *gsh, GenmapScalar *x,
                    GenmapScalar *buf, buffer *bfr);
struct gs_data *get_csr_top(csr_mat M, struct comm *c);

#endif
