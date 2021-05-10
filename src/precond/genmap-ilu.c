#include <genmap-ilu.h>
#include <genmap-impl.h>

static int ilu0(csr_mat M) {
  const uint rn = M->rn;
  assert(rn > 1);

  const uint *off = M->row_off;
  double *v = M->v;
  const double *diag = M->diag;
  const ulong *col = M->col;

  double a_ii, a_ji, a_ik;

  uint i, ii, j, jj, k, kk;
  for (ii = 0; ii < rn - 1; ii++) { /* Go over number of rows */
    i = ii + 1;
    assert(csr_mat_get(&a_ii, M, i, i) == 0);

    for (jj = ii + 1; jj < rn; jj++) {
      j = jj + 1;
      assert(csr_mat_get(&a_ji, M, j, i) == 0);
      if (fabs(a_ji) < 1e-12)
        continue;

      kk = off[jj];
      while (col[kk] < i)
        kk++;
      /* a_ji = a_ji / a_ii */
      if (col[kk] != i) {
        printf("kk = %u, col[kk] = %lu, j = %u, i = %u a_ji = %lf\n", kk,
               col[kk], j, i, a_ji);
        assert(0);
      }
      a_ji = v[kk] = v[kk] / a_ii;

      for (kk = kk + 1; kk < off[jj + 1]; kk++) { /* Go over the columns */
        k = col[kk];
        assert(csr_mat_get(&a_ik, M, i, k) == 0);
        /* a_jk = a_jk - a_ji * a_ik */
        v[kk] = v[kk] - a_ji * a_ik;
      }
    }
  }
}

static int lu_solve(double *x, struct csr_mat_ *A, double *b, buffer *buf) {
  uint n = A->rn;
  buffer_reserve(buf, sizeof(double) * n);
  double *y = buf->ptr;

  sint i, j;
  /* Forward substitution with L - Lower matrix (diagonal = all 1s)
   * L(Ux) = L(y) = b */
  for (i = 0; i < n; i++) {
    y[i] = b[i];
    /* Remove the dot product between L[i, :] and currently determined solution
     */
    for (j = A->row_off[i]; j < A->row_off[i + 1] && A->col[j] < i + 1; j++)
      y[i] -= A->v[j] * y[A->col[j]];
  }

  /* Back substitution with U, Ux = y */
  for (i = n - 1; i >= 0; i--) {
    x[i] = y[i];
    /* Remove the dot product between L[i, :] and currently determined solution
     */
    for (j = A->row_off[i + 1] - 1; j >= A->row_off[i] && A->col[j] > i + 1;
         j--)
      x[i] -= A->v[j] * x[A->col[j]];
    x[i] /= A->v[j];
  }

  return 0;
}

int ilu_setup(genmap_handle h, struct comm *c, struct ilu_data *d) {
  d->M = tmalloc(struct csr_mat_, 1);

  csr_mat_copy(d->M, h->M);

  csr_mat_print(d->M, c);

  ilu0(d->M);

  csr_mat_print(d->M, c);

  return 0;
}

int ilu_apply(GenmapScalar *u, GenmapScalar *rhs, struct ilu_data *d,
              buffer *buf) {
  return lu_solve(u, d->M, rhs, buf);
}

int ilu_free(struct ilu_data *d) {
  csr_mat_free(d->M);
  free(d);

  return 0;
}
