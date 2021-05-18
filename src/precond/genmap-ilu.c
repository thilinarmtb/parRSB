#include <genmap-ilu.h>
#include <genmap-impl.h>

static int ilu0(csr_mat M, csr_mat N, struct comm *c) {
  csr_mat_copy(M, N);

  const uint rn = M->rn;
  assert(rn > 1);

  const uint *off = M->row_off;
  double *v = M->v;
  const ulong *col = M->col;

  double a_kk, a_kj, a_ik;
  uint i, j, k, kk;

  for (i = 2; i <= rn; i++) { /* Go over number of rows */
    for (k = 1; k < i; k++) {
      a_kk = M->diag[k - 1];

      assert(csr_mat_get_local(&a_ik, &kk, M, i, k) == 0);
      if (fabs(a_ik) < 1e-12)
        continue;
      a_ik = v[kk] = v[kk] / a_kk;

      /* a_ij = a_ij - a_ik * a_kj */
      for (kk++; kk < off[i]; kk++) { /* Go over the columns */
        j = col[kk];
        assert(csr_mat_get_local(&a_kj, NULL, M, k, j) == 0);
        v[kk] = v[kk] - a_ik * a_kj;
      }
    }
  }
}

static int ilut(csr_mat M, csr_mat N, struct comm *c, GenmapScalar threshold,
                buffer *buf) {
  assert(fabs(threshold) < 1.0);

  const uint rn = N->rn;
  assert(rn > 1);

  const uint *off = N->row_off;
  double *v = N->v;
  const ulong *col = N->col;

  /* M will have same dimensions as N, we will realloc v and col as we go */
  M->rn = rn;
  GenmapCalloc(rn + 1, &M->row_off);
  GenmapCalloc(rn, &M->diag);
  GenmapCalloc(1, &M->v);
  GenmapCalloc(1, &M->col);
  M->gsh = NULL;

  buffer_reserve(buf, sizeof(GenmapScalar) * rn);
  GenmapScalar *w = buf->ptr;

  double a_kk, t_i;
  uint i, j, jj, k, kk, km1, nnz = 0;

  for (i = 0; i < rn; i++) {
    /* Set w = a_{i, *}  and t_i */
    for (j = 0; j < rn; j++)
      w[j] = 0.0;
    t_i = 0.0;
    for (j = off[i]; j < off[i + 1]; j++) {
      w[col[j] - 1] = v[j];
      t_i += v[j] * v[j];
    }
    t_i = sqrt(t_i) * threshold;

    for (kk = off[i]; kk < off[i + 1] && col[kk] - 1 < i; kk++) {
      k = col[kk];
      km1 = k - 1;
      if (fabs(w[km1]) < 1e-12)
        continue;

      /* To-do: store diagonal separately */
      assert(csr_mat_get_local(&a_kk, NULL, N, k, k) == 0);
      w[km1] = w[km1] / a_kk;

      /* Apply first dropping rule */
      if (fabs(w[km1]) > t_i) {
        jj = M->row_off[km1];
        while (jj < M->row_off[k] && M->col[jj] < k)
          jj++;

        for (; jj < M->row_off[k]; jj++) {
          j = M->col[jj] - 1;
          w[j] = w[j] - w[km1] * M->v[jj];
        }
      } else
        w[km1] = 0.0;
    }

    /* Apply second dropping rule */
    for (k = 0; k < rn; k++)
      if (fabs(w[k]) > t_i)
        nnz++;

    j = M->row_off[i];
    GenmapRealloc(nnz, &M->v);
    GenmapRealloc(nnz, &M->col);
    for (k = 0; k < rn; k++) {
      if (fabs(w[k]) > t_i) {
        M->col[j] = k + 1;
        M->v[j] = w[k];
        j++;
      }
    }
    assert(j == nnz);
    M->row_off[i + 1] = nnz;
  }

  return 0;
}

/* FIXME: Store L and U seprately */
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
      y[i] -= A->v[j] * y[A->col[j] - 1];
  }

  /* Back substitution with U, Ux = y */
  for (i = n - 1; i >= 0; i--) {
    x[i] = y[i];
    /* Remove the dot product between L[i, :] and currently determined solution
     */
    for (j = A->row_off[i + 1] - 1; j >= A->row_off[i] && A->col[j] > i + 1;
         j--)
      x[i] -= A->v[j] * x[A->col[j] - 1];
    x[i] /= A->v[j];
  }

  return 0;
}

int ilu_setup(genmap_handle h, struct comm *c, struct ilu_data *d) {
  d->M = tmalloc(struct csr_mat_, 1);

#if 1
  ilu0(d->M, h->M, c);
#else
  ilut(d->M, h->M, c, 0.1, &h->buf);
#endif

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
