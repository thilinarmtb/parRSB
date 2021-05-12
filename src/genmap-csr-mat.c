#include <genmap-impl.h>
#include <genmap-multigrid.h>

static struct gs_data *csr_get_top(csr_mat M, struct comm *c, buffer *buf) {
  const uint rn = M->rn;
  const uint n = M->row_off[rn];

  buffer_reserve(buf, sizeof(slong) * n);
  slong *ids = buf->ptr;

  uint i, j;
  for (i = 0; i < rn; i++)
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++)
      if (M->row_start + i == M->col[j])
        ids[j] = M->col[j];
      else
        ids[j] = -M->col[j];

  struct gs_data *gsh = gs_setup(ids, n, c, 0, gs_pairwise, 0);

  return gsh;
}

void csr_mat_setup(csr_mat M, struct array *entries, struct comm *c,
                   buffer *buf) {
  csr_entry *ptr = entries->ptr;
  sarray_sort_2(csr_entry, ptr, entries->n, r, 1, c, 1, buf);

  uint i = 0, j, n = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = j, n++;
  }

  M->rn = n;

  slong out[2][1], bf[2][1];
  slong in = M->rn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  M->row_start = out[0][0] + 1;

  GenmapMalloc(M->rn + 1, &M->row_off);

  if (n == 0) {
    M->col = NULL;
    M->v = NULL;
    M->diag = NULL;
  } else {
    GenmapMalloc(entries->n, &M->col);
    GenmapMalloc(entries->n, &M->v);
    GenmapMalloc(M->rn, &M->diag);
  }

  ptr = entries->ptr;
  uint rn = 0;
  for (i = 0; i < entries->n; i++) {
    M->col[i] = ptr[i].c, M->v[i] = ptr[i].v;
    if (ptr[i].r == ptr[i].c)
      M->diag[rn++] = ptr[i].v;
  }
  if (rn != M->rn) {
    printf("rn = %u, M->rn = %u\n", rn, M->rn);
    assert(0);
  }

  M->row_off[0] = 0, i = 0;
  uint nn = 0;
  while (i < entries->n) {
    j = i + 1;
    while (j < entries->n && ptr[i].r == ptr[j].r)
      j++;
    i = M->row_off[++nn] = j;
  }
  assert(n == nn);
  assert(M->row_off[n] == entries->n);

  M->gsh = csr_get_top(M, c, buf);
}

static void csr_mat_gather(csr_mat M, struct gs_data *gsh, GenmapScalar *x,
                           GenmapScalar *buf, buffer *bfr) {
  ulong s = M->row_start;
  sint i, j;
  for (i = 0; i < M->rn; i++)
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++)
      if (M->col[j] == s + i)
        buf[j] = x[i];
      else
        buf[j] = 0.0;

  gs(buf, gs_scalar, gs_add, 0, gsh, bfr);
}

void csr_mat_apply(GenmapScalar *v, csr_mat M, GenmapScalar *u, buffer *buf) {
  const uint rn = M->rn;
  const uint nnz = M->row_off[rn];

  GenmapScalar *x;
  GenmapCalloc(nnz, &x);

  csr_mat_gather(M, M->gsh, u, x, buf);

  if (rn == 0)
    return;

  const uint *offsets = M->row_off;
  const GenmapScalar *val = M->v;

  uint i, j;
  for (i = 0; i < rn; i++) {
    v[i] = 0.0;
    for (j = offsets[i]; j < offsets[i + 1]; j++)
      v[i] += (*val++) * (*x++);
  }

  GenmapFree(x - nnz);
}

void csr_mat_print(csr_mat M, struct comm *c) {
  const sint rn = M->rn;
  const uint *offsets = M->row_off;
  const GenmapScalar *v = M->v;
  const ulong *col = M->col;

  uint i, j, k;

  for (k = 0; k < c->np; k++) {
    comm_barrier(c);
    if (c->id == k) {
      for (i = 0; i < rn; i++) {
        for (j = offsets[i]; j < offsets[i + 1]; j++)
          fprintf(stderr, "(%lu,%lu) -> %.10lf\n", M->row_start + i, col[j],
                  v[j]);
      }
    }
    fflush(stderr);
  }
}

int csr_mat_get(double *val, csr_mat M, uint i, uint j) {
  if (i > M->rn || j > M->rn || i == 0 || j == 0) {
    printf("%s:%d: %u %u %u\n", __FILE__, __LINE__, i, j, M->rn);
    return 1;
  }

  /* TODO: Use binary search */
  uint s;
  for (s = M->row_off[i - 1]; s < M->row_off[i]; s++)
    if (M->col[s] == j) {
      *val = M->v[s];
      break;
    }

  if (s == M->row_off[i])
    *val = 0.0;

  return 0;
}

int csr_mat_copy(csr_mat D, csr_mat S) {
  /* TODO: Check for error */
  int err;

  D->row_start = S->row_start;
  uint rn = D->rn = S->rn;

  D->row_off = calloc(rn + 1, sizeof(uint));
  memcpy(D->row_off, S->row_off, sizeof(uint) * (rn + 1));

  D->diag = calloc(rn, sizeof(double));
  memcpy(D->diag, S->diag, sizeof(double) * rn);

  uint nnz = S->row_off[rn];

  D->col = calloc(nnz, sizeof(ulong));
  memcpy(D->col, S->col, sizeof(ulong) * nnz);

  D->v = calloc(nnz, sizeof(double));
  memcpy(D->v, S->v, sizeof(double) * nnz);

  D->gsh = NULL;

  return 0;
}

int csr_mat_free(csr_mat M) {
  if (M->col)
    GenmapFree(M->col);
  if (M->v)
    GenmapFree(M->v);
  if (M->diag)
    GenmapFree(M->diag);
  if (M->row_off)
    GenmapFree(M->row_off);
  if (M->gsh)
    gs_free(M->gsh);
  GenmapFree(M);

  return 0;
}
