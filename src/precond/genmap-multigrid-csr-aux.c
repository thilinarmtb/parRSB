#include <genmap-impl.h>
#include <genmap-multigrid-csr.h>

#define GETPTR(p, i, off) ((char *)(p) + (off) + (i) * sizeof(entry))

static void csr_to_array(struct array *entries, csr_mat M) {
  uint rn = M->rn;

  uint nnz = M->row_off[rn];
  array_init(entry, entries, nnz);
  entries->n = nnz;

  uint i, j, nn = 0;
  entry *ptr = entries->ptr;
  for (i = 0; i < rn; i++)
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++) {
      ptr[nn].r = M->row_start + i;
      ptr[nn].c = M->col[j];
      ptr[nn].v = M->v[j];
      nn++;
    }

  assert(nn == nnz);
}

static void set_owner(struct array *coarse, size_t in_off, size_t out_off,
                      slong ng, sint np) {

  slong lelg = ng / 2;
  np = (lelg < np) ? lelg : np;

  uint lelt = lelg / np;
  uint nrem = lelg - lelt * np;

  ulong *in_ptr;
  sint *out_ptr;
  sint i;
  slong row;
  for (i = 0; i < coarse->n; i++) {
    in_ptr = (ulong *)GETPTR(coarse->ptr, i, in_off);
    out_ptr = (sint *)GETPTR(coarse->ptr, i, out_off);
    row = *in_ptr - 1;
    // FIXME: Assumes the 'reverse-Nek' element distribution
    // if (row<lelt*(np-nrem)) *out_ptr=(sint) row/lelt;
    // else *out_ptr=np-nrem+(sint) (row-lelt*(np-nrem))/(lelt+1);
    if (nrem == 0)
      *out_ptr = (sint)row / lelt;
    else if (row < (lelt + 1) * nrem)
      *out_ptr = (sint)row / (lelt + 1);
    else
      *out_ptr = nrem + (sint)(row - (lelt + 1) * nrem) / lelt;
  }
}

static void calculate_coarse_indices(slong ng, struct array *entries) {
  uint i;
  entry *ptr = entries->ptr;
  for (i = 0; i < entries->n; i++) {
    ptr[i].rn = (ptr[i].r + 1) / 2;
    ptr[i].cn = (ptr[i].c + 1) / 2;
  }

  if (ng % 2 == 1)
    for (i = 0; i < entries->n; i++) {
      if (ptr[i].c == ng)
        ptr[i].cn -= 1;
      if (ptr[i].r == ng)
        ptr[i].rn -= 1;
    }
}

static void coarsen_col(struct array *coarse, struct array *entries) {
  GenmapScalar v;
  sint i, j;

  i = 0;
  entry *ptr = entries->ptr;
  while (i < entries->n) {
    v = ptr[i].v, j = i + 1;
    while (j < entries->n && ptr[j].r == ptr[i].r && ptr[j].cn == ptr[i].cn)
      v += ptr[j].v, ptr[j].v = 0.0, j++;
    ptr[i].v = v;
    array_cat(entry, coarse, &ptr[i], 1);
    i = j;
  }
}

static csr_mat create_matrix_from_array(struct array *coarse, struct comm *c,
                                          buffer *buf) {
  uint size = coarse->n;
  sarray_sort_2(entry, coarse->ptr, size, rn, 1, cn, 1, buf);

  uint i, j, nn;
  i = j = nn = 0;
  entry *ptr = coarse->ptr;
  while (i < size) {
    while (j < size && ptr[j].rn == ptr[i].rn)
      j++;
    i = j, nn++;
  }

  uint nnz1;
  i = j = nnz1 = 0;
  while (i < size) {
    while (j < size && ptr[j].rn == ptr[i].rn && ptr[j].cn == ptr[i].cn)
      j++;
    i = j, nnz1++;
  }

  /* Create the matrix */
  csr_mat M1;
  GenmapMalloc(1, &M1);
  M1->rn = nn;
  GenmapMalloc(M1->rn + 1, &M1->row_off);

  slong out[2][1], bf[2][1];
  slong in = nn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
  M1->row_start = out[0][0] + 1;

  if (nnz1 == 0) {
    M1->col = NULL;
    M1->v = NULL;
    M1->diag = NULL;
  } else {
    GenmapMalloc(nnz1, &M1->col);
    GenmapMalloc(nnz1, &M1->v);
    GenmapMalloc(M1->rn, &M1->diag);
  }

  uint rn1;
  GenmapScalar v;

  M1->row_off[0] = nn = rn1 = 0;
  i = j = 0;
  while (i < size) {
    v = 0.0;
    while (j < size && ptr[j].rn == ptr[i].rn && ptr[j].cn == ptr[i].cn) {
      v += ptr[j].v;
      j++;
    }
    M1->col[nn] = ptr[i].cn, M1->v[nn] = v, nn++;

    if ((j < size && ptr[j].rn != ptr[i].rn) || j >= size)
      M1->row_off[++rn1] = nn;
    i = j;
  }
  assert(nn == nnz1);    // Sanity check
  assert(rn1 == M1->rn); // Sanity check

  /* Setup gs handle for the mat-vec */
  slong *ids;
  GenmapCalloc(nnz1, &ids);
  for (i = 0; i < M1->rn; i++)
    for (j = M1->row_off[i]; j < M1->row_off[i + 1]; j++)
      if (M1->row_start + i == M1->col[j])
        ids[j] = M1->col[j];
      else
        ids[j] = -M1->col[j];

  M1->gsh = gs_setup(ids, nnz1, c, 0, gs_pairwise, 0);

  GenmapFree(ids);

  return M1;
}

void mg_level_setup(struct mg_data_csr *d, slong *wrk, buffer *buf) {
  struct comm *c = &d->c;

  struct array coarse, entries;
  slong out[2][1], bf[2][1];
  uint lvl = 1;
  while (lvl < d->nlevels) {
    csr_mat M0 = d->levels[lvl - 1]->M;
    uint rn0 = M0->rn;

    slong in = rn0;
    comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
    slong ng = out[1][0];
    if (ng == 1)
      return;

    /* Coarsen the columns - easy since it is local */
    csr_to_array(&entries, M0);
    calculate_coarse_indices(ng, &entries);
    sarray_sort_2(entry, entries.ptr, entries.n, r, 1, cn, 1, buf);
    array_init(entry, &coarse, rn0);
    coarsen_col(&coarse, &entries);

    /* Setup gs ids for fine level (rhs interpolation) */
    uint nnz0 = entries.n;
    entry *ptr = entries.ptr;
    uint i, j;
    for (i = j = 0; i < nnz0; i++)
      if (ptr[i].r == ptr[i].c)
        wrk[j++] = -ptr[i].cn;
    assert(j == rn0);

    /* Send the rows to the correct processor */
    set_owner(&coarse, offsetof(entry, rn), offsetof(entry, p), ng, c->np);

    struct crystal cr;
    crystal_init(&cr, c);
    sarray_transfer(entry, &coarse, p, 1, &cr);
    crystal_free(&cr);

    /* Setup the CSR matrix by coarsening the rows */
    csr_mat M1 = create_matrix_from_array(&coarse, c, buf);

    /* Setup gs ids for coarse level (rhs interpolation ) */
    uint nn;
    for (i = nn = 0; i < M1->rn; i++)
      for (j = M1->row_off[i]; j < M1->row_off[i + 1]; j++)
        if (M1->row_start + i == M1->col[j]) {
          wrk[rn0 + nn] = M1->col[j];
          M1->diag[i] = M1->v[j];
          nn++;
        }
    assert(nn == M1->rn);

    struct gs_data *J = gs_setup(wrk, rn0 + M1->rn, c, 0, gs_crystal_router, 0);
    d->levels[lvl - 1]->J = J;

    GenmapMalloc(1, &d->levels[lvl]);
    mgLevel l = d->levels[lvl];
    l->M = M1;
    l->nsmooth = 2;
    l->sigma = 0.6;
    l->J = NULL;

    d->level_off[lvl + 1] = d->level_off[lvl] + M1->rn;

    array_free(&entries);
    array_free(&coarse);

    lvl++;
  }
}

#undef GETPTR
