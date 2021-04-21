// FIXME: genmap-impl.h is only for GenmapFree
#include <genmap-impl.h>
#include <genmap-multigrid-csr.h>

static int mg_get_nlevels_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;
  return d->nlevels;
}

static uint *mg_get_level_off_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;
  return d->level_off;
}

static int mg_get_nsmooth_csr(struct mg_data *d, int lvl) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  return dd->nsmooth[lvl];
}

static GenmapScalar mg_get_sigma_csr(struct mg_data *d, int lvl) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  return dd->sigma[lvl];
}

static GenmapScalar *mg_get_diagonal_csr(struct mg_data *d, int lvl) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  return dd->M[lvl]->diag;
}

static void mg_restrict_csr(struct mg_data *d, int lvl, GenmapScalar *v,
                            buffer *buf) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  gs(v, gs_double, gs_add, 1, dd->J[lvl], buf);
}

static void mg_interpolate_csr(struct mg_data *d, int lvl, GenmapScalar *v,
                               buffer *buf) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  gs(v, gs_double, gs_add, 0, dd->J[lvl], buf);
}

static void mg_operator_csr(struct mg_data *d, int lvl, GenmapScalar *v,
                            GenmapScalar *u, buffer *buf) {
  struct mg_data_csr *dd = d->data;
  assert(lvl < dd->nlevels);

  csr_mat M = dd->M[lvl];
  csr_mat_gather(M, M->gsh, u, dd->buf, buf);
  csr_mat_apply(v, M, dd->buf);
}

void mg_free_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;

  comm_free(&d->c);

  sint i;
  for (i = 0; i < d->nlevels; i++) {
    if (d->J[i] != NULL)
      gs_free(d->J[i]);
    if (i > 0 && d->M[i] != NULL)
      csr_mat_free(d->M[i]);
  }

  GenmapFree(d->nsmooth);
  GenmapFree(d->sigma);
  GenmapFree(d->J);
  GenmapFree(d->M);

  GenmapFree(d->level_off);
  GenmapFree(d->buf);

  free(d);
}

void mg_setup_csr(genmap_handle h, struct comm *gsc, struct mg_data *dd) {
  dd->get_nlevels = mg_get_nlevels_csr;
  dd->get_level_off = mg_get_level_off_csr;
  dd->get_nsmooth = mg_get_nsmooth_csr;
  dd->get_sigma = mg_get_sigma_csr;
  dd->get_diagonal = mg_get_diagonal_csr;
  dd->G = mg_operator_csr;
  dd->rstrct = mg_restrict_csr;
  dd->intrp = mg_interpolate_csr;
  dd->free = mg_free_csr;

  struct mg_data_csr *d = dd->data = calloc(1, sizeof(struct mg_data_csr));

  // FIXME: May be this is not necessary
  comm_dup(&d->c, gsc);

  uint rn = h->M->rn;

  slong in = rn;
  slong out[2][1], bf[2][1];
  comm_scan(out, &d->c, gs_long, gs_add, &in, 1, bf);
  slong rg = out[1][0];

  d->nlevels = log2i(rg) + 1;

  /* Allocate data structures */
  GenmapMalloc(d->nlevels + 1, &d->level_off);
  GenmapMalloc(d->nlevels, &d->nsmooth);
  GenmapMalloc(d->nlevels, &d->sigma);
  GenmapMalloc(d->nlevels, &d->J);
  GenmapMalloc(d->nlevels, &d->M);

  /* Setup Level 0 */
  d->level_off[0] = 0;
  d->level_off[1] = rn;
  d->nsmooth[0] = 2;
  d->sigma[0] = 0.6;
  d->J[0] = NULL;
  d->M[0] = h->M;

  /* Setup other levels */
  slong *wrk;
  GenmapMalloc(2 * rn, &wrk);
  mg_setup_aux_csr(d, wrk, &h->buf);
  GenmapFree(wrk);

  uint nnz = h->M->row_off[rn];
  GenmapMalloc(nnz, &d->buf);
}
