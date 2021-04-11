#include <genmap-impl.h>
#include <genmap-multigrid-csr.h>

int mg_level_get_nsmooth_csr(struct mg_level *lvl) {
  struct mg_level_csr *l = lvl->data;

  return l->nsmooth;
}

GenmapScalar mg_level_get_sigma_csr(struct mg_level *lvl) {
  struct mg_level_csr *l = lvl->data;

  return l->sigma;
}

void mg_setup_csr(genmap_handle h, struct comm *gsc, struct mg_data *dd) {
  struct mg_data_csr *d = dd->data = calloc(1, sizeof(struct mg_data_csr));

  comm_dup(&d->c, gsc);

  csr_mat M = h->M;
  uint rn = M->rn;

  slong out[2][1], bf[2][1], in = rn;
  comm_scan(out, &d->c, gs_long, gs_add, &in, 1, bf);
  slong rg = out[1][0];

  d->nlevels = log2i(rg) + 1;

  /* Setup Level 0 */
  GenmapMalloc(d->nlevels + 1, &d->level_off);
  d->level_off[0] = 0;
  d->level_off[1] = M->rn;

  GenmapMalloc(d->nlevels, &d->levels);
  GenmapMalloc(1, &d->levels[0]);
  d->levels[0]->M = M;
  d->levels[0]->nsmooth = 2;
  d->levels[0]->sigma = 0.6;
  d->levels[0]->J = NULL;

  /* Setup other levels */
  slong *wrk;
  GenmapMalloc(2 * rn, &wrk);
  mg_level_setup(d, wrk, &h->buf);
  GenmapFree(wrk);

  GenmapMalloc(d->level_off[d->nlevels], &d->x);
  GenmapMalloc(d->level_off[d->nlevels], &d->y);
  GenmapMalloc(d->level_off[d->nlevels], &d->b);
  GenmapMalloc(d->level_off[d->nlevels], &d->u);
  GenmapMalloc(d->level_off[d->nlevels], &d->rhs);

  uint nnz = M->row_off[M->rn];
  GenmapMalloc(nnz, &d->buf);
}

int mg_get_nlevels_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;
  return d->nlevels;
}

uint *mg_get_level_off_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;
  return d->level_off;
}

void mg_free_csr(struct mg_data *dd) {
  struct mg_data_csr *d = dd->data;

  mgLevel *l = d->levels;
  uint i;
  for (i = 0; i < d->nlevels; i++) {
    if (i > 0)
      csr_mat_free(l[i]->M);
    if (i < d->nlevels - 1)
      gs_free(l[i]->J);
    GenmapFree(l[i]);
  }

  comm_free(&d->c);

  GenmapFree(l);
  GenmapFree(d->level_off);
  GenmapFree(d->y);
  GenmapFree(d->x);
  GenmapFree(d->b);
  GenmapFree(d->buf);
  GenmapFree(d->rhs);
  GenmapFree(d->u);

  free(d);
}
