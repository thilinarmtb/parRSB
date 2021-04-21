// FIXME: genmap-impl.h is only for GenmapFree
#include <genmap-impl.h>
#include <genmap-multigrid-gs.h>

static int mg_get_nlevels_gs(struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;
  return d->nlevels;
}

static uint *mg_get_level_off_gs(struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;
  return d->level_off;
}

static int mg_get_nsmooth_gs(struct mg_data *d, int lvl) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  return dd->nsmooth[lvl];
}

static GenmapScalar mg_get_sigma_gs(struct mg_data *d, int lvl) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  return dd->sigma[lvl];
}

static GenmapScalar *mg_get_diagonal_gs(struct mg_data *d, int lvl) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  return NULL;
}

static void mg_restrict_gs(struct mg_data *d, int lvl, GenmapScalar *v,
                           buffer *buf) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  gs(v, gs_double, gs_add, 1, dd->J[lvl], buf);
}

static void mg_interpolate_gs(struct mg_data *d, int lvl, GenmapScalar *v,
                              buffer *buf) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  gs(v, gs_double, gs_add, 0, dd->J[lvl], buf);
}

static void mg_operator_gs(struct mg_data *d, int lvl, GenmapScalar *v,
                           GenmapScalar *u, buffer *buf) {
  if (lvl > 0) {
    // Interpolate
  }

  //gs(....);

  if (lvl > 0) {
    // Restrict
  }
}

static void mg_free_gs(struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;

  comm_free(&d->c);

  GenmapFree(d->nsmooth);
  GenmapFree(d->sigma);
  GenmapFree(d->J);

  GenmapFree(d->level_off);
  GenmapFree(d->buf);

  free(d);
}

void mg_setup_gs(genmap_handle h, struct comm *c, struct mg_data *dd) {
  dd->get_nlevels = mg_get_nlevels_gs;
  dd->get_level_off = mg_get_level_off_gs;
  dd->get_nsmooth = mg_get_nsmooth_gs;
  dd->get_sigma = mg_get_sigma_gs;
  dd->get_diagonal = mg_get_diagonal_gs;
  dd->G = mg_operator_gs;
  dd->rstrct = mg_restrict_gs;
  dd->intrp = mg_interpolate_gs;
  dd->free = mg_free_gs;

  struct mg_data_gs *d = dd->data = calloc(1, sizeof(struct mg_data_gs));

  // FIXME: May be this is not necessary
  comm_dup(&d->c, c);

  uint rn = genmap_get_nel(h);

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
  GenmapMalloc(d->nlevels, &d->R0);

  /* Setup Level 0 */
  d->level_off[0] = 0;
  d->level_off[1] = rn;
  d->nsmooth[0] = 2;
  d->sigma[0] = 0.6;
  d->J[0] = NULL;
  d->G = h->gs;
  d->diagonal = h->diagonal;

  /* Setup other levels */
  slong *wrk;
  GenmapMalloc(2 * rn, &wrk);
  mg_setup_aux_gs(d, wrk, &h->buf);
  GenmapFree(wrk);

  GenmapMalloc(rn, &d->buf);
}
