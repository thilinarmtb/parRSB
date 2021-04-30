// FIXME: genmap-impl.h is only for GenmapFree
#include <genmap-impl.h>
#include <genmap-multigrid-gs.h>
#include <genmap-partition.h>

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

static void mg_diagonal_scaling_gs(struct mg_data *d, int lvl, GenmapScalar *v,
                                   GenmapScalar *u, GenmapScalar sigma) {
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  genmap_handle h = d->h;

  uint nelt = dd->level_off[1] - dd->level_off[0];
  uint size = nelt + dd->level_off[lvl + 1] - dd->level_off[lvl];

  GenmapScalar *vv;
  GenmapCalloc(size, &vv);

  uint i;
  for (i = nelt; i < size; i++)
    vv[i] = u[i - nelt];

  // Interpolate to Level 0
  if (0 < lvl && lvl < dd->nlevels)
    gs(vv, gs_double, gs_add, 0, dd->J0[lvl], &h->buf);

  GenmapScalar *diag = dd->diagonal;

  int nv = genmap_get_nvertices(h);

  if (lvl > 0) {
    for (i = 0; i < nelt; i++)
      vv[i] = sigma * vv[i] / diag[i];
  } else {
    for (i = nelt; i < size; i++)
      vv[i] = sigma * vv[i] / diag[i - nelt];
  }

  // Restrict from Level 0
  if (0 < lvl && lvl < dd->nlevels) {
    for (i = nelt; i < size; i++)
      vv[i] = 0.0;
    gs(vv, gs_double, gs_add, 1, dd->J0[lvl], &h->buf);
  }

  for (i = nelt; i < size; i++)
    v[i - nelt] = vv[i];

  GenmapFree(vv);
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
  struct mg_data_gs *dd = d->data;
  assert(lvl < dd->nlevels);

  uint nelt = dd->level_off[1] - dd->level_off[0];
  uint size = nelt + dd->level_off[lvl + 1] - dd->level_off[lvl];

  GenmapScalar *vv;
  GenmapCalloc(size, &vv);

  uint i;
  for (i = nelt; i < size; i++)
    vv[i] = u[i - nelt];

  // Interpolate to Level 0
  if (0 < lvl && lvl < dd->nlevels)
    gs(vv, gs_double, gs_add, 0, dd->J0[lvl], buf);

  if (lvl > 0)
    genmap_laplacian_gs(d->h, vv, vv);
  else
    genmap_laplacian_gs(d->h, vv + nelt, vv + nelt);

  // Restrict from Level 0
  if (0 < lvl && lvl < dd->nlevels) {
    for (i = nelt; i < size; i++)
      vv[i] = 0.0;
    gs(vv, gs_double, gs_add, 1, dd->J0[lvl], buf);
  }

  for (i = nelt; i < size; i++)
    v[i - nelt] = vv[i];

  GenmapFree(vv);
}

static void mg_coarse_gs(struct mg_data *d, GenmapScalar *u, GenmapScalar *r) {
  struct mg_data_gs *dd = d->data;

  int nlevels = dd->nlevels;
  uint nelt = dd->level_off[1] - dd->level_off[0];
  uint off = dd->level_off[nlevels - 1];
  uint n = dd->level_off[nlevels] - off;
  uint size = nelt + n;

  GenmapScalar *vv;
  GenmapCalloc(size, &vv);

  GenmapScalar *diag = dd->diagonal;

  uint i;
  for (i = 0; i < nelt; i++)
    vv[i] = diag[i] - 2.0;

  buffer buf;
  buffer_init(&buf, 1024);

  // Interpolate to Level 0
  gs(vv, gs_double, gs_add, 1, dd->J0[nlevels - 1], &buf);

  buffer_free(&buf);

  if (n == 1) {
    if (fabs(vv[nelt]) > sqrt(GENMAP_TOL))
      u[off] = r[off] / vv[nelt];
    else
      u[off] = 0.0;
  }
}

static void mg_free_gs(struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;

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
  dd->diagonal_scaling = mg_diagonal_scaling_gs;
  dd->G = mg_operator_gs;
  dd->rstrct = mg_restrict_gs;
  dd->intrp = mg_interpolate_gs;
  dd->free = mg_free_gs;

  struct mg_data_gs *d = dd->data = calloc(1, sizeof(struct mg_data_gs));

  d->c = c;
  uint rn = genmap_get_nel(h);

  slong in = rn;
  slong out[2][1], bf[2][1];
  comm_scan(out, d->c, gs_long, gs_add, &in, 1, bf);
  slong rg = out[1][0];

  d->nlevels = log2i(rg) + 1;

  /* Allocate data structures */
  GenmapMalloc(d->nlevels + 1, &d->level_off);
  GenmapMalloc(d->nlevels, &d->nsmooth);
  GenmapMalloc(d->nlevels, &d->sigma);
  GenmapMalloc(d->nlevels, &d->J);
  GenmapMalloc(d->nlevels, &d->J0);

  /* Setup Level 0 */
  d->level_off[0] = 0;
  d->level_off[1] = rn;
  d->nsmooth[0] = 2;
  d->sigma[0] = 0.6;
  d->J[0] = NULL;
  d->diagonal = h->diagonal;
  d->nv = genmap_get_nvertices(h);

  /* Setup other levels */
  slong *wrk;
  size_t size = 2 * rn;
  GenmapMalloc(size, &wrk);
  mg_setup_aux_gs(d, wrk, &h->buf);
  GenmapFree(wrk);

  GenmapMalloc(rn, &d->buf);
}
