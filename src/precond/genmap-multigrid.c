#include <genmap-impl.h>
#include <genmap-multigrid.h>

int log2i(sint i) {
  sint k = 1, l = 0;
  while (k <= i)
    k *= 2, l++;
  return l - 1;
}

void mg_check(genmap_handle h, struct comm *c) {
  struct mg_data csr, gs;
  mg_setup_csr(h, c, &csr);
  mg_setup_gs(h, c, &gs);

  int nlevels = mg_get_nlevels(&csr);
  assert(nlevels == mg_get_nlevels(&gs));

  uint *csr_off = mg_get_level_off(&csr);
  uint *gs_off = mg_get_level_off(&gs);

  for (int i = 0; i < nlevels + 1; i++)
    assert(csr_off[i] == gs_off[i]);

  GenmapScalar *csr_r, *gs_r;
  GenmapMalloc(csr_off[2] - csr_off[0], &csr_r);
  GenmapMalloc(gs_off[2] - gs_off[0], &gs_r);

  GenmapFree(csr_r);
  GenmapFree(gs_r);
}

void mg_setup(genmap_handle h, struct comm *c, struct mg_data *d) {
  // FIXME: Select between gs and csr
  mg_setup_csr(h, c, d);
}

int mg_get_nlevels(struct mg_data *d) { return d->get_nlevels(d); }

uint *mg_get_level_off(struct mg_data *d) { return d->get_level_off(d); }

int mg_get_nsmooth(struct mg_data *d, int level) {
  return d->get_nsmooth(d, level);
}

GenmapScalar mg_get_sigma(struct mg_data *d, int level) {
  return d->get_sigma(d, level);
}

GenmapScalar *mg_get_diagonal(struct mg_data *d, int level) {
  return d->get_diagonal(d, level);
}

void mg_operator(GenmapScalar *v, GenmapScalar *u, int level, struct mg_data *d,
                 buffer *buf) {
  d->G(d, level, v, u, buf);
}

void mg_restrict(GenmapScalar *v, int level, struct mg_data *d, buffer *buf) {
  d->rstrct(d, level, v, buf);
}

void mg_interpolate(GenmapScalar *v, int level, struct mg_data *d,
                    buffer *buf) {
  d->intrp(d, level, v, buf);
}

void mg_free(struct mg_data *d) { d->free(d); }
