#include <genmap-multigrid-csr.h>
#include <genmap-multigrid-gs.h>
#include <genmap-multigrid.h>

int log2i(sint i) {
  sint k = 1, l = 0;
  while (k <= i)
    k *= 2, l++;
  return l - 1;
}

int mg_level_get_nsmooth(struct mg_level *lvl) { return lvl->get_nsmooth(lvl); }

GenmapScalar mg_level_get_sigma(struct mg_level *lvl) {
  return lvl->get_sigma(lvl);
}

void mg_setup(genmap_handle h, struct comm *c, struct mg_data *d) {
  // Currently only supporting csr
  d->get_nlevels = mg_get_nlevels_csr;
  d->get_level_off = mg_get_level_off_csr;
  d->setup = mg_setup_csr;
  d->free = mg_free_csr;

  d->setup(h, c, d);
}

int mg_get_nlevels(struct mg_data *d) { return d->get_nlevels(d); }

uint *mg_get_level_off(struct mg_data *d) { return d->get_level_off(d); }

void mg_free(struct mg_data *d) { d->free(d); }
