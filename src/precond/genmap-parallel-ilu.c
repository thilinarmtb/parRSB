#include <genmap-ilu.h>
#include <genmap-impl.h>

int parallel_ilu_setup(genmap_handle h, struct comm *c,
                       struct parallel_ilu_data *d) {
  d->M = tmalloc(struct csr_mat_, 1);
  csr_mat_copy(d->M, h->M);
  d->h = h;

  parrsb_ilu0(d->h->nlevels, d->h->level_off, d->h->comms, d->M);

  return 0;
}

int parallel_ilu_apply(GenmapScalar *u, GenmapScalar *rhs,
                       struct parallel_ilu_data *d, buffer *buf) {
  parrsb_lu_solve(u, d->M, rhs, d->h->nlevels, d->h->level_off, d->h->comms);
  return 0;
}

int parallel_ilu_free(struct parallel_ilu_data *d) {
  csr_mat_free(d->M);
  free(d);
  return 0;
}
