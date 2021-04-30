#include <genmap-impl.h>

int genmap_laplacian_init(genmap_handle h, struct comm *c) {
  int gs = h->options->rsb_laplacian_implementation & 1;
  int csr = h->options->rsb_laplacian_implementation & 2;

  if (gs)
    genmap_laplacian_gs_init(h, c);
  if (csr)
    genmap_laplacian_csr_init(h, c);
}

int genmap_laplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  int gs = h->options->rsb_laplacian_implementation & 1;
  int csr = h->options->rsb_laplacian_implementation & 2;

  if (gs)
    genmap_laplacian_gs(h, u, v);
  else if (csr)
    genmap_laplacian_csr(h, u, v);
}
