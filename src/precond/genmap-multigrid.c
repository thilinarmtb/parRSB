#include <genmap-multigrid.h>

int log2i(sint i) {
  sint k = 1, l = 0;
  while (k <= i)
    k *= 2, l++;
  return l - 1;
}

void mg_setup(genmap_handle h, struct comm *c, struct mg_data *d) {
  // Currently only supporting csr
  d->data = calloc(1, sizeof(struct mg_data_csr));
  mg_setup_csr(h, c, (struct mg_data_csr *)d->data);
}

void mg_free(struct mg_data *d) {
  mg_free_csr((struct mg_data_csr *)d->data);
  free(d->data);
}
