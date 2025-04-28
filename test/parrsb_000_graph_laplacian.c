#include "parrsb_test.h"

static int test_brick_00(struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  unsigned n = 1;
  graph_init_brick(&g, n, c->c);

  element_info ei;
  graph_element_info_init(&ei);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c, bfr);
  scalar u = 1, v = 1;
  laplacian_op(&v, l, &u, bfr);
  laplacian_free(&l);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_element_info_free(&ei);
  graph_free(&g);

  parrsb_test_error(fabs(v) > PARRSB_TEST_EPS, c);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  int err = 0;
  err |= test_brick_00(&cr, &bfr);

  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
