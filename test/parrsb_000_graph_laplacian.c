#include "parrsb_test.h"

static void test_base(laplacian *l, const element_info ei, unsigned n,
                      graph_type_t type, struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, n, type, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian_init(l, &nlist, ei, c);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);
}

static int test_empty_path(const element_info ei, struct crystal *cr,
                           buffer *bfr) {
  laplacian l;
  test_base(&l, ei, 0, PATH, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_test(n == 0, c);
}

static int test_empty_ring(const element_info ei, struct crystal *cr,
                           buffer *bfr) {
  laplacian l;
  test_base(&l, ei, 0, RING, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_test(n == 0, c);
}

static int test_empty_complete(const element_info ei, struct crystal *cr,
                               buffer *bfr) {
  laplacian l;
  test_base(&l, ei, 0, COMPLETE, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_test(n == 0, c);
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

  element_info ei;
  graph_element_info_init(&ei);

  int err = 0;
  err |= test_empty_path(ei, &cr, &bfr);
  err |= test_empty_ring(ei, &cr, &bfr);
  err |= test_empty_complete(ei, &cr, &bfr);

  graph_element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
