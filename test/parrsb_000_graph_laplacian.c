#include "parrsb_test.h"

static int test_brick_00(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, 0, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  uint size = laplacian_size(l);
  laplacian_op(NULL, l, NULL);
  laplacian_free(&l);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(size == 0, c);
}

static int test_brick_05(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, 1, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  uint size = laplacian_size(l);
  scalar u = 1, v = 1;
  laplacian_op(&v, l, &u);
  laplacian_free(&l);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(size == 0 || fabs(v) < PARRSB_TEST_EPS, c);
}

static int test_brick_10(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  if (c->id >= c->np / 2)
    graph_init(&g, 1, c->c);
  else
    graph_init(&g, 0, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  uint size = laplacian_size(l);
  scalar u = 1, v = 1;
  laplacian_op(&v, l, &u);
  laplacian_free(&l);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(size == 0 || fabs(v) < PARRSB_TEST_EPS, c);
}

static int test_brick_15(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, c->id % 2, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  uint size = laplacian_size(l);
  scalar u = 1, v = 1;
  laplacian_op(&v, l, &u);
  laplacian_free(&l);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(size == 0 || fabs(v) < PARRSB_TEST_EPS, c);
}

static int test_brick_20(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

#define SIZE 2

  graph_t g;
  graph_init(&g, SIZE, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  scalar u[SIZE];
  scalar v[SIZE];
  for (uint i = 0; i < SIZE; i++) u[i] = 1;
  laplacian_op(v, l, u);
  laplacian_free(&l);

  scalar normv = norm2(v, SIZE);

#undef SIZE

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(normv < PARRSB_TEST_EPS, c);
}

static int test_brick_25(const element_info ei, struct crystal *cr,
                         buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, c->id, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);
  uint size = laplacian_size(l);
  scalar *u = tcalloc(scalar, size);
  scalar *v = tcalloc(scalar, size);
  for (uint i = 0; i < size; i++) u[i] = 1;
  laplacian_op(v, l, u);
  laplacian_free(&l);

  scalar normv = norm2(v, size);
  free(u), free(v);

  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_test(size == 0 || normv < PARRSB_TEST_EPS, c);
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
  err |= test_brick_00(ei, &cr, &bfr);
  err |= test_brick_05(ei, &cr, &bfr);
  err |= test_brick_10(ei, &cr, &bfr);
  err |= test_brick_15(ei, &cr, &bfr);
  err |= test_brick_20(ei, &cr, &bfr);
  err |= test_brick_25(ei, &cr, &bfr);

  graph_element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
