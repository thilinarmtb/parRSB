#include "parrsb_test.h"

static int test_laplacian(const element_info ei, unsigned n, graph_type_t type,
                          struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, n, type, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);

  uint size = laplacian_size(l);
  scalar *u = tcalloc(scalar, size);
  scalar *v = tcalloc(scalar, size);
  for (uint i = 0; i < size; i++) v[i] = u[i] = 1;

  if (size == 0)
    laplacian_op(NULL, l, NULL);
  else
    laplacian_op(v, l, u);

  scalar normv = norm2(v, size, c);

  free(u), free(v);
  laplacian_free(&l);
  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  parrsb_assert(size == 0 || normv < PARRSB_TEST_TOL, c);
}

static int test_path_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 0, PATH, cr, bfr);
}

static int test_ring_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 0, RING, cr, bfr);
}

static int test_complete_00(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_laplacian(ei, 0, COMPLETE, cr, bfr);
}

static int test_path_05(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 1, PATH, cr, bfr);
}

static int test_ring_05(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 1, RING, cr, bfr);
}

static int test_complete_05(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_laplacian(ei, 1, COMPLETE, cr, bfr);
}

static int test_path_10(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 2, PATH, cr, bfr);
}

static int test_ring_10(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_laplacian(ei, 2, RING, cr, bfr);
}

static int test_complete_10(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_laplacian(ei, 2, COMPLETE, cr, bfr);
}

static int test_path_15(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id % 2, PATH, cr, bfr);
}

static int test_ring_15(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id % 2, RING, cr, bfr);
}

static int test_complete_15(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id % 2, COMPLETE, cr, bfr);
}

static int test_path_20(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id >= c->np / 2, PATH, cr, bfr);
}

static int test_ring_20(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id >= c->np / 2, RING, cr, bfr);
}

static int test_complete_20(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id >= c->np / 2, COMPLETE, cr, bfr);
}

static int test_path_25(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id, PATH, cr, bfr);
}

static int test_ring_25(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id, RING, cr, bfr);
}

static int test_complete_25(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_laplacian(ei, c->id, COMPLETE, cr, bfr);
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
  element_info_init(&ei, GRAPH);

  int err = 0;
  err |= test_path_00(ei, &cr, &bfr);
  err |= test_path_05(ei, &cr, &bfr);
  err |= test_path_10(ei, &cr, &bfr);
  err |= test_path_15(ei, &cr, &bfr);
  err |= test_path_20(ei, &cr, &bfr);
  err |= test_path_25(ei, &cr, &bfr);
  err |= test_ring_00(ei, &cr, &bfr);
  err |= test_ring_05(ei, &cr, &bfr);
  err |= test_ring_10(ei, &cr, &bfr);
  err |= test_ring_15(ei, &cr, &bfr);
  err |= test_ring_20(ei, &cr, &bfr);
  err |= test_ring_25(ei, &cr, &bfr);
  err |= test_complete_00(ei, &cr, &bfr);
  err |= test_complete_05(ei, &cr, &bfr);
  err |= test_complete_10(ei, &cr, &bfr);
  err |= test_complete_15(ei, &cr, &bfr);
  err |= test_complete_20(ei, &cr, &bfr);
  err |= test_complete_25(ei, &cr, &bfr);

  element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
