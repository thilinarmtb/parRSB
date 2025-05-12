#include "parrsb_test.h"

static int test_base(const element_info ei, unsigned n, graph_type_t type,
                     struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, n, type, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);

  uint size = laplacian_size(l);
  scalar *u = 0, *v = 0;

  if (size == 0) {
    laplacian_op(NULL, l, NULL);
    goto check_and_finalize;
  } else {
    u = tcalloc(scalar, 2 * size), v = u + size;
    for (unsigned i = 0; i < 2 * size; i++) u[i] = 1;

    laplacian_op(v, l, u);
  }

check_and_finalize:
  scalar normv = norm2(v, size);
  free(u);
  laplacian_free(&l);
  parrsb_test(size == 0 || normv < PARRSB_TEST_EPS, c);
}

static int test_path_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 0, PATH, cr, bfr);
}

static int test_ring_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 0, RING, cr, bfr);
}

static int test_complete_00(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_base(ei, 0, COMPLETE, cr, bfr);
}

static int test_path_05(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 1, PATH, cr, bfr);
}

static int test_ring_05(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 1, RING, cr, bfr);
}

static int test_complete_05(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_base(ei, 1, COMPLETE, cr, bfr);
}

static int test_path_10(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 2, PATH, cr, bfr);
}

static int test_ring_10(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  return test_base(ei, 2, RING, cr, bfr);
}

static int test_complete_10(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  return test_base(ei, 2, COMPLETE, cr, bfr);
}

static int test_path_15(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id % 2, PATH, cr, bfr);
}

static int test_ring_15(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id % 2, RING, cr, bfr);
}

static int test_complete_15(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id % 2, COMPLETE, cr, bfr);
}

static int test_path_20(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id >= c->np / 2, PATH, cr, bfr);
}

static int test_ring_20(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id >= c->np / 2, RING, cr, bfr);
}

static int test_complete_20(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id >= c->np / 2, COMPLETE, cr, bfr);
}

static int test_path_25(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id, PATH, cr, bfr);
}

static int test_ring_25(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id, RING, cr, bfr);
}

static int test_complete_25(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  const struct comm *c = &(cr->comm);
  return test_base(ei, c->id, COMPLETE, cr, bfr);
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

  graph_element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
