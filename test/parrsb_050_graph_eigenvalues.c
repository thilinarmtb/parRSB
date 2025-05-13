#include "parrsb_test.h"

static scalar test_base(const element_info ei, unsigned n, graph_type_t type,
                        struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  graph_t g;
  graph_init(&g, n, type, c->c);

  struct array nlist;
  graph_load_balance(&nlist, g->n, g->nodes, g->offsets, g->neighbors, cr, bfr);

  laplacian l;
  laplacian_init(&l, &nlist, ei, c);

  uint size = laplacian_size(l);

  slong out[2][1], wrk[2][1], in = size;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1;

  scalar *f = tcalloc(scalar, 2 * size), *v = f + size;
  for (uint i = 0; i < size; i++) v[i] = start + i;

  parrsb_options opts;
  parrsb_options_get_default(&opts);
  lanczos(f, l, v, c, opts);
  parrsb_options_free(&opts);

  laplacian_op(v, l, f);
  scalar normf = norm2(f, size, c);
  scalar normv = norm2(v, size, c);

  free(f);
  laplacian_free(&l);
  graph_restore(NULL, &nlist, cr, bfr);
  graph_free(&g);

  if (fabs(normf) < PARRSB_TEST_EPS)
    return 0;
  else
    return normv / normf;
}

static int test_path_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  sint n = 3;
  scalar evc = test_base(ei, n, PATH, cr, bfr);

  const struct comm *c = &(cr->comm);
  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);

  scalar eva = 2 * (1 - cos((M_PI * 1) / ng));

  parrsb_test(((fabs(evc - eva) / eva) < PARRSB_TEST_EPS), c);
}

static int test_ring_00(const element_info ei, struct crystal *cr,
                        buffer *bfr) {
  sint n = 3;
  scalar evc = test_base(ei, n, RING, cr, bfr);

  const struct comm *c = &(cr->comm);
  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);

  scalar eva = 2 * (1 - cos((2 * M_PI * 1) / ng));

  parrsb_test(((fabs(evc - eva) / eva) < PARRSB_TEST_EPS), c);
}

static int test_complete_00(const element_info ei, struct crystal *cr,
                            buffer *bfr) {
  sint n = 3;
  scalar evc = test_base(ei, n, COMPLETE, cr, bfr);

  const struct comm *c = &(cr->comm);
  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);
  scalar eva = ng;

  parrsb_test(((fabs(evc - eva) / eva) < PARRSB_TEST_EPS), c);
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
  err |= test_ring_00(ei, &cr, &bfr);
  err |= test_complete_00(ei, &cr, &bfr);

  graph_element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
