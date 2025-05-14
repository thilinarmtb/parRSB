#include "parrsb_test.h"

static void test_laplacian_init(laplacian *l, const element_info ei, unsigned n,
                                graph_type_t type, struct crystal *cr,
                                buffer *bfr) {
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
  test_laplacian_init(&l, ei, 0, PATH, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_assert(n == 0, c);
}

static int test_empty_ring(const element_info ei, struct crystal *cr,
                           buffer *bfr) {
  laplacian l;
  test_laplacian_init(&l, ei, 0, RING, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_assert(n == 0, c);
}

static int test_empty_complete(const element_info ei, struct crystal *cr,
                               buffer *bfr) {
  laplacian l;
  test_laplacian_init(&l, ei, 0, COMPLETE, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  const struct comm *c = &(cr->comm);
  parrsb_assert(n == 0, c);
}

static int test_single_node_path(const element_info ei, struct crystal *cr,
                                 buffer *bfr) {
  const struct comm *c = &(cr->comm);

  laplacian l;
  test_laplacian_init(&l, ei, c->id == 0, PATH, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);
  parrsb_assert(ng == 1, c);
}

static int test_single_node_ring(const element_info ei, struct crystal *cr,
                                 buffer *bfr) {
  const struct comm *c = &(cr->comm);

  laplacian l;
  test_laplacian_init(&l, ei, c->id == 0, RING, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);
  parrsb_assert(ng == 1, c);
}

static int test_single_node_complete(const element_info ei, struct crystal *cr,
                                     buffer *bfr) {
  const struct comm *c = &(cr->comm);

  laplacian l;
  test_laplacian_init(&l, ei, c->id == 0, COMPLETE, cr, bfr);
  uint n = laplacian_size(l);
  laplacian_free(&l);

  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);
  parrsb_assert(ng == 1, c);
}

static int test_multiple_node_path(const element_info ei, struct crystal *cr,
                                   buffer *bfr) {
  const struct comm *c = &(cr->comm);

  sint n = 2;

  laplacian l;
  test_laplacian_init(&l, ei, n, PATH, cr, bfr);

  uint m = laplacian_size(l);

  slong out[2][1], wrk[2][1], in = m;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1, ng = out[1][0];

  scalar *u = tcalloc(scalar, m);
  scalar *v = tcalloc(scalar, m);

  int err = 0;
  for (slong i = 1; i <= ng; i++) {
    for (uint j = 0; j < m; j++) {
      if ((start + j) == i)
        u[j] = 1;
      else
        u[j] = 0;
    }

    laplacian_op(v, l, u);

    for (uint j = 0; j < m; j++) {
      slong li = start + j;
      if ((i > 1) && (i < ng)) {
        if (li == (i - 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
        if (li == i) err |= (fabs((v[j] - 2) / 2.0) > PARRSB_TEST_EPS);
        if (li == (i + 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
      }

      if (i == 1) {
        if (li == i) err |= (fabs(v[j] - 1) > PARRSB_TEST_EPS);
        if (li == (i + 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
      }

      if (i == ng) {
        if (li == (i - 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
        if (li == i) err |= (fabs(v[j] - 1) > PARRSB_TEST_EPS);
      }
    }
  }

  free(u), free(v);
  laplacian_free(&l);
  parrsb_assert(err == 0, c);
}

static int test_multiple_node_ring(const element_info ei, struct crystal *cr,
                                   buffer *bfr) {
  const struct comm *c = &(cr->comm);

  sint n = 3;

  laplacian l;
  test_laplacian_init(&l, ei, n, RING, cr, bfr);

  uint m = laplacian_size(l);

  slong out[2][1], wrk[2][1], in = m;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1, ng = out[1][0];

  scalar *u = tcalloc(scalar, m);
  scalar *v = tcalloc(scalar, m);

  int err = 0;
  for (slong i = 1; i <= ng; i++) {
    for (uint j = 0; j < m; j++) {
      if ((start + j) == i)
        u[j] = 1;
      else
        u[j] = 0;
    }

    laplacian_op(v, l, u);

    for (uint j = 0; j < m; j++) {
      slong li = start + j;
      if ((i > 1) && (i < ng)) {
        if (li == (i - 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
        if (li == i) err |= (fabs((v[j] - 2) / 2.0) > PARRSB_TEST_EPS);
        if (li == (i + 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
      }

      if (i == 1) {
        if (li == ng) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
        if (li == i) err |= (fabs((v[j] - 2) / 2.0) > PARRSB_TEST_EPS);
        if (li == (i + 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
      }

      if (i == ng) {
        if (li == (i - 1)) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
        if (li == i) err |= (fabs((v[j] - 2) / 2.0) > PARRSB_TEST_EPS);
        if (li == 1) err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
      }
    }
  }

  free(u), free(v);
  laplacian_free(&l);
  parrsb_assert(err == 0, c);
}

static int test_multiple_node_complete(const element_info ei,
                                       struct crystal *cr, buffer *bfr) {
  const struct comm *c = &(cr->comm);

  sint n = 2;

  laplacian l;
  test_laplacian_init(&l, ei, n, COMPLETE, cr, bfr);

  uint m = laplacian_size(l);

  slong out[2][1], wrk[2][1], in = m;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1, ng = out[1][0];

  scalar *u = tcalloc(scalar, m);
  scalar *v = tcalloc(scalar, m);

  int err = 0;
  for (slong i = 1; i <= ng; i++) {
    for (uint j = 0; j < m; j++) {
      if ((start + j) == i)
        u[j] = 1;
      else
        u[j] = 0;
    }

    laplacian_op(v, l, u);

    for (uint j = 0; j < m; j++) {
      slong li = start + j;
      if (li == i)
        err |= (fabs((v[j] - (ng - 1)) / (ng - 1)) > PARRSB_TEST_EPS);
      else
        err |= (fabs(v[j] + 1) > PARRSB_TEST_EPS);
    }
  }

  free(u), free(v);
  laplacian_free(&l);
  parrsb_assert(err == 0, c);
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
  err |= test_empty_path(ei, &cr, &bfr);
  err |= test_empty_ring(ei, &cr, &bfr);
  err |= test_empty_complete(ei, &cr, &bfr);
  err |= test_single_node_path(ei, &cr, &bfr);
  err |= test_single_node_ring(ei, &cr, &bfr);
  err |= test_single_node_complete(ei, &cr, &bfr);
  err |= test_multiple_node_path(ei, &cr, &bfr);
  err |= test_multiple_node_ring(ei, &cr, &bfr);
  err |= test_multiple_node_complete(ei, &cr, &bfr);

  element_info_free(&ei);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
