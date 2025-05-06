#ifndef _PARRSB_TEST_H_
#define _PARRSB_TEST_H_

#include <math.h>

#include <parrsb_impl.h>

#define PARRSB_TEST_TOL 1e-12
#define PARRSB_TEST_EPS 1e-15

#define parrsb_test(cond, c)                                                   \
  do {                                                                         \
    sint err = !(cond);                                                        \
    sint wrk;                                                                  \
    comm_allreduce((c), gs_int, gs_add, &err, 1, &wrk);                        \
    return err;                                                                \
  } while (0)

struct graph_t {
  unsigned n;
  long long *nodes;
  unsigned *offsets;
  long long *neighbors;
};

typedef struct graph_t *graph_t;

/*
 * brick is a $N x 1 x 1$ graph or a mesh with N vertices or elements.
 * Here, $N * = \sum_{i=1}^{P}{n}$ and $P$ is the number of processes.
 */
static inline void graph_init_brick(graph_t *graph, unsigned n, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  slong out[2][1], wrk[2][1], in = n;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1, ng = out[1][0];

  graph_t g = *graph = tcalloc(struct graph_t, 1);
  g->n = n;

  g->nodes = tcalloc(long long, n);
  for (uint i = 0; i < n; i++) g->nodes[i] = start + i;

  g->offsets = tcalloc(unsigned, n + 1);
  for (uint i = 0; i < n; i++) {
    slong node = start + i;
    g->offsets[i + 1] = g->offsets[i] + (node > 1) + (node < ng);
  }

  g->neighbors = tcalloc(long long, g->offsets[n]);
  uint nn = 0;
  for (uint i = 0; i < n; i++) {
    slong node = start + i;
    if (node > 1) g->neighbors[nn++] = node - 1;
    if (node < ng) g->neighbors[nn++] = node + 1;
  }

  comm_free(&c);
}

static inline void graph_free(graph_t *graph) {
  if (!graph || !(*graph)) return;

  graph_t g = *graph;
  free(g->nodes), free(g->offsets), free(g->neighbors);
  free(*graph), *graph = 0;
}

static scalar norm2(const scalar *u, unsigned n) {
  scalar v = 0;
  for (uint i = 0; i < n; i++) v += u[i] * u[i];
  return sqrt(v);
}

#endif
