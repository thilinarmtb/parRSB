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

typedef enum { PATH = 0, RING, COMPLETE } graph_type_t;

/*
 * Path is a $N x 1 x 1$ graph or a mesh with N vertices or elements.
 * Here, $N = \sum_{i=1}^{P}{n}$ and $P$ is the number of processes.
 */
static inline void graph_init(graph_t *graph, unsigned n, graph_type_t type,
                              MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  slong out[2][1], wrk[2][1], in = n;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0] + 1, ng = out[1][0];

  graph_t g = *graph = tcalloc(struct graph_t, 1);
  g->n = n;

  if (ng == 0) return;

  g->nodes = tcalloc(long long, n);
  for (uint i = 0; i < n; i++) g->nodes[i] = start + i;

  g->offsets = tcalloc(unsigned, n + 1);
  if (type == PATH) {
    for (uint i = 0; i < n; i++) {
      slong node = start + i;
      g->offsets[i + 1] = g->offsets[i] + (node > 1) + (node < ng);
    }

    g->neighbors = tcalloc(long long, g->offsets[n]);
    uint count = 0;
    for (uint i = 0; i < n; i++) {
      slong node = start + i;
      if (node > 1) g->neighbors[count++] = node - 1;
      if (node < ng) g->neighbors[count++] = node + 1;
    }
  } else if (type == RING) {
    for (uint i = 0; i < n; i++) g->offsets[i + 1] = g->offsets[i] + 2;

    g->neighbors = tcalloc(long long, g->offsets[n]);
    uint count = 0;
    for (uint i = 0; i < n; i++) {
      slong node = start + i;
      g->neighbors[count++] = ((node - 1) == 0) ? ng : (node - 1);
      g->neighbors[count++] = ((node + 1) > ng) ? 1 : (node + 1);
    }
  } else if (type == COMPLETE) {
    for (uint i = 0; i < n; i++) g->offsets[i + 1] = g->offsets[i] + ng - 1;

    g->neighbors = tcalloc(long long, g->offsets[n]);
    slong count = 0;
    for (uint i = 0; i < n; i++) {
      ulong node = start + i;
      for (ulong j = 1; j < node; j++) g->neighbors[count++] = j;
      for (ulong j = node + 1; (slong)j <= ng; j++) g->neighbors[count++] = j;
    }
  }

  comm_free(&c);
}

static inline void graph_free(graph_t *graph) {
  if (!graph || !(*graph)) return;

  graph_t g = *graph;
  free(g->nodes), free(g->offsets), free(g->neighbors);
  free(*graph), *graph = 0;
}

static scalar norm2(const scalar *u, unsigned n, const struct comm *c) {
  scalar v = 0;
  for (uint i = 0; i < n; i++) v += u[i] * u[i];
  scalar wrk;
  comm_allreduce(c, gs_scalar, gs_add, &v, 1, &wrk);
  return sqrt(v);
}

#endif
