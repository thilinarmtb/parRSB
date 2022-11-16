#include <coarse.h>

// Find compressed/unique ids (cids) and the map (u2c) from user ids (ids) to
// cids. Also, return number of compress ids (cn).
static uint unique_ids(sint *u2c, ulong *cids, uint n, const ulong *ids,
                       buffer *bfr) {
  struct id_t {
    ulong id;
    uint idx;
    sint u2c;
  };

  struct array arr;
  array_init(struct id_t, &arr, n);

  struct id_t t = {.u2c = -1};
  for (uint i = 0; i < n; i++) {
    t.id = ids[i], t.idx = i;
    array_cat(struct id_t, &arr, &t, 1);
  }

  sarray_sort(struct id_t, arr.ptr, arr.n, id, 1, bfr);

  // Ignore the ids numbered zero
  sint cn = 0;
  ulong last = 0;
  struct id_t *pa = (struct id_t *)arr.ptr;
  for (uint i = 0; i < arr.n; i++) {
    if (pa[i].id != last)
      last = cids[cn] = pa[i].id, cn++;
    pa[i].u2c = cn - 1;
  }

  sarray_sort(struct id_t, pa, n, idx, 0, bfr);
  pa = (struct id_t *)arr.ptr;
  for (uint i = 0; i < n; i++)
    u2c[i] = pa[i].u2c;

  array_free(&arr);
  return cn;
}

struct array *assembly(uint n, const ulong *ids, uint nz, const uint *Ai,
                       const uint *Aj, const double *A, unsigned null_space,
                       unsigned type, const struct comm *c) {
  comm_barrier(c);
  double t0 = comm_time();

  // May be initialize buffer with better initial size.
  buffer bfr;
  buffer_init(&bfr, 1024);

  // TODO: Reorder dofs

  // TODO: Setup u2c and cid here.

  struct array *mijs = tcalloc(struct array, 1);

  double t1 = comm_time() - t0, min = t1, max = t1, wrk;
  comm_allreduce(c, gs_double, gs_max, &max, 1, &wrk);
  comm_allreduce(c, gs_double, gs_min, &min, 1, &wrk);
  if (c->id == 0) {
    printf("assembly: %g %g (min max)\n", min, max);
    fflush(stdout);
  }

  buffer_free(&bfr);

  return mijs;
}
