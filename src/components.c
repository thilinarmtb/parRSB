#include "genmap-impl.h"
#include "metrics.h"
#include "sort.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// Find the number of disconnected components
uint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv) {
  struct unmarked {
    uint index;
  };

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];

  if (nelg == 0)
    return 0;

  slong *p = tcalloc(slong, nelt * nv);
  slong *ids = tcalloc(slong, nelt * nv);

  int null_input = 0;
  if (null_input = (component == NULL))
    component = tcalloc(sint, nelt);

  uint e;
  for (e = 0; e < nelt; e++)
    component[e] = -1;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;
  struct unmarked u;
  uint d, count = 0;
  slong nnz1, nnzg, nnzg0, nnzb, nmarked = 0;
  do {
    // count unmarked elements
    arr.n = 0;
    for (e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        u.index = e;
        array_cat(struct unmarked, &arr, &u, 1);
      }
    }

    int bin = (arr.n > 0);
    comm_split(c, bin, c->id, &cc);

    nnz1 = nnzg = nnzg0 = 0;
    if (bin == 1) {
      // Initialize p
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          p[e * nv + d] = 0;

      // Mark the first non-marked element as seed
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index, mfirst = first;
      comm_allreduce(&cc, gs_long, gs_min, &mfirst, 1, wrk);

      if (mfirst == first) {
        for (d = 0; d < nv; d++)
          p[0 * nv + d] = 1;
      }

      // Setup gs
      for (e = 0; e < arr.n; e++)
        for (d = 0; d < nv; d++)
          ids[e * nv + d] = elements[ptr[e].index].vertices[d];

      struct gs_data *gsh = gs_setup(ids, arr.n * nv, &cc, 0, gs_pairwise, 0);

      do {
        gs(p, gs_long, gs_add, 0, gsh, buf);

        nnz1 = 0;
        for (e = 0; e < arr.n; e++) {
          for (d = 0; d < nv; d++) {
            if (p[e * nv + d] > 0) {
              nnz1++;
              component[ptr[e].index] = count;
              break;
            }
          }
          // There was one non-zero vertex in the element
          if (d < nv) {
            for (d = 0; d < nv; d++)
              p[e * nv + d] = 1;
          }
        }

        nnzg0 = nnzg, nnzg = nnz1;
        comm_allreduce(&cc, gs_long, gs_add, &nnzg, 1, &nnzb);
      } while (nnzg > nnzg0);
      gs_free(gsh);
    }
    comm_free(&cc);

    comm_allreduce(c, gs_long, gs_add, &nnz1, 1, &nnzb);
    nmarked += nnz1, count++;
  } while (nmarked < nelg);

  free(p), free(ids);
  if (null_input == 1)
    free(component);

  return count;
}

struct cmp_t {
  uint c, p, uid;
};

static sint find_or_insert(struct array *cids, struct cmp_t *t) {
  // If there are no elements in the array, insert and exit
  if (cids->n == 0) {
    array_cat(struct cmp_t, cids, t, 1);
    return -1;
  }

  // Otherwise, we will do a binary search
  struct cmp_t *pc = (struct cmp_t *)cids->ptr;
  sint s = 0, e = cids->n - 1, mid = 0;
  while (s <= e) {
    mid = (s + e) / 2;
    if (t->c == pc[mid].c)
      return pc[mid].uid;
    else if (t->c < pc[mid].c)
      e = mid - 1;
    else // t->c > pc[mid].c
      s = mid + 1;
  }

  // Okay, not found -- insert at `mid` or `mid + 1`
  uint max = cids->max;
  if (max == cids->n) {
    max += max / 2 + 1;
    pc = (struct cmp_t *)array_reserve(struct cmp_t, cids, max);
  }

  uint n = mid;
  if (t->c > pc[mid].c)
    n = mid + 1;

  struct cmp_t t0 = *t, t1;
  for (; n < cids->n; n++) {
    t1 = pc[n];
    pc[n] = t0;
    t0 = t1;
  }
  pc[n] = t0, cids->n++;

  return -1;
}

uint get_components_v2(sint *component, uint nelt, unsigned nv,
                       struct rsb_element *elements, const struct comm *ci,
                       buffer *bfr) {
  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, ci, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];

  if (nelg == 0)
    return 0;

  uint nev = nelt * nv;
  sint *p = tcalloc(sint, nev);
  slong *ids = tcalloc(slong, nev);
  uint *inds = tcalloc(uint, nev);

  int null_input = (component == NULL);
  if (null_input)
    component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++)
    component[e] = -1;

  struct comm c;
  slong nmkd = 0, nc = 0;
  do {
    // Copy unmarked elements to ids
    uint unmkd = 0;
    for (uint e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        assert(unmkd >= 0 && unmkd < nelt);
        inds[unmkd] = e;
        for (uint v = 0; v < nv; v++)
          ids[unmkd * nv + v] = elements[e].vertices[v];
        unmkd++;
      }
    }
    assert(unmkd <= nelt);

    int bin = (unmkd > 0);
    comm_split(ci, bin, ci->id, &c);

    slong nnzg = 0, nnzg0, ncg = 0;
    if (bin == 1) {
      // Setup gs
      struct gs_data *gsh = gs_setup(ids, unmkd * nv, &c, 0, gs_pairwise, 0);

      // Mark the first unmarked element as seed for the component c.id
      for (uint v = 0; v < nv; v++)
        p[0 * nv + v] = c.id;

      // Initialize the rest of p
      for (uint e = 1; e < unmkd; e++)
        for (uint v = 0; v < nv; v++)
          p[e * nv + v] = -1;

      do {
        nnzg0 = nnzg, nnzg = 0;

        gs(p, gs_int, gs_max, 0, gsh, bfr);

        for (uint e = 0; e < unmkd; e++) {
          uint v = 0;
          for (; v < nv; v++) {
            if (p[e * nv + v] > -1) {
              nnzg++;
              break;
            }
          }

          // There was one non-zero vertex in the element
          if (v < nv) {
            sint c = p[e * nv + v];
            for (v = 0; v < nv; v++)
              p[e * nv + v] = c;
          }
        }

        comm_allreduce(&c, gs_long, gs_add, &nnzg, 1, wrk);
      } while (nnzg0 < nnzg);
      gs_free(gsh);

      // Find unique local components and then use them to find unique
      // global components
      struct array cids;
      array_init(struct cmp_t, &cids, 100);

      struct cmp_t t;
      for (uint e = 0; e < unmkd; e++) {
        if (p[e * nv + 0] > -1) {
          t.c = p[e * nv + 0], t.p = t.c % c.np;
          find_or_insert(&cids, &t);
        }
      }
      assert(cids.n <= c.np);

      struct crystal cr;
      crystal_init(&cr, &c);
      sarray_transfer(struct cmp_t, &cids, p, 1, &cr);

      // find unique components and number them
      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);
      uint cnt = 0;
      if (cids.n > 0) {
        cnt++;
        struct cmp_t *pc = (struct cmp_t *)cids.ptr;
        for (uint i = 1; i < cids.n; i++) {
          if (pc[i].c > pc[i - 1].c)
            cnt++;
        }
      }

      in = cnt;
      comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
      ulong s = out[0][0];
      ncg = out[1][0];

      if (cids.n > 0) {
        struct cmp_t *pc = (struct cmp_t *)cids.ptr;
        pc[0].uid = s;
        for (uint i = 1; i < cids.n; i++) {
          if (pc[i].c > pc[i - 1].c)
            s++;
          pc[i].uid = s;
        }
      }

      sarray_transfer(struct cmp_t, &cids, p, 0, &cr);
      crystal_free(&cr);

      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);
      for (uint e = 0; e < unmkd; e++) {
        if (p[e * nv + 0] > -1) {
          t.c = p[e * nv + 0];
          sint uid = find_or_insert(&cids, &t);
          assert(uid > -1);
          component[inds[e]] = nc + uid;
        }
      }

      array_free(&cids);
    }
    comm_free(&c);

    comm_allreduce(ci, gs_long, gs_max, &nnzg, 1, &wrk);
    nmkd += nnzg;
    comm_allreduce(ci, gs_long, gs_max, &ncg, 1, &wrk);
    nc += ncg;
  } while (nmkd < nelg);

  free(p), free(ids), free(inds);
  if (null_input == 1)
    free(component);

  return nc;
}
