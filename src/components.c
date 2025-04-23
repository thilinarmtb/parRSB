#include "metrics.h"
#include "parrsb_impl.h"
#include "sort.h"

static slong mesh_components_v1(sint *component, struct array *elems, uint nv,
                                struct comm *c, buffer *buf) {
  slong out[2][1], wrk[2][1], in = elems->n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];
  if (nelg == 0) return 0;

  const uint nelt = elems->n;
  const uint nev = nelt * nv;
  const int null_input = (component == NULL);

  slong *p = tcalloc(slong, nev);
  slong *ids = tcalloc(slong, nev);
  if (null_input) component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++) component[e] = -1;

  struct unmarked {
    uint index;
  };
  struct unmarked u;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  struct comm cc;
  slong count = 0;
  slong nnz1, nnzg, nnzg0, nnzb;
  ulong nmarked = 0;
  do {
    // Count unmarked elements
    arr.n = 0;
    for (uint e = 0; e < nelt; e++) {
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
      for (uint e = 0; e < arr.n; e++)
        for (uint d = 0; d < nv; d++) p[e * nv + d] = 0;

      // Mark the first non-marked element as seed
      struct unmarked *ptr = (struct unmarked *)arr.ptr;
      slong first = start + ptr[0].index;
      slong mfirst = first;
      comm_allreduce(&cc, gs_long, gs_min, &mfirst, 1, wrk);

      if (mfirst == first) {
        for (uint d = 0; d < nv; d++) p[0 * nv + d] = 1;
      }

      // Setup gs
      for (uint e = 0; e < arr.n; e++)
        for (uint d = 0; d < nv; d++)
          ids[e * nv + d] = pe[ptr[e].index].vertices[d];

      struct gs_data *gsh = gs_setup(ids, arr.n * nv, &cc, 0, gs_pairwise, 0);

      do {
        gs(p, gs_long, gs_add, 0, gsh, buf);

        nnz1 = 0;
        for (uint e = 0; e < arr.n; e++) {
          uint d = 0;
          for (; d < nv; d++) {
            if (p[e * nv + d] > 0) {
              nnz1++;
              component[ptr[e].index] = count;
              break;
            }
          }
          // There was one non-zero vertex in the element
          if (d < nv) {
            for (d = 0; d < nv; d++) p[e * nv + d] = 1;
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

  array_free(&arr);

  free(p), free(ids);
  if (null_input == 1) free(component);

  return count;
}

struct cmp_t {
  uint c, p, uid;
};

static sint find_or_insert(struct array *cids, struct cmp_t *t) {
  // If there are no elements in the array, insert and exit.
  if (cids->n == 0) {
    array_cat(struct cmp_t, cids, t, 1);
    return -1;
  }

  // Otherwise, we will do a binary search.
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

  // Okay, not found -- insert at `mid` or `mid + 1`.
  uint max = cids->max;
  if (max == cids->n) {
    max += max / 2 + 1;
    pc = (struct cmp_t *)array_reserve(struct cmp_t, cids, max);
  }

  uint n = mid;
  if (t->c > pc[mid].c) n = mid + 1;

  struct cmp_t t0 = *t, t1;
  for (; n < cids->n; n++) {
    t1 = pc[n];
    pc[n] = t0;
    t0 = t1;
  }
  pc[n] = t0, cids->n++;

  // Sanity check.
  for (uint i = 1; i < cids->n; i++) assert(pc[i - 1].c < pc[i].c);

  return -1;
}

static slong mesh_components_v2(sint *component, const struct array *elems,
                                uint nv, const struct comm *ci, buffer *bfr) {
  slong out[2][1], wrk[2][1], in = elems->n;
  comm_scan(out, ci, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0];
  if (nelg == 0) return 0;

  metric_tic(ci, RSB_COMPONENTS);

  const uint nelt = elems->n;
  const uint nev = nelt * nv;
  const int null_input = (component == NULL);

  sint *p0 = tcalloc(sint, nev);
  sint *p = tcalloc(sint, nev);
  slong *ids = tcalloc(slong, nev);
  uint *inds = tcalloc(uint, nev);
  if (null_input) component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++) component[e] = -1;

  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  struct comm c;
  ulong nmkd = 0;
  slong nc = 0;
  do {
    // Copy unmarked elements to ids.
    uint unmkd = 0;
    for (uint e = 0; e < nelt; e++) {
      if (component[e] == -1) {
        inds[unmkd] = e;
        for (uint v = 0; v < nv; v++) ids[unmkd * nv + v] = pe[e].vertices[v];
        unmkd++;
      }
    }

    int bin = (unmkd > 0);
    comm_split(ci, bin, ci->id, &c);

    slong nnzg = 0, ncg = 0;
    if (bin == 1) {
      // Mark the first unmarked element as seed for the component c.id.
      for (uint v = 0; v < nv; v++) p[0 * nv + v] = c.id;

      // Initialize the rest of p.
      for (uint e = 1; e < unmkd; e++)
        for (uint v = 0; v < nv; v++) p[e * nv + v] = -1;

      // Setup gather-scatter to do BFS.
      struct gs_data *gsh = gs_setup(ids, unmkd * nv, &c, 0, gs_pairwise, 0);

      // Perform BFS.
      sint changed;
      do {
        for (uint i = 0; i < unmkd * nv; i++) p0[i] = p[i];

        gs(p, gs_int, gs_max, 0, gsh, bfr);

        changed = 0;
        sint nnz = 0;
        for (uint e = 0; e < unmkd; e++) {
          sint v0 = -1;
          for (uint v = 0; v < nv; v++) {
            if (p[e * nv + v] > -1) {
              if (v0 == -1)
                v0 = v;
              else if (p[e * nv + v0] < p[e * nv + v])
                v0 = v;
            }
          }

          // If there was at least one non-zero vertex in the element, we mark
          // the element with that value.
          if (v0 > -1) {
            sint c = p[e * nv + v0];
            for (uint v = 0; v < nv; v++) p[e * nv + v] = c;
            nnz++;
          }

          // Check if the component id changed.
          for (uint v = 0; v < nv; v++) {
            if (p[e * nv + v] != p0[e * nv + v]) {
              changed = 1;
              break;
            }
          }
        }

        nnzg = nnz;
        comm_allreduce(&c, gs_long, gs_add, &nnzg, 1, wrk);
        comm_allreduce(&c, gs_int, gs_add, &changed, 1, wrk);
      } while (changed);

      gs_free(gsh);

      // Find unique local components and then use them to find unique
      // global components.
      struct array cids;
      array_init(struct cmp_t, &cids, 100);

      struct cmp_t t;
      for (uint e = 0; e < unmkd; e++) {
        if (p[e * nv + 0] > -1) {
          t.c = p[e * nv + 0], t.p = t.c % c.np;
          find_or_insert(&cids, &t);
        }
      }

      struct crystal cr;
      crystal_init(&cr, &c);

      // Send the component id `C` to `C % P` where `P` is the number of
      // processors.
      sarray_transfer(struct cmp_t, &cids, p, 1, &cr);
      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);

      // Find unique components and number them globally.
      uint cnt = 0;
      if (cids.n > 0) {
        cnt++;
        struct cmp_t *pc = (struct cmp_t *)cids.ptr;
        for (uint i = 1; i < cids.n; i++) {
          if (pc[i].c > pc[i - 1].c) cnt++;
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
          if (pc[i].c > pc[i - 1].c) s++;
          pc[i].uid = s;
        }
      }

      sarray_transfer(struct cmp_t, &cids, p, 0, &cr);
      crystal_free(&cr);
      sarray_sort(struct cmp_t, cids.ptr, cids.n, c, 0, bfr);

      // Now assign the global component id to the marked elements.
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

  if (null_input == 1) free(component);
  free(p0), free(p), free(ids), free(inds);

  metric_toc(ci, RSB_COMPONENTS);
  return nc;
}

slong get_components(sint *component, const struct array *elems,
                     const element_info ei, const struct comm *c, buffer *buf) {
  if (ei->nv > 0)
    return mesh_components_v2(component, elems, ei->nv, c, buf);
  else
    return 0;
}

static sint repair_mesh(struct array *elements, uint nv, const struct comm *lc,
                        const struct comm *gc, int bin, buffer *bfr) {
  metric_tic(lc, RSB_BALANCE);

  struct ielem_t {
    uint index, orig;
    sint dest;
    scalar fiedler;
  };

  // Calculate expected # of elements per processor.
  uint ne = elements->n;
  slong nelgt = ne, nglob = ne, wrk[2][1];
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, wrk);
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, wrk);

  sint ne_ = nglob / gc->np, nrem = nglob - ne_ * gc->np;
  slong nelgt_exp = ne_ * lc->np + nrem / 2 + (nrem % 2) * (1 - bin);
  slong send_cnt = (nelgt - nelgt_exp) > 0 ? (nelgt - nelgt_exp) : 0;

  // Setup gather-scatter.
  uint size = ne * nv;
  slong *ids = tcalloc(slong, size);
  struct rsb_element *elems = (struct rsb_element *)elements->ptr;
  for (uint e = 0; e < ne; e++)
    for (uint v = 0; v < nv; v++) ids[e * nv + v] = elems[e].vertices[v];
  struct gs_data *gsh = gs_setup(ids, size, gc, 0, gs_pairwise, 0);

  sint *input = (sint *)ids;
  if (send_cnt > 0)
    for (uint e = 0; e < size; e++) input[e] = 0;
  else
    for (uint e = 0; e < size; e++) input[e] = 1;

  gs(input, gs_int, gs_add, 0, gsh, bfr);

  for (uint e = 0; e < ne; e++) elems[e].proc = gc->id;

  sint sid = (send_cnt == 0) ? gc->id : INT_MAX;
  comm_allreduce(gc, gs_int, gs_min, &sid, 1, &wrk);

  sint balanced = 0;
  if (send_cnt == 0) goto allreduce;

  struct array ielems;
  array_init(struct ielem_t, &ielems, 10);

  struct ielem_t ielem = {.orig = lc->id, .dest = -1};
  int mul = (sid == 0) ? 1 : -1;
  for (uint e = 0; e < ne; e++) {
    for (uint v = 0; v < nv; v++) {
      if (input[e * nv + v] > 0) {
        ielem.index = e, ielem.fiedler = mul * elems[e].fiedler;
        array_cat(struct ielem_t, &ielems, &ielem, 1);
        break;
      }
    }
  }

  // Sort based on fiedler value and sets `orig` field
  parallel_sort(struct ielem_t, &ielems, fiedler, gs_scalar, 0, 1, lc, bfr);

  slong out[2][1], nielems = ielems.n;
  comm_scan(out, lc, gs_long, gs_add, &nielems, 1, wrk);
  slong start = out[0][0];

  sint P = gc->np - lc->np;
  sint part_size = (send_cnt + P - 1) / P;

  if (out[1][0] >= send_cnt) {
    balanced = 1;
    struct ielem_t *ptr = ielems.ptr;
    for (uint e = 0; start + e < send_cnt && e < ielems.n; e++)
      ptr[e].dest = sid + (start + e) / part_size;

    struct crystal cr;
    crystal_init(&cr, lc);
    sarray_transfer(struct ielem_t, &ielems, orig, 0, &cr);
    crystal_free(&cr);

    ptr = ielems.ptr;
    for (uint e = 0; e < ielems.n; e++)
      if (ptr[e].dest != -1) elems[ptr[e].index].proc = ptr[e].dest;
  }

  array_free(&ielems);

allreduce:
  comm_allreduce(gc, gs_int, gs_max, &balanced, 1, &wrk);

  free(ids), gs_free(gsh);

  metric_toc(lc, RSB_BALANCE);
  return balanced;
}

static void repair(struct array *elements, const element_info ei,
                   const struct comm *lc, struct comm *gc, int bin,
                   buffer *bfr) {
  // Return if there is only one processor (or partition).
  if (gc->np == 1 || gc->np == lc->np || element_info_type(ei) == GRAPH) return;

  sint balanced = repair_mesh(elements, ei->nv, lc, gc, bin, bfr);
  if (!balanced) return;

  struct crystal cr;
  crystal_init(&cr, gc);
  sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
  crystal_free(&cr);

  parallel_sort_(elements, ei->size, ei->align, 0, 1, lc, bfr, 1, gs_scalar,
                 offsetof(struct base_element, fiedler));
}
