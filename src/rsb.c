#include "metrics.h"
#include "parrsb_impl.h"
#include "sort.h"

/*
 * Find the partition id. A partition is a group of processors sharing the same
 * local communicator.
 */
static uint get_partition(const struct comm *const gc,
                          const struct comm *const lc) {
  sint out[2][1], wrk[2][1], root = (lc->id == 0);
  comm_scan(out, gc, gs_int, gs_add, &root, 1, wrk);
  sint part = out[0][0] * (lc->id == 0);
  comm_allreduce(lc, gs_int, gs_max, &part, 1, wrk);
  return part;
}

static uint get_neighbors(const struct array *const elems, const uint nv,
                          const struct comm *const gc,
                          const struct comm *const lc, buffer *bfr) {
  const uint n = elems->n;
  const uint size = elems->n * nv;

  struct vertex_t {
    ulong v;
    uint p, partition;
  };

  struct array vertices;
  array_init(struct vertex_t, &vertices, size);

  const struct rsb_element *const pe =
      (const struct rsb_element *const)elems->ptr;
  struct vertex_t vt = {.partition = get_partition(gc, lc)};
  for (uint i = 0; i < n; i++) {
    for (uint v = 0; v < nv; v++) {
      vt.v = pe[i].vertices[v], vt.p = vt.v % gc->np;
      array_cat(struct vertex_t, &vertices, &vt, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, gc);

  sarray_transfer(struct vertex_t, &vertices, p, 1, &cr);
  sarray_sort(struct vertex_t, vertices.ptr, vertices.n, v, 1, bfr);

  struct array neighbors;
  array_init(struct vertex_t, &neighbors, vertices.n * 27);

  const struct vertex_t *const pv = (const struct vertex_t *const)vertices.ptr;
  uint s = 0;
  while (s < vertices.n) {
    uint e = s + 1;
    while (e < vertices.n && pv[s].v == pv[e].v) e++;
    for (uint i = s; i < e; i++) {
      struct vertex_t vt = pv[i];
      for (uint j = s; j < e; j++) {
        vt.partition = pv[j].partition;
        array_cat(struct vertex_t, &neighbors, &vt, 1);
      }
    }
    s = e;
  }
  array_free(&vertices);

  sarray_transfer(struct vertex_t, &neighbors, p, 0, &cr);
  crystal_free(&cr);
  sarray_sort(struct vertex_t, neighbors.ptr, neighbors.n, partition, 0, bfr);

  // Now, extract out different partition ids found locally into an array.
  struct unique_t {
    uint p, partition;
  };

  struct array unique;
  array_init(struct unique_t, &unique, 27);

  if (neighbors.n > 0) {
    const struct vertex_t *const pn =
        (const struct vertex_t *const)neighbors.ptr;
    struct unique_t ut = {.partition = pn[0].partition,
                          .p = pn[0].partition % lc->np};
    array_cat(struct unique_t, &unique, &ut, 1);
    for (uint i = 1; i < neighbors.n; i++) {
      if (pn[i].partition > pn[i - 1].partition) {
        ut.partition = pn[i].partition, ut.p = ut.partition % lc->np;
        array_cat(struct unique_t, &unique, &ut, 1);
      }
    }
  }
  array_free(&neighbors);

  crystal_init(&cr, lc);
  sarray_transfer(struct unique_t, &unique, p, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct unique_t, unique.ptr, unique.n, partition, 0, bfr);
  sint un = 0;
  if (unique.n > 0) {
    un = 1;
    struct unique_t *pu = (struct unique_t *)unique.ptr;
    for (uint i = 1; i < unique.n; i++) {
      if (pu[i].partition > pu[un - 1].partition) pu[un] = pu[i], un++;
    }
  }
  array_free(&unique);

  sint wrk;
  comm_allreduce(lc, gs_int, gs_add, &un, 1, &wrk);
  assert(un >= 1);

  return un - 1;
}

/*
 * Find the number of disconnected components.
 */
static uint get_components(sint *component, struct array *elems, uint nv,
                           struct comm *c, buffer *buf) {
  uint nelt = elems->n;
  struct rsb_element *pe = (struct rsb_element *)elems->ptr;

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0], start = out[0][0];

  if (nelg == 0) return 0;

  uint nev = nelt * nv;
  slong *p = tcalloc(slong, nev);
  slong *ids = tcalloc(slong, nev);

  int null_input = (component == NULL);
  if (null_input) component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++) component[e] = -1;

  struct unmarked {
    uint index;
  };
  struct unmarked u;

  struct array arr;
  array_init(struct unmarked, &arr, nelt);

  struct comm cc;
  uint count = 0;
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

static slong get_components_v2(sint *component, struct array *elems, uint nv,
                               const struct comm *ci, buffer *bfr) {
  metric_tic(ci, RSB_COMPONENTS);

  slong nc = 0;

  slong out[2][1], wrk[2][1], in = elems->n;
  comm_scan(out, ci, gs_long, gs_add, &in, 1, wrk);
  ulong nelg = out[1][0];
  if (nelg == 0) goto exit_early;

  const uint nelt = elems->n;
  const uint nev = nelt * nv;
  sint *p0 = tcalloc(sint, nev);
  sint *p = tcalloc(sint, nev);
  slong *ids = tcalloc(slong, nev);
  uint *inds = tcalloc(uint, nev);

  int null_input = (component == NULL);
  if (null_input) component = tcalloc(sint, nelt);

  for (uint e = 0; e < nelt; e++) component[e] = -1;

  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  struct comm c;
  ulong nmkd = 0;
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

exit_early:
  metric_toc(ci, RSB_COMPONENTS);
  return nc;
}

static void check_rsb_partition(const struct comm *gc,
                                const parrsb_options opts) {
  sint max_levels = metric_get_levels();
  uint miter = opts->rsb_max_iter;
  uint mpass = opts->rsb_max_passes;

  slong wrk[4];
  comm_allreduce(gc, gs_int, gs_max, &max_levels, 1, (void *)wrk);

  for (sint i = 0; i < max_levels; i++) {
    sint converged = 1;
    uint val = (uint)metric_get_value(i, RSB_FIEDLER_CALC_NITER);
    if (opts->rsb_algo == 0) {
      if (val == miter * mpass) converged = 0;
    } else if (opts->rsb_algo == 1) {
      if (val == mpass) converged = 0;
    }

    struct comm c;
    comm_split(gc, converged, gc->id, &c);
    if (converged == 1) goto print_components;

    if (opts->rsb_algo == 0) {
      double init = metric_get_value(i, TOL_INIT);
      comm_allreduce(&c, gs_double, gs_min, &init, 1, (void *)wrk);

      double target = metric_get_value(i, TOL_TGT);
      comm_allreduce(&c, gs_double, gs_min, &target, 1, (void *)wrk);

      double final = metric_get_value(i, TOL_FNL);
      comm_allreduce(&c, gs_double, gs_min, &final, 1, (void *)wrk);

      if (c.id == 0) {
        fprintf(stderr,
                "Warning: Lanczos reached a residual of %lf (target: %lf) "
                "after %u x %u iterations in Level=%d!\n",
                final, target, mpass, miter, i);
        fflush(stderr);
      }
    } else if (opts->rsb_algo == 1) {
      if (c.id == 0) {
        fprintf(stderr,
                "Warning: Inverse iteration didn't converge after %d "
                "iterations in Level = %d\n",
                mpass, i);
        fflush(stderr);
      }
    }
    comm_free(&c);

  print_components:
    if (opts->find_disconnected_comps == 0) continue;

    slong minc = (slong)metric_get_value(i, RSB_COMPONENTS_NCOMP);
    slong maxc = minc;
    comm_allreduce(gc, gs_int, gs_min, &minc, 1, (void *)wrk);
    comm_allreduce(gc, gs_int, gs_max, &maxc, 1, (void *)wrk);

    if (maxc > 1 && gc->id == 0) {
      fprintf(stderr,
              "Warning: Partition created %lld/%lld (min/max) disconnected "
              "components in Level=%d!\n",
              minc, maxc, i);
      fflush(stderr);
    }
  }
}

static int balance_partitions(struct array *elements, uint nv, struct comm *lc,
                              struct comm *gc, int bin, buffer *bfr) {
  metric_tic(lc, RSB_BALANCE);

  // Return if there is only one processor (or partition).
  if (gc->np == 1 || gc->np == lc->np) goto early_exit;

  struct ielem_t {
    uint index, orig;
    sint dest;
    scalar fiedler;
  };

  // Calculate expected # of elements per processor.
  uint ne = elements->n;
  slong nelgt = ne, nglob = ne, wrk;
  comm_allreduce(lc, gs_long, gs_add, &nelgt, 1, &wrk);
  comm_allreduce(gc, gs_long, gs_add, &nglob, 1, &wrk);

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

  struct crystal cr;
  sint balanced = 0;

  if (send_cnt > 0) {
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
    parallel_sort(struct ielem_t, &ielems, fiedler, gs_double, 0, 1, lc, bfr);

    slong out[2][1], bfr[2][1], nielems = ielems.n;
    comm_scan(out, lc, gs_long, gs_add, &nielems, 1, bfr);
    slong start = out[0][0];

    sint P = gc->np - lc->np;
    sint part_size = (send_cnt + P - 1) / P;

    if (out[1][0] >= send_cnt) {
      balanced = 1;
      struct ielem_t *ptr = ielems.ptr;
      for (uint e = 0; start + e < send_cnt && e < ielems.n; e++)
        ptr[e].dest = sid + (start + e) / part_size;

      crystal_init(&cr, lc);
      sarray_transfer(struct ielem_t, &ielems, orig, 0, &cr);
      crystal_free(&cr);

      ptr = ielems.ptr;
      for (uint e = 0; e < ielems.n; e++)
        if (ptr[e].dest != -1) elems[ptr[e].index].proc = ptr[e].dest;
    }

    array_free(&ielems);
  }

  comm_allreduce(gc, gs_int, gs_max, &balanced, 1, &wrk);
  if (balanced == 1) {
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    // Do a load balanced sort in each partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, lc,
                  bfr);
  } else {
    // Forget about disconnected components, just do a load balanced partition
    parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, gc,
                  bfr);
  }

  free(ids), gs_free(gsh);
early_exit:
  metric_toc(lc, RSB_BALANCE);
  return 0;
}

static sint get_bin(const struct comm *const lc, const uint level,
                    const uint levels, const struct comm *comms) {
  sint psize = lc->np, pid = lc->id;
  if (level < levels - 1) {
    sint out[2][1], wrk[2][1], in = (comms[level + 1].id == 0);
    comm_scan(out, lc, gs_int, gs_add, &in, 1, wrk);
    psize = out[1][0], pid = (comms[level + 1].id == 0) * out[0][0];
    comm_allreduce(&comms[level + 1], gs_int, gs_max, &pid, 1, wrk);
  }

  return (pid >= (psize + 1) / 2);
}

static uint get_level_cuts(const uint level, const uint levels,
                           const struct comm *comms) {
  uint n = comms[level].np;
  if (level < levels - 1) {
    sint size = (comms[level + 1].id == 0), wrk;
    comm_allreduce(&comms[level], gs_int, gs_add, &size, 1, &wrk);
    n = size;
  }

  sint cuts = 0;
  uint pow2 = 1;
  while (pow2 < n) pow2 <<= 1, cuts++;

  sint wrk;
  comm_allreduce(&comms[0], gs_int, gs_max, &cuts, 1, &wrk);
  return cuts;
}

static void run_pre_partitioner(struct array *elements, int ndim,
                                const struct comm *lc,
                                const parrsb_options options, buffer *bfr) {
  metric_tic(lc, RSB_PRE);
  switch (options->rsb_pre) {
  case 0:
    parallel_sort(struct rsb_element, elements, globalId, gs_long, 0, 1, lc,
                  bfr);
    break;
  case 1: rcb(elements, sizeof(struct rsb_element), ndim, lc, bfr); break;
  case 2: rib(elements, sizeof(struct rsb_element), ndim, lc, bfr); break;
  default: break;
  }
  metric_toc(lc, RSB_PRE);
}

void rsb(struct array *elements, int nv, const parrsb_options options,
         const struct comm *comms, buffer *bfr) {
  const uint levels = options->levels;
  const uint ndim = nv_to_ndim(nv);
  const struct comm *gc = &comms[0];
  for (uint level = 0; level < levels; level++) {
    // Find the maximum number of RSB cuts in current level.
    uint ncuts = get_level_cuts(level, levels, comms);

    struct comm lc;
    comm_dup(&lc, &comms[level]);
    for (uint cut = 0; cut < ncuts; cut++) {
      run_pre_partitioner(elements, ndim, &lc, options, bfr);

      struct rsb_element *const pe = (struct rsb_element *const)elements->ptr;
      for (uint i = 0; i < elements->n; i++) pe[i].proc = lc.id;

      // Find the Fiedler vector.
      fiedler(elements, nv, options, &lc, bfr);

      // Sort by Fiedler value.
      metric_tic(&lc, RSB_SORT);
      parallel_sort(struct rsb_element, elements, fiedler, gs_double, 0, 1, &lc,
                    bfr);
      metric_toc(&lc, RSB_SORT);

      // Get the bin of the current process and create a temporary communicator.
      sint bin = get_bin(&lc, level, levels, comms);
      struct comm tc;
      comm_split(&lc, bin, lc.id, &tc);

      // Find the number of disconnected components.
      if (options->find_disconnected_comps == 1) {
        slong nc = get_components_v2(NULL, elements, nv, &tc, bfr);
        metric_acc(RSB_COMPONENTS_NCOMP, nc);
      }

      // Bisect and balance.
      balance_partitions(elements, nv, &tc, &lc, bin, bfr);

      // Split the communicator and recurse on the sub-problems.
      comm_free(&lc), comm_dup(&lc, &tc), comm_free(&tc);

      const uint nbrs = get_neighbors(elements, nv, gc, &lc, bfr);
      metric_acc(RSB_NEIGHBORS, nbrs);

      metric_push_level();
    }
    comm_free(&lc);
  }

  check_rsb_partition(gc, options);
}
