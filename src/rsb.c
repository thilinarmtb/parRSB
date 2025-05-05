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

element_t element_info_type(element_info ei) {
  if (ei->nv > 0) return MESH;
  return GRAPH;
}

static void prepartition(struct array *arr, const element_info ei,
                         const parrsb_options options, const struct comm *c,
                         buffer *bfr) {
  if (element_info_type(ei) != MESH) return;

  switch (options->rsb_pre) {
  case 0:
    parallel_sort_(arr, ei->size, ei->align, 0, 1, c, bfr, 1, gs_long,
                   offsetof(struct base_element, globalId));
    break;
  case 1: rcb(arr, ei, c, bfr); break;
  case 2: rib(arr, ei, c, bfr); break;
  default: break;
  }
}

static sint find_bin(const struct comm *const c, const uint level,
                     const uint levels, const struct comm *comms) {
  sint size = c->np, id = c->id;
  if (level < levels - 1) {
    sint out[2][1], wrk[2][1], in = (comms[level + 1].id == 0);
    comm_scan(out, c, gs_int, gs_add, &in, 1, wrk);
    size = out[1][0], id = (comms[level + 1].id == 0) * out[0][0];
    comm_allreduce(&comms[level + 1], gs_int, gs_max, &id, 1, wrk);
  }

  return (id >= (size + 1) / 2);
}

static void distribute_mesh(struct array *arr, const element_info ei,
                            const scalar *f, const struct comm *c,
                            buffer *bfr) {
  char *p = (char *)arr->ptr;
  for (uint i = 0; i < arr->n; i++)
    ((struct base_element *)(p + ei->size * i))->fiedler = f[i];

  parallel_sort_(arr, ei->size, ei->align, 0, 1, c, bfr, 1, gs_scalar,
                 offsetof(struct base_element, fiedler));
}

static void distribute_graph(struct array *arr, const element_info ei,
                             const scalar *f, const struct comm *c,
                             buffer *bfr) {
  sarray_sort_2(struct graph_element, arr->ptr, arr->n, u, 1, v, 1, bfr);
  graph_element pe = (const graph_element)arr->ptr;
  uint n = 0, nr = 0;
  while (n < arr->n) {
    uint i = n;
    while (i < arr->n && pe[i].u == pe[n].u) pe[i].fiedler = f[nr], i++;
    nr++, n = i;
  }

  // Sort based on the fiedler vector.
  parallel_sort_(arr, ei->size, ei->align, 0, 1, c, bfr, 1, gs_scalar,
                 offsetof(struct base_element, fiedler));

  // Setup a gs handle to find the smallest processor id which owns a
  // given row.
  buffer_reserve(bfr, arr->n * sizeof(slong));
  slong *ids = (slong *)bfr->ptr;

  pe = (const graph_element)arr->ptr;
  for (uint i = 0; i < arr->n; i++) ids[i] = pe[i].u;
  struct gs_data *gsh = gs_setup(ids, arr->n, c, 0, gs_crystal_router, 0);

  // Find the minimum processor id where a given row id resides.
  sint *proc = tcalloc(sint, arr->n);
  for (uint i = 0; i < arr->n; i++) proc[i] = c->id;
  gs(proc, gs_sint, gs_min, 0, gsh, bfr);
  for (uint i = 0; i < arr->n; i++) pe[i].proc = proc[i];
  free(proc), gs_free(gsh);

  // Send a given row to the minumum processor id.
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct graph_element, arr, proc, 0, &cr);

  // Re-calculate the number of rows and load balance.
  sarray_sort_2(struct graph_element, arr->ptr, arr->n, u, 1, v, 1, bfr);
  pe = (const graph_element)arr->ptr;
  n = nr = 0;
  while (n < arr->n) {
    uint i = n;
    while (i < arr->n && pe[i].u == pe[n].u) i++;
    nr++, n = i;
  }

  slong out[2][1], wrk[2][1], in = nr;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nrg = out[1][0];

  uint nstar = nrg / c->np;
  uint nrem = nrg - nstar * c->np;
  slong lower = (nstar + 1) * nrem;

  n = nr = 0;
  while (n < arr->n) {
    slong rg = start + nr;
    uint p = 0;
    if (nstar == 0)
      p = rg;
    else if (rg <= lower)
      p = rg / (nstar + 1);
    else
      p = (rg - lower) / nstar + nrem;

    uint i = n;
    while (i < arr->n && pe[i].u == pe[n].u) pe[i].proc = p, i++;
    n = i, nr++;
  }

  sarray_transfer(struct graph_element, arr, proc, 0, &cr);
  crystal_free(&cr);
}

static void bisect(struct comm *c, struct array *elements, const scalar *f,
                   sint bin, const element_info ei, buffer *bfr) {
  if (element_info_type(ei) == MESH) distribute_mesh(elements, ei, f, c, bfr);
  if (element_info_type(ei) == GRAPH) distribute_graph(elements, ei, f, c, bfr);

  // Split the communicator based on new partitions.
  struct comm tc;
  comm_split(c, bin, c->id, &tc);
  comm_free(c), comm_dup(c, &tc);
  comm_free(&tc);
}

static void calc_stats(const struct array *elements, const struct comm *gc,
                       const struct comm *lc, const parrsb_options options,
                       const element_info ei, buffer *bfr) {
  // Find the number of disconnected components.
  if (options->find_disconnected_comps == 1) {
    slong nc = get_components(NULL, elements, ei, lc, bfr);
    metric_acc(RSB_COMPONENTS_NCOMP, nc);
  }

  // Find the nymber of neighbors.
  uint nbrs = get_neighbors(elements, ei->nv, gc, lc, bfr);
  metric_acc(RSB_NEIGHBORS, nbrs);

  metric_push_level();
}

static void check_partition(const struct comm *gc, const parrsb_options opts) {
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

void rsb(struct array *elements, const element_info ei,
         const parrsb_options options, const struct comm *comms, buffer *bfr) {
  arena_t arena;
  arena_init(&arena);

  scalar *f = arena_tstart(scalar, arena, elements->n + 1);
  const struct comm *gc = &comms[0];
  const uint levels = options->levels;
  for (uint level = 0; level < levels; level++) {
    struct comm lc;
    comm_dup(&lc, &comms[level]);

    // Find the maximum number of RSB cuts in current level.
    uint ncuts = get_level_cuts(level, levels, comms);
    for (uint cut = 0; cut < ncuts; cut++) {
      // Pre-partition using RCB, RIB or simply by sorting. Only applicable for
      // mesh partitioning.
      prepartition(elements, ei, options, &lc, bfr);

      // Setup the laplacian and find the Fiedler vector.
      laplacian l;
      laplacian_init(&l, elements, ei, &lc);
      fiedler(f, l, options, &lc);
      laplacian_free(&l);

      // Bisect the elements by Fiedler value.
      sint bin = find_bin(&lc, level, levels, comms);
      bisect(&lc, elements, f, bin, ei, bfr);

      // Calculate partition stats.
      calc_stats(elements, gc, &lc, options, ei, bfr);
    }
    comm_free(&lc);
  }

  check_partition(gc, options);
  arena_stop(arena);
  arena_free(&arena);
}
