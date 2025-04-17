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

static sint balance_partitions_mesh(struct array *elements, uint nv,
                                    struct comm *lc, struct comm *gc, int bin,
                                    buffer *bfr) {
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
  parallel_sort(struct ielem_t, &ielems, fiedler, gs_double, 0, 1, lc, bfr);

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

static void balance_partitions(struct array *elements, const element_info ei,
                               struct comm *lc, struct comm *gc, int bin,
                               buffer *bfr) {
  // Return if there is only one processor (or partition).
  if (gc->np == 1 || gc->np == lc->np) return;

  sint balanced = 0;
  // Try to balance partitions if the input is a mesh.
  if (ei->nv > 0)
    balanced = balance_partitions_mesh(elements, ei->nv, lc, gc, bin, bfr);

  if (balanced && ei->nv > 0) {
    struct crystal cr;
    crystal_init(&cr, gc);
    sarray_transfer(struct rsb_element, elements, proc, 0, &cr);
    crystal_free(&cr);

    parallel_sort_(elements, ei->size, ei->align, 0, 1, lc, bfr, 1, gs_double,
                   offsetof(struct base_element, fiedler));
  } else {

    parallel_sort_(elements, ei->size, ei->align, 0, 1, gc, bfr, 1, gs_double,
                   offsetof(struct base_element, fiedler));
  }
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

static void pre_partitioner(struct array *elements, const element_info ei,
                            const struct comm *c, const parrsb_options options,
                            buffer *bfr) {
  switch (options->rsb_pre) {
  case 0:
    parallel_sort_(elements, ei->size, ei->align, 0, 1, c, bfr, 1, gs_long,
                   offsetof(struct base_element, globalId));
    break;
  case 1: rcb(elements, ei, c, bfr); break;
  case 2: rib(elements, ei, c, bfr); break;
  default: break;
  }
}

static void set_proc(struct array *elements, size_t esize,
                     const struct comm *c) {
  char *p = (char *)elements->ptr;
  for (uint i = 0; i < elements->n; i++)
    ((struct base_element *)(p + esize * i))->proc = c->id;
}

static void set_fiedler(struct array *elements, size_t esize, const scalar *f) {
  char *p = (char *)elements->ptr;
  for (uint i = 0; i < elements->n; i++)
    ((struct base_element *)(p + esize * i))->fiedler = f[i];
}

void rsb(struct array *elements, const element_info ei,
         const parrsb_options options, const struct comm *comms, buffer *bfr) {
  const struct comm *gc = &comms[0];

  scalar *f = tcalloc(scalar, elements->n + 1);
  const uint levels = options->levels;
  for (uint level = 0; level < levels; level++) {
    struct comm lc;
    comm_dup(&lc, &comms[level]);

    // Find the maximum number of RSB cuts in current level.
    uint ncuts = get_level_cuts(level, levels, comms);
    for (uint cut = 0; cut < ncuts; cut++) {
      pre_partitioner(elements, ei, &lc, options, bfr);

      set_proc(elements, ei->size, &lc);

      // Setup the laplacian and find the Fiedler vector.
      laplacian wl;
      laplacian_init(&wl, elements, ei, &lc, bfr);
      fiedler(f, wl, options, &lc, bfr);
      laplacian_free(&wl);

      // Sort the elements by Fiedler value.
      set_fiedler(elements, ei->size, f);
      parallel_sort_(elements, ei->size, ei->align, 0, 1, &lc, bfr, 1,
                     gs_double, offsetof(struct base_element, fiedler));

      // Get the bin of the current process and create a temporary communicator.
      struct comm tc;
      sint bin = get_bin(&lc, level, levels, comms);
      comm_split(&lc, bin, lc.id, &tc);

      // Find the number of disconnected components.
      if (options->find_disconnected_comps == 1) {
        slong nc = get_components(NULL, elements, ei, &tc, bfr);
        metric_acc(RSB_COMPONENTS_NCOMP, nc);
      }

      // Bisect and balance.
      balance_partitions(elements, ei, &tc, &lc, bin, bfr);

      // Split the communicator and recurse on the sub-problems.
      comm_free(&lc), comm_dup(&lc, &tc), comm_free(&tc);

      const uint nbrs = get_neighbors(elements, ei->nv, gc, &lc, bfr);
      metric_acc(RSB_NEIGHBORS, nbrs);

      metric_push_level();
    }
    comm_free(&lc);
  }

  free(f);

  check_rsb_partition(gc, options);
}
