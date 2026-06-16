#include "metrics.h"
#include "parrsb_impl.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

static void mesh_load_balance(struct array *elist, uint nel,
                              const double *const xyz,
                              const long long *const vtx, const element_info ei,
                              struct crystal *cr, buffer *bfr) {
  struct comm *c = &cr->comm;

  slong out[2][1], wrk[2][1], in = nel;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];

  uint nstar = nelg / c->np;
  uint nrem = nelg - nstar * c->np;
  slong lower = (nstar + 1) * nrem;

  size_t unit_size = ei->size;
  int ndim = ei->nd;
  int nv = ei->nv;

  array_init_(elist, nel, unit_size, __FILE__, __LINE__);

  struct rcb_element *pe = (struct rcb_element *)calloc(1, unit_size);
  pe->origin = c->id;
  for (uint e = 0; e < nel; ++e) {
    slong eg = pe->globalId = start + e + 1;
    if (nstar == 0)
      pe->proc = eg - 1;
    else if (eg <= lower)
      pe->proc = (eg - 1) / (nstar + 1);
    else
      pe->proc = (eg - 1 - lower) / nstar + nrem;

    pe->coord[0] = pe->coord[1] = pe->coord[2] = 0.0;
    if (xyz != NULL) {
      for (int v = 0; v < nv; v++)
        for (int n = 0; n < ndim; n++)
          pe->coord[n] += xyz[e * ndim * nv + v * ndim + n];
      for (int n = 0; n < ndim; n++) pe->coord[n] /= nv;
    }

    array_cat_(unit_size, elist, pe, 1, __FILE__, __LINE__);
  }

  if (vtx != NULL) {
    struct rsb_element *pr = (struct rsb_element *)elist->ptr;
    for (uint e = 0; e < nel; e++)
      for (int v = 0; v < nv; v++) pr[e].vertices[v] = vtx[e * nv + v];
  }

  sarray_transfer_(elist, unit_size, offsetof(struct rcb_element, proc), 0, cr);
  if (vtx == NULL)
    sarray_sort(struct rcb_element, elist->ptr, elist->n, globalId, 1, bfr);
  else
    sarray_sort(struct rsb_element, elist->ptr, elist->n, globalId, 1, bfr);

  free(pe);
}

static void mesh_restore(int *part, struct array *elist, size_t usize,
                         struct crystal *cr, buffer *bfr) {
  sarray_transfer_(elist, usize, offsetof(struct rcb_element, origin), 1, cr);

  if (!part) goto free_array;

  uint nel = elist->n;
  if (usize == sizeof(struct rsb_element))
    sarray_sort(struct rsb_element, elist->ptr, nel, globalId, 1, bfr);
  else if (usize == sizeof(struct rcb_element))
    sarray_sort(struct rcb_element, elist->ptr, nel, globalId, 1, bfr);

  struct rcb_element *element = 0;
  for (uint e = 0; e < nel; e++) {
    element = (struct rcb_element *)((char *)elist->ptr + e * usize);
    part[e] = element->origin;
  }

free_array:
  array_free(elist);
}

static void parrsb_part_mesh_v0(int *part, const long long *const vtx,
                                const double *const xyz, const uint nel,
                                const unsigned nv,
                                const parrsb_options_t options,
                                const struct comm *const c,
                                struct crystal *const cr, buffer *const bfr) {
  element_info ei = tcalloc(struct element_info, 1);
  ei->nv = nv;
  ei->nd = nv_to_ndim(nv);
  ei->size = sizeof(struct rsb_element);
  ei->align = ALIGNOF(struct rsb_element);
  if (vtx == NULL) {
    ei->size = sizeof(struct rcb_element);
    ei->align = ALIGNOF(struct rcb_element);
  }

  struct array elist;
  mesh_load_balance(&elist, nel, xyz, vtx, ei, cr, bfr);

  struct comm active;
  comm_split(c, elist.n > 0, c->id, &active);

  const int verbose = options->verbose_level;

  // Setup topology aware partitioning.
  struct comm comms[8];
  mpi_topology(&options->levels, comms, &active, verbose);

  // Run the partitioner.
  parrsb_info(c, verbose, "parrsb_part_mesh_v0: running partitioner ...");
  if (elist.n > 0) {
    switch (options->partitioner) {
    case 0: rsb(&elist, ei, options, comms, bfr); break;
    case 1: rcb(&elist, ei, &active, bfr); break;
    case 2: rib(&elist, ei, &active, bfr); break;
    default: break;
    }
  }

  for (int l = 0; l < options->levels; l++) comm_free(&comms[l]);
  comm_free(&active);

  parrsb_info(c, verbose, "parrsb_part_mesh_v0: restore original input");
  mesh_restore(part, &elist, ei->size, cr, bfr);

  free(ei);
}

static void update_frontier(sint *const target, sint *const hop,
                            sint *const frontier, const unsigned nv,
                            const unsigned hid, buffer *const bfr) {
  // If target is already set, we don't update either target or hop.
  // We simply update frontier to previous target value and return.
  if (*target >= 0) {
    // Check invariant: *hop < INT_MAX
    assert(*hop < INT_MAX);
    for (uint i = 0; i < nv; i++) frontier[i] = *target;
    return;
  }

  struct dest_t {
    uint target;
  };

  struct array dests;
  array_init(struct dest_t, &dests, nv);
  {
    struct dest_t dt;
    for (uint i = 0; i < nv; i++) {
      if (frontier[i] >= 0) {
        dt.target = frontier[i];
        array_cat(struct dest_t, &dests, &dt, 1);
      }
    }
  }

  if (dests.n > 0) {
    sarray_sort(struct dest_t, dests.ptr, dests.n, target, 0, bfr);

    const struct dest_t *const pd = (const struct dest_t *const)dests.ptr;
    uint current_target = pd[0].target, current_count = 1;
    uint final_target = current_target, final_count = 1;
    for (uint i = 1; i < dests.n; i++) {
      if (pd[i].target == current_target) {
        current_count++;
      } else {
        if (current_count > final_count)
          final_count = current_count, final_target = current_target;
        current_target = pd[i].target, current_count = 1;
      }
    }
    if (current_count > final_count) final_target = current_target;

    // Update frontier, target and hop.
    for (uint j = 0; j < nv; j++) frontier[j] = final_target;
    *target = final_target, *hop = hid + 1;
  }

  array_free(&dests);
}

int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int ne, const int nv,
                     const parrsb_options_t options, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  options_update(options);
  // TODO: parrsb_options_validate(options, c);

  const int verbose = options->verbose_level;
  if (c.id == 0 && verbose) parrsb_options_print(options);

  if (vtx == NULL && xyz == NULL) {
    parrsb_info(&c, verbose,
                "parRSB: Both vertices and coordinates can't be NULL");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (vtx == NULL && options->partitioner == 0) {
    parrsb_info(&c, verbose,
                "parRSB: Vertices can't be NULL if the partitioner is RSB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (xyz == NULL && options->partitioner > 0) {
    parrsb_info(&c, verbose,
                "parRSB: Coordinates can't be NULL if the partitioner is "
                "RCB or RIB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (xyz == NULL && options->rsb_pre != 0) {
    parrsb_info(&c, verbose,
                "parRSB: Coordinates can't be NULL if the pre-partitioner is "
                "RCB or RIB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  slong ng = ne, wrk;
  comm_allreduce(&c, gs_long, gs_add, &ng, 1, &wrk);
  parrsb_info(&c, verbose, "parRSB: # elements = %lld, # vertices/element = %d",
              ng, nv);

  buffer bfr;
  buffer_init(&bfr, (ne + 1) * 72);

  struct crystal cr;
  crystal_init(&cr, &c);

  metric_init();

  parrsb_barrier(&c);
  const double t = comm_time();

  parrsb_part_mesh_v0(part, vtx, xyz, ne, nv, options, &c, &cr, &bfr);

  parrsb_info(&c, verbose, "par%s finished in %g seconds.",
              ALGO[options->partitioner], comm_time() - t);
  metric_rsb_print(&c, options->profile_level);

  metric_finalize();
  crystal_free(&cr);
  buffer_free(&bfr);
  comm_free(&c);

  return 0;
}

int graph_load_balance(struct array *nlist, uint nn, const long long *nodes,
                       const unsigned *offsets, const long long *neighbors,
                       struct crystal *cr, buffer *bfr) {
  const struct comm *c = &cr->comm;

  slong out[2][1], wrk[2][1], in = nn;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], ng = out[1][0];

  uint nstar = ng / c->np;
  uint nrem = ng - (slong)nstar * (slong)c->np;
  slong threshold = (nstar + 1) * nrem;

  array_init(struct graph_element, nlist, nn * 27);

  struct graph_element g;
  g.origin = c->id;
  for (uint i = 0; i < nn; i++) {
    slong ig = g.globalId = start + i;
    if (nstar == 0)
      g.proc = ig;
    else if (ig < threshold)
      g.proc = ig / (nstar + 1);
    else
      g.proc = nrem + (ig - threshold) / nstar;

    g.u = nodes[i];
    for (uint j = offsets[i], je = offsets[i + 1]; j < je; j++) {
      g.v = neighbors[j];
      array_cat(struct graph_element, nlist, &g, 1);
    }
    // Add a self edge.
    g.v = g.u;
    array_cat(struct graph_element, nlist, &g, 1);
  }

  sarray_transfer(struct graph_element, nlist, proc, 0, cr);
  sarray_sort_2(struct graph_element, nlist->ptr, nlist->n, u, 1, v, 1, bfr);

  return 0;
}

void element_info_init(element_info *ei_, element_type_t type) {
  element_info ei = *ei_ = tcalloc(struct element_info, 1);
  ei->nv = 0;
  ei->nd = 0;

  if (type == GRAPH) {
    ei->size = sizeof(struct graph_element);
    ei->align = ALIGNOF(struct graph_element);
  }
}

void element_info_free(element_info *ei) {
  if (!ei || !(*ei)) return;

  free(*ei), *ei = 0;
}

void graph_restore(int *part, struct array *nlist, struct crystal *cr,
                   buffer *bfr) {
  sarray_transfer(struct graph_element, nlist, origin, 1, cr);

  if (!part) goto free_array;

  uint n = nlist->n;
  struct graph_element *pg = (struct graph_element *)nlist->ptr;
  sarray_sort(struct graph_element, pg, n, globalId, 1, bfr);
  for (uint i = 0; i < n; i++) part[i] = pg[i].origin;

free_array:
  array_free(nlist);
}

int parrsb_part_graph(int *part, unsigned num_nodes, const long long *nodes,
                      const unsigned *offsets, const long long *neighbors,
                      const parrsb_options_t options, const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  options_update(options);
  int verbose = options->verbose_level;
  if (c.id == 0 && verbose) parrsb_options_print(options);

  slong n = num_nodes, wrk;
  comm_allreduce(&c, gs_long, gs_add, &n, 1, &wrk);
  parrsb_info(&c, verbose, "parRSB: # vertices = %lld\n", n);

  parrsb_barrier(&c);
  double t = comm_time();

  buffer bfr;
  buffer_init(&bfr, (num_nodes + 1) * 27);

  struct crystal cr;
  crystal_init(&cr, &c);

  metric_init();

  element_info ei;
  element_info_init(&ei, GRAPH);

  struct array nlist;
  graph_load_balance(&nlist, num_nodes, nodes, offsets, neighbors, &cr, &bfr);

  struct comm active;
  comm_split(&c, nlist.n > 0, c.id, &active);

  struct comm comms[8];
  mpi_topology(&options->levels, comms, &active, verbose);

  parrsb_info(&c, verbose, "parrsb_part_graph: running partitioner ...");
  if (nlist.n > 0) {
    switch (options->partitioner) {
    case 0: rsb(&nlist, ei, options, comms, &bfr); break;
    default: break;
    }
  }

  for (int i = 0; i < options->levels; i++) comm_free(&comms[i]);
  comm_free(&active);

  parrsb_info(&c, verbose, "parrsb_part_graph: restore original input");
  metric_rsb_print(&c, options->profile_level);

  graph_restore(part, &nlist, &cr, &bfr);
  element_info_free(&ei);
  metric_finalize();
  crystal_free(&cr);
  buffer_free(&bfr);

  t = comm_time() - t;
  parrsb_info(&c, verbose, "par%s finished in %g seconds.",
              ALGO[options->partitioner], t);

  comm_free(&c);

  return 0;
}
