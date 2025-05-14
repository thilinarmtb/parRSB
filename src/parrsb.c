#include "metrics.h"
#include "parrsb_impl.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static char *ALGO[3] = {"RSB", "RCB", "RIB"};

static void options_update(parrsb_options options) {
#define OPT_UPDATE(opt, var, is_int)                                           \
  do {                                                                         \
    const char *val = getenv(var);                                             \
    if (val != NULL) {                                                         \
      if (is_int)                                                              \
        options->opt = atoi(val);                                              \
      else                                                                     \
        options->opt = atof(val);                                              \
    }                                                                          \
  } while (0)

  OPT_UPDATE(partitioner, "PARRSB_PARTITIONER", 1);
  OPT_UPDATE(tagged, "PARRSB_TAGGED", 1);
  OPT_UPDATE(levels, "PARRSB_LEVELS", 1);
  OPT_UPDATE(find_disconnected_comps, "PARRSB_FIND_DISCONNECTED_COMPONENTS", 1);
  OPT_UPDATE(repair, "PARRSB_REPAIR", 1);
  OPT_UPDATE(verbose_level, "PARRSB_VERBOSE_LEVEL", 1);
  OPT_UPDATE(profile_level, "PARRSB_PROFILE_LEVEL", 1);
  OPT_UPDATE(rsb_algo, "PARRSB_RSB_ALGO", 1);
  OPT_UPDATE(rsb_pre, "PARRSB_RSB_PRE", 1);
  OPT_UPDATE(rsb_max_iter, "PARRSB_RSB_MAX_ITER", 1);
  OPT_UPDATE(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", 1);
  OPT_UPDATE(rsb_tol, "PARRSB_RSB_TOL", 0);
  OPT_UPDATE(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", 1);
  OPT_UPDATE(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", 1);

#undef OPT_UPDATE
}

static void initialize_node_comm(struct comm *c, const struct comm *const gc) {
  MPI_Comm node;
  MPI_Comm_split_type(gc->c, MPI_COMM_TYPE_SHARED, gc->id, MPI_INFO_NULL,
                      &node);
  comm_init(c, node);
  MPI_Comm_free(&node);
}

static void initialize_levels(struct comm *const comms, parrsb_options options,
                              const struct comm *const c) {
  const int verbose = options->verbose_level;

  // Level 1 communicator is the global communicator.
  comm_dup(&comms[0], c);
  // Node level communicator is the last level communicator.
  struct comm nc;
  initialize_node_comm(&nc, c);

  // Find the number of nodes under the global communicator.
  sint in = (nc.id == 0), wrk;
  comm_allreduce(c, gs_int, gs_add, &in, 1, &wrk);
  uint nn = in;

  // Number of MPI ranks in the node level communicator.
  uint rpn = nc.np;

  // Check invariant: rpn should be the same across all the nodes.
  sint max = rpn, min = rpn;
  comm_allreduce(&comms[0], gs_int, gs_max, &max, 1, &wrk);
  comm_allreduce(&comms[0], gs_int, gs_min, &min, 1, &wrk);
  assert(max == min);

  // Check invariant: rpn must be larger than 0.
  assert(rpn > 0);
  parrsb_print(c, verbose, "parRSB: nodes = %u, ranks per node = %u", nn, rpn);

  // Hardcode the maximum number of levels to two for now.
  sint levels = MIN(2, options->levels);
  if (levels > 1) comm_dup(&comms[levels - 1], &nc);
  comm_free(&nc);

  parrsb_print(c, verbose, "parRSB: levels requested %u, enabled = %u",
               options->levels, levels);
  options->levels = levels;
}

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
                                const unsigned nv, const parrsb_options options,
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

  struct comm ca;
  comm_split(c, elist.n > 0, c->id, &ca);

  // Setup communicators for each level of the partitioning.
  struct comm comms[8];
  initialize_levels(comms, options, &ca);

  const int verbose = options->verbose_level;
  parrsb_print(c, verbose, "parrsb_part_mesh_v0: running partitioner ...");
  if (elist.n > 0) {
    switch (options->partitioner) {
    case 0: rsb(&elist, ei, options, comms, bfr); break;
    case 1: rcb(&elist, ei, &ca, bfr); break;
    case 2: rib(&elist, ei, &ca, bfr); break;
    default: break;
    }
  }

  for (int l = 0; l < options->levels; l++) comm_free(&comms[l]);
  comm_free(&ca);

  parrsb_print(c, verbose, "parrsb_part_mesh_v0: restore original input");
  mesh_restore(part, &elist, ei->size, cr, bfr);

  free(ei);
}

static void check_tagged_partitions(const long long *const eids,
                                    const long long *const vtx, const uint nel,
                                    const unsigned nv, const uint ntags,
                                    const struct comm *const c,
                                    const int verbose) {
  parrsb_print(c, verbose, "Check if the input elements are sorted locally.");
  {
    sint sorted = 1;
    for (uint i = 1; i < nel; i++) {
      if (eids[i] < eids[i - 1]) {
        sorted = 0;
        break;
      }
    }

    sint wrk;
    comm_allreduce(c, gs_int, gs_min, &sorted, 1, &wrk);
    if (!sorted) {
      if (c->id == 0) {
        fprintf(stderr, "Input elements are not sorted.\n");
        fflush(stderr);
      }
      exit(EXIT_FAILURE);
    }
  }

  // Number the elements within the each tag id and setup a gs handle based on
  // 2D element id.
  parrsb_print(c, verbose, "Number elements within each layer.");
  const uint tag_id = c->id / ntags;
  struct comm lc;
  struct gs_data *gse = NULL;
  {
    comm_split(c, tag_id, c->id, &lc);

    slong out[2][1], wrk[2][1], in = nel;
    comm_scan(out, &lc, gs_long, gs_add, &in, 1, wrk);
    slong start = out[0][0];

    slong *lids = tcalloc(slong, nel);
    for (uint i = 0; i < nel; i++) lids[i] = start + i;

    gse = gs_setup(lids, nel, c, 0, gs_pairwise, 0);
    free(lids);
  }

  // Setup a local gs handle based on the original gs vertex ids.
  parrsb_print(c, verbose, "Setup multiplicity.");
  const size_t size = nel * nv;
  buffer bfr;
  buffer_init(&bfr, size);
  sint *mul = tcalloc(sint, size);
  {
    struct gs_data *gsl = gs_setup(vtx, size, &lc, 0, gs_pairwise, 0);
    for (uint i = 0; i < size; i++) mul[i] = 1;
    gs(mul, gs_int, gs_add, 0, gsl, &bfr);
    gs_free(gsl);
  }

  // Now let's compare the multiplicity across the layers.
  parrsb_print(c, verbose, "Check multiplicity across the layers.");
  {
    sint *lmin = tcalloc(sint, nel);
    sint *lmax = tcalloc(sint, nel);
    for (uint v = 0; v < nv; v++) {
      for (uint e = 0; e < nel; e++) {
        lmin[e] = mul[e * nv + v];
        lmax[e] = mul[e * nv + v];
      }

      gs(lmin, gs_int, gs_min, 0, gse, &bfr);
      gs(lmax, gs_int, gs_max, 0, gse, &bfr);

      for (uint e = 0; e < nel; e++) assert(lmin[e] == lmax[e]);
    }

    free(lmin), free(lmax);
  }

  free(mul);
  buffer_free(&bfr);
  gs_free(gse);
  comm_free(&lc);

  return;
}

static void parrsb_part_mesh_v1(int *part, const long long *const vtx,
                                const double *const xyz, const int *const tag,
                                const uint nel, const unsigned nv,
                                const parrsb_options options,
                                const struct comm *const c,
                                struct crystal *const cr, buffer *const bfr) {
  const int verbose = options->verbose_level;
  parrsb_print(c, verbose, "Find number of tags in the mesh ...");

  struct tag_t {
    uint p, tag, seq, tagn;
  };

  struct array tags;
  array_init(struct tag_t, &tags, nel);

  {
    struct tag_t tt;
    for (uint i = 0; i < nel; i++) {
      tt.seq = i, tt.tag = tag[i], tt.p = tt.tag % c->np;
      array_cat(struct tag_t, &tags, &tt, 1);
    }
    sarray_sort(struct tag_t, tags.ptr, tags.n, tag, 0, bfr);
  }

  struct array unique;
  array_init(struct tag_t, &unique, 1024);

  if (tags.n > 0) {
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    array_cat(struct tag_t, &unique, &pt[0], 1);
    for (uint i = 1; i < tags.n; i++) {
      if (pt[i].tag > pt[i - 1].tag)
        array_cat(struct tag_t, &unique, &pt[i], 1);
    }
  }

  sint out[2][1];
  {
    sarray_transfer(struct tag_t, &unique, p, 1, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);

    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    sint in = 0;
    if (unique.n > 0) {
      in = 1;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag) in++;
      }
    }

    sint wrk[2][1];
    comm_scan(out, c, gs_int, gs_add, &in, 1, wrk);
  }
  const uint num_tags = out[1][0], tag_start = out[0][0];

  parrsb_print(c, verbose, "Num tags: %d", num_tags);
  if (c->np % num_tags != 0) {
    if (c->id == 0) {
      fprintf(stderr,
              "Number of processes must be a multiple of number of tags: "
              "processes = %d, tags = %d.\n",
              c->np, num_tags);
    }
    exit(EXIT_FAILURE);
  }

  {
    struct tag_t *const pu = (struct tag_t *const)unique.ptr;
    uint start = tag_start;
    if (unique.n > 0) {
      pu[0].tagn = start;
      for (uint i = 1; i < unique.n; i++) {
        if (pu[i].tag > pu[i - 1].tag) start++;
        pu[i].tagn = start;
      }
    }

    sarray_transfer(struct tag_t, &unique, p, 0, cr);
    sarray_sort(struct tag_t, unique.ptr, unique.n, tag, 0, bfr);
  }

  const uint chunk_size = c->np / num_tags;
  parrsb_print(c, verbose, "Processes per tag: %d", chunk_size);
  {
    struct tag_t *const pt = (struct tag_t *const)tags.ptr;
    const struct tag_t *const pu = (const struct tag_t *const)unique.ptr;
    for (uint i = 0, s = 0; i < unique.n; i++) {
      uint e = s + 1;
      assert(pt[s].tag == pu[i].tag);
      while (e < tags.n && pt[e].tag == pu[i].tag) e++;
      for (uint j = s; j < e; j++)
        pt[j].p = chunk_size * pu[i].tagn + pt[i].seq % chunk_size;
      s = e;
    }

    sarray_sort(struct tag_t, tags.ptr, tags.n, seq, 0, bfr);
  }
  array_free(&unique);

  struct element_t {
    uint proc, part, seq;
    scalar xyz[MAXDIM * MAXNV];
    slong vertices[MAXNV];
  };

  struct array elements;
  array_init(struct element_t, &elements, nel);

  parrsb_print(c, verbose,
               "Pack element data for transfering. tags.n=%u, nel=%u", tags.n,
               nel);
  const unsigned ndim = (nv == 8) ? 3 : 2;
  {
    assert(tags.n == nel);
    const struct tag_t *const pt = (const struct tag_t *const)tags.ptr;
    struct element_t et;
    for (uint i = 0; i < tags.n; i++) {
      et.proc = pt[i].p, et.seq = i;
      for (uint j = 0; j < nv; j++) {
        et.vertices[j] = vtx[i * nv + j];
        for (uint k = 0; k < ndim; k++)
          et.xyz[j * ndim + k] = xyz[i * nv * ndim + j * ndim + k];
      }
      array_cat(struct element_t, &elements, &et, 1);
    }

    sarray_transfer(struct element_t, &elements, proc, 1, cr);
  }
  array_free(&tags);

  parrsb_print(c, verbose, "Copy element data for feeding to parRSB.");
  long long *lvtx = tcalloc(long long, (elements.n + 1) * nv);
  double *lxyz = tcalloc(double, (elements.n + 1) * nv * ndim);
  {
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      for (uint j = 0; j < nv; j++) {
        lvtx[e * nv + j] = pe[e].vertices[j];
        for (uint k = 0; k < ndim; k++)
          lxyz[e * nv * ndim + j * ndim + k] = pe[e].xyz[j * ndim + k];
      }
    }
  }

  parrsb_print(c, verbose, "Run parRSB locally within a tag now.");
  {
    int *lpart = tcalloc(int, elements.n + 1);

    struct comm lc;
    comm_split(c, c->id / chunk_size, c->id, &lc);

    struct crystal lcr;
    crystal_init(&lcr, &lc);

    options->verbose_level = 0;
    options->profile_level = 0;
    parrsb_part_mesh_v0(lpart, lvtx, lxyz, elements.n, nv, options, &lc, &lcr,
                        bfr);
    crystal_free(&lcr), comm_free(&lc);

    struct element_t *const pe = (struct element_t *const)elements.ptr;
    for (uint e = 0; e < elements.n; e++) {
      pe[e].part = lpart[e] + (c->id / chunk_size) * chunk_size;
      assert(pe[e].part < c->np);
    }
    free(lpart);

    sarray_transfer(struct element_t, &elements, proc, 0, cr);
    assert(nel == elements.n);
  }
  free(lvtx), free(lxyz);

  {
    sarray_sort(struct element_t, elements.ptr, elements.n, seq, 0, bfr);
    const struct element_t *const pe =
        (const struct element_t *const)elements.ptr;
    for (uint i = 0; i < nel; i++) part[i] = pe[i].part;
  }

  array_free(&elements);
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

void parrsb_part_solid(int *part, const long long *const vtx2,
                       const unsigned nel2, const long long *const vtx1,
                       const unsigned nel1, const unsigned nv,
                       const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);
  parrsb_print(&c, 1, "Running greedy solid ... nel1 = %d nel2 = %d", nel1,
               nel2);

  for (uint i = 0; i < nel2; i++) part[i] = -1;

  buffer bfr;
  buffer_init(&bfr, 1024);

  struct crystal cr;
  crystal_init(&cr, &c);

  // Return if global size is 0.
  const uint nelt = nel1 + nel2;
  slong nelg = nelt;
  {
    slong wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelg, 1, &wrk);
    if (nelg == 0) {
      parrsb_print(&c, 1, "Mesh is empty ...");
      crystal_free(&cr);
      buffer_free(&bfr);
      comm_free(&c);
      return;
    }
  }

  const size_t size1 = nel1 * nv;
  const size_t size2 = nel2 * nv;
  const size_t size = size1 + size2;

  // Setup the gather-scatter handle to find connectivity through BFS.
  parrsb_print(&c, 1, "Setup gather-scatter handle ...");
  struct gs_data *gsh = NULL;
  {
    slong *vtx = tcalloc(slong, size);
    for (size_t i = 0; i < size1; i++) vtx[i] = vtx1[i];
    for (size_t i = 0; i < size2; i++) vtx[size1 + i] = vtx2[i];

    gsh = gs_setup(vtx, size, &c, 0, gs_pairwise, 0);
    free(vtx);
  }

  // Check if the solid + fluid mesh is connected. Otherwise, we cannot use
  // the greedy solid partitioner.
  parrsb_print(&c, 1, "Check if fluid + solid is connected ...");
  {
    slong wrk;
    sint idmin = (c.id + 1) * (size > 0);
    comm_allreduce(&c, gs_int, gs_min, &idmin, 1, &wrk);
    assert(idmin > 0);

    sint *const component = tcalloc(sint, size);
    if (c.id + 1 == (uint)idmin) {
      for (uint i = 0; i < nv; i++) component[i] = 1;
    }

    slong marked0 = 0, marked1 = 1;
    sint epoch = 0;
    while (marked1 > marked0) {
      gs(component, gs_int, gs_max, 0, gsh, &bfr);

      marked0 = marked1, marked1 = 0;
      for (uint i = 0; i < nel1 + nel2; i++) {
        sint v = 0;
        for (uint j = 0; j < nv; j++) v += component[i * nv + j];
        if (v > 0) {
          for (uint j = 0; j < nv; j++) component[i * nv + j] = 1;
          marked1 += 1;
        }
      }

      comm_allreduce(&c, gs_long, gs_add, &marked1, 1, &wrk);
      parrsb_print(&c, 1, "\tepoch = %d marked0 = %lld marked1 = %lld", epoch,
                   marked0, marked1);
      epoch++;
    }
    free(component);

    if (marked1 != nelg) {
      if (c.id == 0) {
        fprintf(stderr, "Fluid + Solid mesh is not connected.\n");
        fflush(stderr);
      }
      exit(EXIT_FAILURE);
    }
  }

  // Calculate the global number of elements in solid mesh and expected number
  // of elements in each partition.
  parrsb_print(&c, 1, "Calculate expected number of elements ...");
  slong nelgt2 = nel2;
  uint nexp2;
  {
    slong wrk;
    comm_allreduce(&c, gs_long, gs_add, &nelgt2, 1, &wrk);
    nexp2 = nelgt2 / c.np;
    nexp2 += (c.id < (nelgt2 - nexp2 * c.np));
    // Check for invariant: (min(nexp2) -  max(nexp2)) <= 1.
    slong nexp2_min = nexp2, nexp2_max = nexp2;
    comm_allreduce(&c, gs_long, gs_min, &nexp2_min, 1, &wrk);
    comm_allreduce(&c, gs_long, gs_max, &nexp2_max, 1, &wrk);
    assert(nexp2_max - nexp2_min <= 1);
    // Check for invariant: (sum(nexp2) == nelgt2).
    slong nexp2_sum = nexp2;
    comm_allreduce(&c, gs_long, gs_add, &nexp2_sum, 1, &wrk);
    assert(nexp2_sum == nelgt2);
  }

  // Initialize array of elements to be sent to each partition.
  struct elem_t {
    sint part;
    uint target, hop, sequence;
  };

  struct array arr;
  array_init(struct elem_t, &arr, nel2);

  // Allocate space for work arrays: frontier, target, and hop.
  sint *const frontier = tcalloc(sint, size);
  sint *const target = tcalloc(sint, nelt);
  sint *const hop = tcalloc(sint, nelt);

  uint nrecv2 = 0;
  slong nrem2 = nelgt2;
  while (nrem2 > 0) {
    parrsb_print(&c, 1, "nrem2 = %lld", nrem2);

    // Check for invariant: nrecv2 <= nexp2.
    assert(nrecv2 <= nexp2);

    // If the partition does not have enough elements, we keep it under
    // consideration for accepting new solid elements. If the partition
    // already has enough elements, we take that partition out of
    // consideration (by setting the frontier to -1). We always initialize solid
    // elements as unassigned (-1) although they may be already assigned. We
    // check for that later when we actually assign the elements to partitions.
    {
      sint id = c.id, hid = 0;
      if (nrecv2 == nexp2) id = -1, hid = INT_MAX;

      // Max id should be >= 0;
      sint wrk, idmax = id;
      comm_allreduce(&c, gs_int, gs_max, &idmax, 1, &wrk);
      assert(idmax >= 0);

      // Initialize frontier, target, and hop.
      for (uint i = 0; i < size1; i++) frontier[i] = id;
      for (uint i = size1; i < size; i++) frontier[i] = -1;
      for (uint i = 0; i < nel1; i++) target[i] = id, hop[i] = hid;
      for (uint i = nel1; i < nelt; i++) target[i] = -1, hop[i] = INT_MAX;
    }

    // Then perform a BFS till we assign all the elements in the solid mesh with
    // a potential partition id.
    parrsb_print(&c, 1, "Assign partition id ...");
    {
      sint assigned = 0;
      slong wrk;
      for (uint hid = 0; !assigned; hid++) {
        gs(frontier, gs_int, gs_max, 0, gsh, &bfr);

        assigned = 1;
        slong unassigned = 0;
        for (uint i = 0; i < nelt; i++) {
          update_frontier(&target[i], &hop[i], &frontier[i * nv], nv, hid,
                          &bfr);
          assigned = assigned && (target[i] >= 0);
          unassigned += (target[i] < 0);
        }

        comm_allreduce(&c, gs_int, gs_min, &assigned, 1, &wrk);
        comm_allreduce(&c, gs_long, gs_add, &unassigned, 1, &wrk);
        parrsb_print(&c, 1, "hid = %d, assigned = %d unassigned = %d", hid,
                     assigned, unassigned);
      }
    }

    // Pack unassigned solid elements and send them to the target partition.
    arr.n = 0;
    {
      struct elem_t et = {.part = -1};
      for (uint i = 0; i < nel2; i++) {
        if (part[i] >= 0) continue;
        et.sequence = i, et.target = target[nel1 + i], et.hop = hop[nel1 + i];
        array_cat(struct elem_t, &arr, &et, 1);
      }

      parrsb_print(&c, 1, "Send elemenets to the target partition ...");
      sarray_transfer(struct elem_t, &arr, target, 1, &cr);
    }

    // Assign elements if the partition still doesn't have enough elements.
    if (nrecv2 < nexp2) {
      // We sort by hop value. Elements with lower hop value are assigned first
      // since they are technically closer to the partition.
      sarray_sort(struct elem_t, arr.ptr, arr.n, hop, 1, &bfr);
      struct elem_t *const pa = (struct elem_t *const)arr.ptr;
      uint keep = MIN(nexp2 - nrecv2, arr.n);
      for (uint i = 0; i < keep; i++) pa[i].part = c.id;
      nrecv2 += keep;
      // Check for invariant: nrecv2 <= nexp2.
      assert(nrecv2 <= nexp2);
    }

    // Send everything back with updated partition id and update the part array.
    {
      parrsb_print(&c, 1, "Send everything back ...");
      sarray_transfer(struct elem_t, &arr, target, 0, &cr);

      const struct elem_t *const pa = (const struct elem_t *const)arr.ptr;
      for (uint j = 0; j < arr.n; j++) part[pa[j].sequence] = pa[j].part;
      arr.n = 0;
    }

    // Update the number of elements remaining.
    {
      slong wrk;
      nrem2 = nexp2 - nrecv2;
      comm_allreduce(&c, gs_long, gs_add, &nrem2, 1, &wrk);
    }
  }

  gs_free(gsh);
  free(frontier), free(target), free(hop);
  array_free(&arr);
  crystal_free(&cr);
  buffer_free(&bfr);
  comm_free(&c);
}

int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int *const tag,
                     const int ne, const int nv, const parrsb_options options,
                     MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  options_update(options);
  // TODO: parrsb_options_validate(options, c);

  const int verbose = options->verbose_level;
  if (c.id == 0 && verbose) parrsb_options_print(options);

  if (vtx == NULL && xyz == NULL) {
    parrsb_print(&c, verbose,
                 "parRSB: Both vertices and coordinates can't be NULL");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (vtx == NULL && options->partitioner == 0) {
    parrsb_print(&c, verbose,
                 "parRSB: Vertices can't be NULL if the partitioner is RSB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (xyz == NULL && options->partitioner > 0) {
    parrsb_print(&c, verbose,
                 "parRSB: Coordinates can't be NULL if the partitioner is "
                 "RCB or RIB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (xyz == NULL && options->rsb_pre != 0) {
    parrsb_print(&c, verbose,
                 "parRSB: Coordinates can't be NULL if the pre-partitioner is "
                 "RCB or RIB");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  if (options->tagged == 1 && !tag) {
    parrsb_print(
        &c, verbose,
        "parRSB: Tagged partitioning requested but the tag array is NULL");
    MPI_Abort(c.c, EXIT_FAILURE);
  }

  slong ng = ne, wrk;
  comm_allreduce(&c, gs_long, gs_add, &ng, 1, &wrk);
  parrsb_print(&c, verbose,
               "parRSB: # elements = %lld, # vertices/element = %d", ng, nv);

  buffer bfr;
  buffer_init(&bfr, (ne + 1) * 72);

  struct crystal cr;
  crystal_init(&cr, &c);

  metric_init();

  parrsb_barrier(&c);
  const double t = comm_time();

  if (options->tagged == 1)
    parrsb_part_mesh_v1(part, vtx, xyz, tag, ne, nv, options, &c, &cr, &bfr);
  if (options->tagged == 0)
    parrsb_part_mesh_v0(part, vtx, xyz, ne, nv, options, &c, &cr, &bfr);

  parrsb_print(&c, verbose, "par%s finished in %g seconds.",
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

void element_info_init(element_info *ei_, element_t type) {
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
                      const parrsb_options options, const MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  options_update(options);
  int verbose = options->verbose_level;
  if (c.id == 0 && verbose) parrsb_options_print(options);

  slong n = num_nodes, wrk;
  comm_allreduce(&c, gs_long, gs_add, &n, 1, &wrk);
  parrsb_print(&c, verbose, "parRSB: # vertices = %lld\n", n);

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

  struct comm ca;
  comm_split(&c, nlist.n > 0, c.id, &ca);

  struct comm comms[8];
  initialize_levels(comms, options, &ca);

  parrsb_print(&c, verbose, "parrsb_part_graph: running partitioner ...");
  if (nlist.n > 0) {
    switch (options->partitioner) {
    case 0: rsb(&nlist, ei, options, comms, &bfr); break;
    default: break;
    }
  }

  for (int i = 0; i < options->levels; i++) comm_free(&comms[i]);
  comm_free(&ca);

  parrsb_print(&c, verbose, "parrsb_part_graph: restore original input");
  metric_rsb_print(&c, options->profile_level);

  graph_restore(part, &nlist, &cr, &bfr);
  element_info_free(&ei);
  metric_finalize();
  crystal_free(&cr);
  buffer_free(&bfr);

  t = comm_time() - t;
  parrsb_print(&c, verbose, "par%s finished in %g seconds.",
               ALGO[options->partitioner], t);

  comm_free(&c);

  return 0;
}
