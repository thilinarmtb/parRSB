#include <parRSB.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-sort.h>

static int number_elements_aux(uint **level_off_, sint nlevels, uint nelt,
                               struct rsb_element *elems, struct comm *gc) {

  sint bfr;
  comm_allreduce(gc, gs_int, gs_max, &nlevels, 1, &bfr);

  slong *start = tcalloc(slong, nlevels);
  sint e;
  for (e = 0; e < nelt; e++)
    start[nlevels - 1 - elems[e].level]++;

  uint *level_off = *level_off_ = tcalloc(uint, nlevels + 1);
  level_off[0] = 0;
  for (e = 0; e < nlevels; e++)
    level_off[e + 1] = level_off[e] + start[e];

  slong *out = tcalloc(slong, 2 * nlevels);
  slong *buf = tcalloc(slong, 2 * nlevels);
  comm_scan(out, gc, gs_long, gs_add, start, nlevels, buf);

  for (e = nlevels; e < 2 * nlevels - 1; e++)
    out[e + 1] += out[e];

  for (e = 0; e < nelt; e++) {
    int lvl = nlevels - 1 - elems[e].level;
    elems[e].globalId = ++out[lvl] + out[nlevels + lvl];
  }

  free(start);
  free(out);
  free(buf);

  return nlevels;
}

static slong id_separator_levels(int level, uint nelt, int nv,
                                 struct rsb_element *elems, int bin,
                                 struct comm *lc, buffer *buf,
                                 struct comm *gc) {
  sint np = lc->np;
  sint id = lc->id;

  if (lc->np == 1)
    return 0;

  uint n = nelt * nv;
  slong *ids = tcalloc(slong, n);
  sint *interface = tcalloc(sint, n);

  uint e;
  int v;
  for (e = 0; e < nelt; e++)
    for (v = 0; v < nv; v++)
      ids[e * nv + v] = elems[e].vertices[v];

  struct gs_data *gsh = gs_setup(ids, n, lc, 0, gs_pairwise, 0);

  if (bin == 0) {
    for (e = 0; e < nelt; e++) {
      if (elems[e].level < 0)
        for (v = 0; v < nv; v++)
          interface[e * nv + v] = 1;
      else
        for (v = 0; v < nv; v++)
          interface[e * nv + v] = 0;
    }
  } else {
    for (e = 0; e < n; e++)
      interface[e] = 0;
  }

  gs(interface, gs_int, gs_add, 0, gsh, buf);

  if (bin == 1) {
    for (e = 0; e < nelt; e++) {
      if (elems[e].level < 0) { /* Not numbered before */
        for (v = 0; v < nv; v++)
          if (interface[e * nv + v] > 0) {
            elems[e].level = level;
            break;
          }
      } else { /* Previously numbered */
        for (v = 0; v < nv; v++)
          interface[e * nv + v] = 0;
      }
    }
  } else {
    for (e = 0; e < n; e++)
      interface[e] = 0;
  }

  gs(interface, gs_int, gs_add, 0, gsh, buf);

  if (bin == 0) {
    for (e = 0; e < nelt; e++) {
      if (elems[e].level < 0) { /* Not numbered before */
        for (v = 0; v < nv; v++)
          if (interface[e * nv + v] > 0) {
            elems[e].level = level;
            break;
          }
      }
    }
  }

  gs_free(gsh);
  free(ids);
  free(interface);
}

static int number_elements(struct comm **comms_, uint **level_off,
                           genmap_handle h) {
  int max_iter = 50;
  int max_pass = 50;

  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_number_faces_and_edges(h, gc);
  genmap_comm_scan(h, gc);

  struct rsb_element *e = genmap_get_elements(h);
  uint nelt = genmap_get_nel(h);

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int max_levels = ceil(log2(gc->np)) + 1;
  struct comm *comms = *comms_ = tcalloc(struct comm, max_levels);

  int level = 0;

  while (lc->np > 1) {
    /* Run RCB, RIB pre-step or just sort by global id */
    if (h->options->rsb_prepartition == 1) // RCB
      rcb(lc, h->elements, ndim, &h->buf);
    else if (h->options->rsb_prepartition == 2) // RIB
      rib(lc, h->elements, ndim, &h->buf);

    /* Run fiedler */
    int iter;
    int ipass = 0;
    do {
      genmap_vector ivec;
      genmap_init_vector(&ivec, ipass == 0, lc, h);

      GenmapFiedler(h, lc, max_iter, ivec);

      genmap_destroy_vector(ivec);
    } while (++ipass < max_pass && iter == max_iter);

    /* Sort by Fiedler vector */
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, bin_sort,
                  1, lc, &h->buf);

    /* Bisect */
    int bin = 1;
    if (lc->id < (lc->np + 1) / 2)
      bin = 0;
    repair_partitions(h, bin, level, lc, gc);

    nelt = genmap_get_nel(h);
    e = genmap_get_elements(h);
    id_separator_levels(level, nelt, nv, e, bin, lc, &h->buf, gc);

    comm_dup(&comms[level], lc);

    struct comm tc;
    genmap_comm_split(lc, bin, lc->id, &tc);
    comm_free(lc);
    comm_dup(lc, &tc);
    comm_free(&tc);

    genmap_comm_scan(h, lc);

    level++;
  }
  comm_dup(&comms[level], lc);

  e = genmap_get_elements(h);
  nelt = genmap_get_nel(h);
  uint i;
  for (i = 0; i < nelt; i++) {
    e[i].globalId = 0;
    if (e[i].level < 0)
      e[i].level = level;
  }

  return number_elements_aux(level_off, level + 1, nelt, e, gc);
}

struct csr_mat_ *ilu_laplacian(unsigned int *nlevels_, unsigned int **level_off,
                               struct comm **comms, genmap_handle h,
                               buffer *bfr) {
  *nlevels_ = number_elements(comms, level_off, h);
  uint nelt = genmap_get_nel(h);

  /* TODO get rid of gid */
  ulong *gid = tcalloc(ulong, nelt);
  struct rsb_element *e = genmap_get_elements(h);
  uint i;
  for (i = 0; i < nelt; i++)
    gid[i] = e[i].globalId;

  struct csr_mat_ *M = NULL;
  genmap_laplacian_csr_init(&M, gid, h, h->global);

  free(gid);

  return M;
}
