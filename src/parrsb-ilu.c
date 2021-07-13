#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-sort.h>

#include <parRSB.h>

/*
 * Numbering
 */
static int number_elements_aux(MPI_Comm **comms_, uint **level_off_,
                               sint nlevel, uint nelt,
                               struct rsb_element *elems, struct comm *gc) {
  sint bfr;
  comm_allreduce(gc, gs_int, gs_max, &nlevel, 1, &bfr);

  slong *start = tcalloc(slong, nlevel);
  sint e;
  for (e = 0; e < nelt; e++)
    start[elems[e].level]++;

  uint *level_off = *level_off_ = tcalloc(uint, nlevel + 1);
  slong current = 0;
  for (e = nlevel - 1; e >= 0; e--) {
    current += start[e];
    level_off[e + 1] = current - start[e];
  }
  level_off[0] = current;

  slong *out = tcalloc(slong, 2 * nlevel);
  slong *buf = tcalloc(slong, 2 * nlevel);
  comm_scan(out, gc, gs_long, gs_add, start, nlevel, buf);

  current = 0;
  for (e = 2 * nlevel - 1; e >= nlevel; e--) {
    current += out[e];
    out[e] = current - out[e];
  }

  for (e = 0; e < nelt; e++)
    elems[e].globalId = ++out[elems[e].level] + out[nlevel + elems[e].level];

  /* FIXME: comms[e] needs to be set in RSB numbering phase */
  MPI_Comm *comms = *comms_ = tcalloc(MPI_Comm, nlevel);
  for (e = 0; e < nlevel; e++) {
    int bin = start[e] > 0 ? 1 : 0;
    MPI_Comm_split(gc->c, bin, gc->id, &comms[e]);
  }

  free(out);
  free(buf);
  free(start);

  return nlevel;
}

static slong id_interface_elements(int level, uint nelt, int nv,
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
    for (e = 0; e < nelt; e++)
      for (v = 0; v < nv; v++)
        interface[e * nv + v] = 1;
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

static int number_elements(MPI_Comm **comms, uint **level_off,
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
  int level = 0;

  while (lc->np > 1) {
    /* Run RCB, RIB pre-step or just sort by global id */
    if (h->options->rsb_prepartition == 1) // RCB
      rcb(lc, h->elements, ndim, &h->buf);
    else if (h->options->rsb_prepartition == 2) // RIB
      rib(lc, h->elements, ndim, &h->buf);

    /* Run fiedler */
    int ipass = 0, iter;
    do {
      genmap_vector ivec;
      genmap_init_vector(&ivec, ipass == 0, lc, h);

      GenmapFiedler(h, lc, max_iter, ivec);

      genmap_destroy_vector(ivec);
    } while (++ipass < max_pass && iter == max_iter);

    /* Sort by Fiedler vector */
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, lc,
                  &h->buf);

    /* Bisect */
    int bin = 1;
    if (lc->id < (lc->np + 1) / 2)
      bin = 0;
    repair_partitions(h, bin, level, lc, gc);

    nelt = genmap_get_nel(h);
    e = genmap_get_elements(h);
    id_interface_elements(level, nelt, nv, e, bin, lc, &h->buf, gc);

    struct comm tc;
    genmap_comm_split(lc, bin, lc->id, &tc);
    comm_free(lc);
    comm_dup(lc, &tc);
    comm_free(&tc);

    genmap_comm_scan(h, lc);

    level++;
  }

  e = genmap_get_elements(h);
  nelt = genmap_get_nel(h);
  uint i;
  int unmarked = 0;
  for (i = 0; i < nelt; i++) {
    e[i].globalId = 0;
    if (e[i].level < 0) {
      e[i].level = level;
      unmarked = 1;
    }
  }
  level += unmarked;

  return number_elements_aux(comms, level_off, level, nelt, e, gc);
}

struct csr_mat_ *parrsb_numbering(unsigned int *nelt_, unsigned int *nlevels_,
                                  unsigned int **level_off, MPI_Comm **comms,
                                  long long *vl, double *coord, int nv,
                                  MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  /* Load balance input data */
  uint nelt = *nelt_;
  struct array elems;
  genmap_load_balance(&elems, nelt, nv, coord, vl, &cr, &bfr);

  /* FIXME: If MG, we need both CSR and gs implementation now */
  parRSB_options options = parrsb_default_options;
  options.rsb_laplacian_implementation = 2;

  genmap_handle h;
  genmap_init(&h, comm, &options);

  genmap_set_elements(h, &elems);
  genmap_set_nvertices(h, nv);
  genmap_comm_scan(h, genmap_global_comm(h));

  int nlevels = *nlevels_ = number_elements(comms, level_off, h);
  nelt = *nelt_ = genmap_get_nel(h);

  ulong *gid = tcalloc(ulong, nelt);
  struct rsb_element *e = genmap_get_elements(h);
  uint i;
  for (i = 0; i < nelt; i++)
    gid[i] = e[i].globalId;

  struct csr_mat_ *M = NULL;
  genmap_laplacian_csr_init(&M, gid, h, &c);

  free(gid);
  genmap_finalize(h);
  array_free(&elems);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return M;
}

genmap_handle parrsb_numbering_w_handle(unsigned int *nelt, long long *vl,
                                        double *coord, int nv, MPI_Comm comm) {
  parRSB_options options = parrsb_default_options;
  options.rsb_laplacian_implementation = 2;

  genmap_handle h;
  genmap_init(&h, comm, &options);

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer bfr;
  buffer_init(&bfr, 1024);

  struct array elems;
  genmap_load_balance(&elems, *nelt, nv, coord, vl, &cr, &bfr);

  genmap_set_elements(h, &elems);
  genmap_set_nvertices(h, nv);

  struct comm *gc = genmap_global_comm(h);
  genmap_comm_scan(h, gc);

  h->nlevels = number_elements(&h->comms, &h->level_off, h);

  *nelt = genmap_get_nel(h);
  ulong *gid = tcalloc(ulong, *nelt);
  struct rsb_element *e = genmap_get_elements(h);
  uint i;
  for (i = 0; i < *nelt; i++)
    gid[i] = e[i].globalId;

  genmap_laplacian_csr_init(&h->M, gid, h, &c);

  free(gid);
  array_free(&elems);
  buffer_free(&bfr);
  crystal_free(&cr);
  comm_free(&c);

  return h;
}

/*
 * ILU0
 */
struct row_info {
  ulong gid;
  uint dest;
  uint owner;
  uint proc;
};

struct csr_entry {
  ulong r, c;
  uint proc;
  GenmapScalar v;
};

/* TODO: Do a binary search */
static int find_gid(ulong id, ulong *gid, uint n) {
  uint i;
  for (i = 0; i < n; i++)
    if (gid[i] == id)
      return 1;

  return 0;
}

static int find_dependent_procs(uint **offsets_, uint **procs_,
                                struct csr_mat_ *M, struct comm *c, int fw,
                                buffer *buf) {
  uint rn = M->rn;
  ulong *col = M->col;

  struct array rows;
  array_init(struct row_info, &rows, rn);

  struct row_info ri;
  ri.dest = c->id;

  sint i, j;
  if (fw == 1) {
    for (i = 0; i < rn; i++) {
      ulong id = M->row_id[i];
      for (j = M->row_off[i]; j < M->row_off[i + 1] && col[j] <= id; j++) {
        ri.gid = col[j];
        ri.proc = col[j] % c->np;
        ri.owner = find_gid(col[j], M->row_id, rn);
        array_cat(struct row_info, &rows, &ri, 1);
      }
    }
  } else {
    for (i = 0; i < rn; i++) {
      ulong id = M->row_id[i];
      for (j = M->row_off[i + 1] - 1; j >= (sint)M->row_off[i] && col[j] >= id;
           j--) {
        ri.gid = col[j];
        ri.proc = col[j] % c->np;
        ri.owner = find_gid(col[j], M->row_id, rn);
        array_cat(struct row_info, &rows, &ri, 1);
      }
    }
  }

  /* Remove local duplicates in rows */
  struct array rows_compact;
  array_init(struct row_info, &rows_compact, rn);

  uint n = rows.n;

  if (n > 0) {
    struct row_info *ptr = rows.ptr;
    sarray_sort(struct row_info, ptr, n, gid, 1, buf);

    uint cur_idx = 0;
    int cur_owner = ptr[0].owner > 0 ? 1 : 0;

    for (i = 1; i < n; i++) {
      if (ptr[i].gid != ptr[cur_idx].gid) {
        ri = ptr[cur_idx];
        ri.owner = cur_owner;
        array_cat(struct row_info, &rows_compact, &ri, 1);
        cur_idx = i;
        cur_owner = ptr[i].owner;
      } else if (ptr[i].owner > 0) {
        cur_owner = 1;
      }
    }
    ri = ptr[cur_idx];
    ri.owner = cur_owner;
    array_cat(struct row_info, &rows_compact, &ri, 1);
  }

  /* Send to  work-proc and identify the owner */
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct row_info, &rows_compact, proc, 0, &cr);

  n = rows_compact.n;
  if (n > 0) {
    struct row_info *ptr = rows_compact.ptr;
    sarray_sort(struct row_info, ptr, n, gid, 1, buf);

    uint cur_idx = 0;
    int cur_proc = ptr[0].owner > 0 ? ptr[0].dest : -1;

    for (i = 1; i < n; i++) {
      if (ptr[i].gid != ptr[cur_idx].gid) {
        for (j = cur_idx; j < i; j++)
          ptr[j].proc = cur_proc;
        cur_idx = i;
        cur_proc = ptr[i].owner > 0 ? ptr[i].dest : -1;
      } else if (ptr[i].owner > 0) {
        cur_proc = ptr[i].dest;
      }
    }
    for (j = cur_idx; j < n; j++)
      ptr[j].proc = cur_proc;
  }

  /* Send to owner */
  sarray_transfer(struct row_info, &rows_compact, proc, 0, &cr);

  struct row_info *ptr = rows_compact.ptr;
  n = rows_compact.n;
  sarray_sort_2(struct row_info, ptr, n, gid, 1, dest, 0, buf);

  uint *offsets = *offsets_ = tcalloc(uint, M->rn + 1);
  uint *procs = *procs_ = tcalloc(uint, n);

  uint ngid = 0;
  offsets[0] = 0;
  for (i = 0; i < n; i++) {
    if (ptr[i].gid != ptr[offsets[ngid]].gid) {
      ngid++;
      offsets[ngid] = i;
    }
    procs[i] = ptr[i].dest;
  }
  ngid++;
  offsets[ngid] = n;

  array_free(&rows);
  array_free(&rows_compact);

  return 0;
}

static int append_dependent_rows(struct array *entries, int level,
                                 uint *level_off, uint *offsets, uint *procs,
                                 struct csr_mat_ *M, buffer *buf) {
  uint s = level_off[level + 1];
  uint e = level_off[level];

  uint rn = M->rn;
  ulong *col = M->col;

  struct csr_entry t;
  uint i, j, p;
  for (i = s; i < e; i++) {
    t.r = M->row_id[i];
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++) {
      t.c = M->col[j];
      t.v = M->v[j];
      for (p = offsets[i]; p < offsets[i + 1]; p++) {
        t.proc = procs[p];
        array_cat(struct csr_entry, entries, &t, 1);
      }
    }
  }

  return 0;
}

static int ilu0_level(struct csr_mat_ *M, struct csr_mat_ *N, int lvl,
                      int nlevels, unsigned int *lvl_off) {
  const uint rn = M->rn;
  const uint *off = M->row_off;
  double *v = M->v;
  const ulong *col = M->col;

  double a_kk, a_kj, a_ik;
  uint i, ii, j, k, kk, ik, ij;

  uint start = lvl_off[lvl + 1];
  uint end = lvl_off[lvl];

  ii = start;
  if (lvl == nlevels - 1)
    ii = start + 1;

  for (; ii < end; ii++) { /* Go over number of rows */
    i = M->row_id[ii];

    for (kk = M->row_off[ii]; kk < M->row_off[ii + 1] && M->col[kk] < i; kk++) {
      k = M->col[kk];
      if (csr_mat_get_global(&a_kk, NULL, M, k, k) != 0)
        csr_mat_get_global(&a_kk, NULL, N, k, k);
      if (fabs(a_kk) < 1e-12)
        continue;

      csr_mat_get_local(&a_ik, &ik, M, ii + 1, k);
      if (fabs(a_ik) < 1e-12)
        continue;

      /* a_ik = a_ik / a_kk */
      a_ik = v[ik] = v[ik] / a_kk;

      /* a_ij = a_ij - a_ik * a_kj */
      for (ij = ik + 1; ij < off[ii + 1]; ij++) { /* Go over the columns */
        j = col[ij];
        if (csr_mat_get_global(&a_kj, NULL, M, k, j) != 0)
          csr_mat_get_global(&a_kj, NULL, N, k, j);
        v[ij] = v[ij] - a_ik * a_kj;
      }
    }
  }

  return 0;
}

static int ilu0_aux(uint nlevels, uint *level_off, MPI_Comm *comms,
                    struct csr_mat_ *M, struct comm *gc, buffer *buf) {
  struct crystal cr;
  crystal_init(&cr, gc);

  uint *offsets, *procs;
  find_dependent_procs(&offsets, &procs, M, gc, 1, buf);

  struct array rows_ext;
  array_init(struct csr_entry, &rows_ext, 10);

  int i;
  for (i = nlevels - 1; i >= 0; i--) {
    int active = (level_off[i + 1] != level_off[i]) ? 1 : 0;

    sint np = 0;
    struct comm c;
    if (active) {
      comm_init(&c, comms[i]);
      np = c.np;
    }

    sint bfr;
    comm_allreduce(gc, gs_int, gs_max, &np, 1, &bfr);

    int p;
    for (p = 0; p < np; p++) {
      /* Do the local ilu */
      if (c.id == p && active) {
        struct csr_mat_ *N = NULL;
        csr_mat_setup(&N, &rows_ext, NULL, buf);
        ilu0_level(M, N, i, nlevels, level_off);
        csr_mat_free(N);
        append_dependent_rows(&rows_ext, i, level_off, offsets, procs, M, buf);
      }

      sarray_transfer(struct csr_entry, &rows_ext, proc, 0, &cr);
    }

    if (active)
      comm_free(&c);
  }

  if (offsets != NULL)
    free(offsets);
  if (procs != NULL)
    free(procs);

  array_free(&rows_ext);
  crystal_free(&cr);

  return 0;
}

void parrsb_ilu0(unsigned int nlevels, unsigned int *level_off, MPI_Comm *comms,
                 struct csr_mat_ *M, MPI_Comm world) {
  buffer buf;
  buffer_init(&buf, 1024);

  struct comm gc;
  comm_init(&gc, world);

  ilu0_aux(nlevels, level_off, comms, M, &gc, &buf);

  comm_free(&gc);
  buffer_free(&buf);
}

/*
 * ILUT
 */
static int update_w(struct array *w, double wk, ulong k, struct csr_mat_ *M,
                    struct csr_mat_ *N, buffer *buf) {
  /* Find k in M and then in N, do a binary search */
  int found = 0;
  struct csr_mat_ *U;

  uint i;
  for (i = 0; i < M->rn; i++) {
    if (k == M->row_id[i]) {
      found = 1;
      U = M;
      break;
    }
  }

  if (found == 0) {
    for (i = 0; i < N->rn; i++)
      if (k == N->row_id[i]) {
        found = 1;
        U = N;
        break;
      }
  }

  /* Update w */
  if (found == 1) {
    sarray_sort(struct csr_entry, w->ptr, w->n, c, 1, buf);

    uint s = U->row_off[i];
    uint e = U->row_off[i + 1];

    struct csr_entry t;
    t.r = k;

    uint ai = 0;
    struct csr_entry *ptr = w->ptr;

    for (i = s; i < e; i++) {
      t.c = U->col[i];
      t.v = -wk * U->v[i];
      while (ai < w->n && ptr[ai].c < t.c)
        ai++;
      if (ai < w->n && ptr[ai].c == t.c)
        ptr[ai].v += t.v;
      else {
        array_cat(struct csr_entry, w, &t, 1);
        sarray_sort(struct csr_entry, w->ptr, w->n, c, 1, buf);
        ai = 0;
        ptr = w->ptr;
      }
    }
  }

  return found == 0 ? 1: 0;
}

static int ilut_level(struct csr_mat_ *M, struct csr_mat_ *U, int lvl,
                      int nlevels, unsigned int *lvl_off, buffer *buf) {
  const uint rn = M->rn;
  const uint *off = M->row_off;
  double *v = M->v;
  const ulong *col = M->col;

  double a_kk, a_kj, a_ik;
  uint i, ii, j, k, kk, ik, ij;

  uint start = lvl_off[lvl + 1];
  uint end = lvl_off[lvl];

  double tau = 0.1, taui;
  double norm;
  double wk;

  struct csr_entry t;

  struct array w;
  array_init(struct csr_entry, &w, 10);

  /* Go over number of rows */
  for (ii = start; ii < end; ii++) {
    /* Initialize w and calculate norm of row i */
    w.n = 0;
    i = t.r = M->row_id[ii];

    norm = 0.0;
    for (kk = M->row_off[ii]; kk < M->row_off[ii + 1]; kk++) {
      t.c = M->col[kk];
      t.v = M->v[kk];
      norm += t.v * t.v;
      array_cat(struct csr_entry, &w, &t, 1);
    }
    norm = sqrt(norm);
    taui = norm * tau;

    for (kk = M->row_off[ii]; kk < M->row_off[ii + 1] && M->col[kk] < i; kk++) {
      k = M->col[kk];
      if (csr_mat_get_global(&a_kk, NULL, M, k, k) != 0)
        csr_mat_get_global(&a_kk, NULL, U, k, k);
      if (fabs(a_kk) < 1e-10)
        continue;

      wk = fabs(M->v[kk]) / a_kk;
      if (wk >= taui)
        update_w(&w, wk, k, M, U, buf); /* w = w - w_k * u_{k,*} */
    }
    /* Apply dropping rule to w */
  }

  return 0;
}

static int ilut_aux(uint nlevels, uint *level_off, MPI_Comm *comms,
                    struct csr_mat_ *M, struct comm *gc, buffer *buf) {
  struct crystal cr;
  crystal_init(&cr, gc);

  uint *offsets, *procs;
  find_dependent_procs(&offsets, &procs, M, gc, 1, buf);

  struct array rows_ext;
  array_init(struct csr_entry, &rows_ext, 10);

  int i;
  for (i = nlevels - 1; i >= 0; i--) {
    int active = (level_off[i + 1] != level_off[i]) ? 1 : 0;

    sint np = 0;
    struct comm c;
    if (active) {
      comm_init(&c, comms[i]);
      np = c.np;
    }

    sint bfr;
    comm_allreduce(gc, gs_int, gs_max, &np, 1, &bfr);

    int p;
    for (p = 0; p < np; p++) {
      /* Do the local ilu */
      if (c.id == p && active) {
        struct csr_mat_ *N = NULL;
        csr_mat_setup(&N, &rows_ext, NULL, buf);
        ilut_level(M, N, i, nlevels, level_off, buf);
        csr_mat_free(N);
        append_dependent_rows(&rows_ext, i, level_off, offsets, procs, M, buf);
      }

      sarray_transfer(struct csr_entry, &rows_ext, proc, 0, &cr);
    }

    if (active)
      comm_free(&c);
  }

  if (offsets != NULL)
    free(offsets);
  if (procs != NULL)
    free(procs);

  array_free(&rows_ext);
  crystal_free(&cr);

  return 0;
}

void parrsb_ilut(unsigned int nlevels, unsigned int *level_off, MPI_Comm *comms,
                 struct csr_mat_ *M, MPI_Comm world) {
  buffer buf;
  buffer_init(&buf, 1024);

  struct comm gc;
  comm_init(&gc, world);

  ilut_aux(nlevels, level_off, comms, M, &gc, &buf);

  comm_free(&gc);
  buffer_free(&buf);
}

/*
 * LU-solve
 */
struct rhs_entry {
  ulong r;
  uint proc;
  GenmapScalar v;
};

/* TODO: Use binary search */
static int get_rhs_ext(double *x, ulong gid, struct array *rhs) {
  int found = 0;

  struct rhs_entry *ptr = rhs->ptr;
  uint i;
  for (i = 0; i < rhs->n; i++)
    if (ptr[i].r == gid) {
      *x = ptr[i].v;
      found = 1;
    }

  return found == 1 ? 0 : 1;
}

static int get_rhs(double *x, ulong gid, ulong *ids, int level,
                   unsigned int *level_off, double *y) {
  int found = 0;

  uint s = level_off[level + 1];
  uint e = level_off[level];
  uint i;
  for (i = s; i < e; i++)
    if (gid == ids[i]) {
      *x = y[i];
      found = 1;
      break;
    }

  return found == 1 ? 0 : 1;
}

static int append_dependent_rhs(struct array *rhs, int level, uint *level_off,
                                uint *offsets, uint *procs, double *y,
                                ulong *row_id, buffer *buf) {
  uint s = level_off[level + 1];
  uint e = level_off[level];

  struct rhs_entry t;
  uint i, p;
  for (i = s; i < e; i++) {
    t.r = row_id[i];
    t.v = y[i];
    for (p = offsets[i]; p < offsets[i + 1]; p++) {
      t.proc = procs[p];
      array_cat(struct rhs_entry, rhs, &t, 1);
    }
  }

  return 0;
}

static int fw_solve_level(double *y, struct csr_mat_ *A, double *b,
                          struct array *rhs, int level, unsigned int *level_off,
                          struct comm *gc) {
  /* Forward substitution with L - Lower matrix (diagonal = all 1s)
   * L(Ux) = L(y) = b */
  double r;
  sint i, j;
  for (i = level_off[level + 1]; i < level_off[level]; i++) {
    y[i] = b[i];
    ulong id = A->row_id[i];
    /* Remove the dot product: L[i, j] * y[j], j < i */
    for (j = A->row_off[i]; j < A->row_off[i + 1] && A->col[j] < id; j++) {
      if (get_rhs(&r, A->col[j], A->row_id, level, level_off, y) != 0)
        get_rhs_ext(&r, A->col[j], rhs);
      y[i] -= A->v[j] * r;
    }
    // printf("0: id = %lld b[%d] = %lf y[%d] = %lf\n", id, i, b[i], i, y[i]);
    // fflush(stdout);
  }

  return 0;
}

static int fw_solve_aux(double *y, struct csr_mat_ *A, double *b, int nlevels,
                        unsigned int *level_off, MPI_Comm *comms,
                        struct comm *gc, struct crystal *cr, buffer *buf) {
  uint *offsets, *procs;
  find_dependent_procs(&offsets, &procs, A, gc, 1, buf);

  struct array rhs_ext;
  array_init(struct rhs_entry, &rhs_ext, 10);

  int i;
  for (i = nlevels - 1; i >= 0; i--) {
    int active = (level_off[i + 1] != level_off[i]) ? 1 : 0;

    sint np = 0;
    struct comm c;
    if (active) {
      comm_init(&c, comms[i]);
      np = c.np;
    }

    sint bfr;
    comm_allreduce(gc, gs_int, gs_max, &np, 1, &bfr);

    /* Do LU solve in serial among the processors in the current level */
    int p;
    for (p = 0; p < np; p++) {
      if (c.id == p && active) {
        fw_solve_level(y, A, b, &rhs_ext, i, level_off, gc);
        append_dependent_rhs(&rhs_ext, i, level_off, offsets, procs, y,
                             A->row_id, buf);
      }

      sarray_transfer(struct rhs_entry, &rhs_ext, proc, 0, cr);
    }

    if (active)
      comm_free(&c);
  }

  if (offsets != NULL)
    free(offsets);
  if (procs != NULL)
    free(procs);

  array_free(&rhs_ext);

  return 0;
}

static int bw_solve_level(double *x, struct csr_mat_ *A, double *y,
                          struct array *rhs, int level, unsigned int *level_off,
                          struct comm *gc) {
  /* Back substitution with U, Ux = y */
  double r;
  sint i, j;
  for (i = (sint)level_off[level] - 1; i >= (sint)level_off[level + 1]; i--) {
    x[i] = y[i];
    ulong id = A->row_id[i];
    /* Remove the dot product: L[i, j] * x[j], j > i */
    for (j = A->row_off[i + 1] - 1; j >= (sint)A->row_off[i] && A->col[j] > id;
         j--) {
      if (get_rhs(&r, A->col[j], A->row_id, level, level_off, x) != 0)
        get_rhs_ext(&r, A->col[j], rhs);
      x[i] -= A->v[j] * r;
    }
    if (fabs(A->v[j]) > 1e-10)
      x[i] /= A->v[j];
    // printf("1: id = %lld y[%d] = %lf x[%d] = %lf\n", id, i, y[i], i, x[i]);
    // fflush(stdout);
  }

  return 0;
}

static int bw_solve_aux(double *y, struct csr_mat_ *A, double *b, int nlevels,
                        unsigned int *level_off, MPI_Comm *comms,
                        struct comm *gc, struct crystal *cr, buffer *buf) {
  uint *offsets, *procs;
  find_dependent_procs(&offsets, &procs, A, gc, 0, buf);

  struct array rhs_ext;
  array_init(struct rhs_entry, &rhs_ext, 10);

  int i;
  for (i = 0; i < nlevels; i++) {
    int active = (level_off[i + 1] != level_off[i]) ? 1 : 0;

    sint np = 0;
    struct comm c;
    if (active) {
      comm_init(&c, comms[i]);
      np = c.np;
    }

    sint bfr;
    comm_allreduce(gc, gs_int, gs_max, &np, 1, &bfr);

    /* Do LU solve in serial among the processors in the current level */
    int p;
    for (p = np - 1; p >= 0; p--) {
      if (c.id == p && active) {
        bw_solve_level(y, A, b, &rhs_ext, i, level_off, gc);
        append_dependent_rhs(&rhs_ext, i, level_off, offsets, procs, y,
                             A->row_id, buf);
      }

      sarray_transfer(struct rhs_entry, &rhs_ext, proc, 0, cr);
    }

    if (active)
      comm_free(&c);
  }

  if (offsets != NULL)
    free(offsets);
  if (procs != NULL)
    free(procs);

  array_free(&rhs_ext);

  return 0;
}

void parrsb_lu_solve(double *x, struct csr_mat_ *M, double *b,
                     unsigned int nlevels, unsigned int *level_off,
                     MPI_Comm *comms, MPI_Comm world) {
  buffer buf;
  buffer_init(&buf, 1024);

  struct comm gc;
  comm_init(&gc, world);

  struct crystal cr;
  crystal_init(&cr, &gc);

  double *y = tcalloc(double, M->rn);
  fw_solve_aux(y, M, b, nlevels, level_off, comms, &gc, &cr, &buf);
  bw_solve_aux(x, M, y, nlevels, level_off, comms, &gc, &cr, &buf);
  free(y);

  crystal_free(&cr);
  comm_free(&gc);
  buffer_free(&buf);
}
