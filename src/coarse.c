#include "coarse-impl.h"
#include "metrics.h"
#include "sort.h"
#include <float.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static uint unique_ids(sint *perm, ulong *uid, uint n, const ulong *ids,
                       buffer *bfr);

//------------------------------------------------------------------------------
// Setup coarse grid system. Initial dumb API.
//
// Number rows, local first then interface. Returns global number of local
// elements.
struct rcb_t {
  uint i, s;
  double coord[3];
  slong vtx[8];
};

static void nmbr_local_rcb(struct array *a, uint s, uint e, const unsigned nc,
                           const unsigned ndim, const unsigned level,
                           struct comm *c, buffer *bfr) {
  sint size = e - s;
  if (size <= 1)
    return;

  double max[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX},
         min[3] = {DBL_MAX, DBL_MAX, DBL_MAX};

  struct rcb_t *pa = (struct rcb_t *)a->ptr;
  for (uint i = s; i < e; i++) {
    for (int j = 0; j < ndim; j++) {
      if (pa[i].coord[j] < min[j])
        min[j] = pa[i].coord[j];
      if (pa[i].coord[j] > max[j])
        max[j] = pa[i].coord[j];
    }
  }

  double len = max[0] - min[0];
  int axis = 0;
  for (int j = 1; j < ndim; j++) {
    if (max[j] - min[j] > len)
      axis = j, len = max[j] - min[j];
  }

  struct rcb_t *ps = pa + s;
  switch (axis) {
  case 0:
    sarray_sort(struct rcb_t, ps, size, coord[0], 3, bfr);
    break;
  case 1:
    sarray_sort(struct rcb_t, ps, size, coord[1], 3, bfr);
    break;
  case 2:
    sarray_sort(struct rcb_t, ps, size, coord[2], 3, bfr);
    break;
  default:
    break;
  }

  // Number the elements in the interface
  uint npts = size * nc;
  slong *vtx = tcalloc(slong, npts);
  for (uint i = s, k = 0; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      vtx[k] = pa[i].vtx[j];
  }

  struct gs_data *gsh = gs_setup(vtx, npts, c, 0, gs_pairwise, 0);

  sint *dof = tcalloc(sint, npts);
  uint mid = (s + e) / 2;
  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 1;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid, k = (mid - s) * nc; i < e; i++) {
    for (int j = 0; j < nc; j++, k++)
      dof[k] = 0;
  }

  gs(dof, gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s, k = 0; i < e; i++, k++) {
    for (int j = 0; j < nc; j++) {
      if (dof[k * nc + j] > 0 && pa[i].s == INT_MAX) {
        pa[i].s = level;
        break;
      }
    }
  }

  gs_free(gsh);
  free(dof), free(vtx);

  nmbr_local_rcb(a, s, mid, nc, ndim, level + 1, c, bfr);
  nmbr_local_rcb(a, mid, e, nc, ndim, level + 1, c, bfr);
}

// Number the DOFs internal first, faces second and all the rest (wire basket)
// next. This keeps zeros as is and renumber the positive entries in `ids`
// array.
static void number_dual_graph_dofs(ulong *dofs, struct coarse *crs, uint n,
                                   const slong *ids, uint nelt, unsigned ndim,
                                   const scalar *coord, buffer *bfr) {
  int nnz = (n > 0);
  struct comm c;
  comm_split(&crs->c, nnz, crs->c.id, &c);

  unsigned nc = n / nelt;
  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, n);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < n; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < n; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < nelt; i++) {
        for (j = 0; j < nc; j++) {
          if (dof[i * nc + j] > 0 && !dofs[i]) {
            dofs[i] = level;
            break;
          }
        }
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      level++;
    }
    free(dof);
  }

  for (i = crs->n[0] = crs->n[1] = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      crs->n[1]++;
    else
      crs->n[0]++;
  }

  slong in[2] = {crs->n[0], crs->n[1]}, out[2][2], wrk[2][2];
  comm_scan(out, &crs->c, gs_long, gs_add, in, 2, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];
  crs->s[1] = out[0][1] + 1, crs->ng[1] = out[1][1];

  struct array local;
  array_init(struct rcb_t, &local, crs->n[0]);

  struct rcb_t t = {.s = INT_MAX};
  ulong s = crs->ng[0] + crs->s[1];
  for (uint i = 0; i < nelt; i++) {
    if (dofs[i] > 0)
      dofs[i] = s++;
    else {
      t.i = i;
      memcpy(t.coord, &coord[i * ndim], ndim * sizeof(scalar));
      memcpy(t.vtx, &ids[i * nc], nc * sizeof(slong));
      array_cat(struct rcb_t, &local, &t, 1);
    }
  }

  if (local.n > 0) {
    nmbr_local_rcb(&local, 0, local.n, nc, ndim, 1, &c, bfr);
    sarray_sort(struct rcb_t, local.ptr, local.n, s, 0, bfr);
    struct rcb_t *pl = (struct rcb_t *)local.ptr;
    ulong s = crs->s[0];
    for (sint i = local.n - 1; i >= 0; i--)
      dofs[pl[i].i] = s++;
  }

  comm_free(&c);
  array_free(&local);
}

struct coarse *coarse_setup(unsigned n, unsigned nc, const long long *vl,
                            const scalar *coord, unsigned null_space,
                            unsigned type, struct comm *c) {
  struct coarse *crs = tcalloc(struct coarse, 1);
  crs->type = type, crs->null_space = null_space, crs->un = n;

  // Setup the buffer and duplicate the communicator.
  buffer_init(&crs->bfr, 1024);
  comm_dup(&crs->c, c);

  uint size = n * nc;
  slong *tid = tcalloc(slong, size);
  for (uint i = 0; i < size; i++)
    tid[i] = vl[i];

  ulong *nid = tcalloc(ulong, n);
  unsigned ndim = (nc == 8) ? 3 : 2;
  number_dual_graph_dofs(nid, crs, size, tid, n, ndim, coord, &crs->bfr);

  // Find unique ids and user vector to compressed vector mapping.
  // In the case of dual-graph Laplacian, all the ids are unique.
  // But here we arrange them in the sorted order.
  ulong *uid = tcalloc(ulong, n);
  crs->u2c = tcalloc(sint, n);
  crs->cn = unique_ids(crs->u2c, uid, n, nid, &crs->bfr);
  assert(crs->cn == crs->un);
  crs->an = crs->cn;

  struct crystal cr;
  crystal_init(&cr, &crs->c);

  struct array nbrs, eij;
  find_nbrs(&nbrs, nid, tid, n, nc, &cr, &crs->bfr);
  // Convert `struct nbr` -> `struct mij` and compress entries which share the
  // same (r, c) values. Set the diagonal element to have zero row sum
  compress_nbrs(&eij, &nbrs, &crs->bfr);
  array_free(&nbrs);

  switch (type) {
  case 0:
    schur_setup(crs, &eij, &cr, &crs->bfr);
    break;
  default:
    break;
  }

  array_free(&eij), crystal_free(&cr);
  free(tid), free(nid), free(uid);

  return crs;
}

void coarse_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol) {
  metric_init();

  scalar *rhs = tcalloc(scalar, 2 * crs->an), *xx = rhs + crs->an;
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      rhs[crs->u2c[i]] += b[i];
  }

  switch (crs->type) {
  case 0:
    schur_solve(xx, crs, rhs, tol, &crs->bfr);
    break;
  default:
    break;
  }

  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      x[i] = xx[crs->u2c[i]];
  }
  free(rhs);

  metric_push_level();
  metric_crs_print(&crs->c, 1);
}

void coarse_free(struct coarse *crs) {
  if (crs != NULL) {
    switch (crs->type) {
    case 0:
      schur_free(crs);
      break;
    default:
      break;
    }
    if (crs->u2c)
      free(crs->u2c);
    comm_free(&crs->c), buffer_free(&crs->bfr);
    free(crs), crs = NULL;
  }
}

//------------------------------------------------------------------------------
// Better API for coarse grid system.
//
static uint unique_ids(sint *perm, ulong *uid, uint n, const ulong *ids,
                       buffer *bfr) {
  struct id_t {
    ulong id;
    uint idx;
    sint perm;
  };

  struct array arr;
  array_init(struct id_t, &arr, n);

  uint i;
  struct id_t t = {.id = 0, .idx = 0, .perm = -1};
  for (i = 0; i < n; i++) {
    t.id = ids[i], t.idx = i;
    array_cat(struct id_t, &arr, &t, 1);
  }

  sarray_sort(struct id_t, arr.ptr, arr.n, id, 1, bfr);
  struct id_t *pa = (struct id_t *)arr.ptr;

  // Ignore the ids numbered zero
  for (i = 0; i < arr.n && pa[i].id == 0; i++)
    ;

  uint un = 0;
  ulong last = 0;
  for (; i < arr.n; i++) {
    ulong v = pa[i].id;
    if (v != last)
      last = uid[un] = v, un++;
    pa[i].perm = un - 1;
  }

  sarray_sort(struct id_t, pa, n, idx, 0, bfr);
  pa = (struct id_t *)arr.ptr;
  for (i = 0; i < n; i++)
    perm[i] = pa[i].perm;

  array_free(&arr);
  return un;
}

// Number rows, local first then interface. Returns global number of local
// elements.
struct rsb_t {
  uint i, s;
  slong vtx[8];
};

static void number_dofs(slong *nid, struct coarse *crs, const slong *ids,
                        const ulong *uid) {
  uint un = crs->un;
  buffer *bfr = &crs->bfr;
  struct comm *ci = &crs->c;
  sint *u2c = crs->u2c;

  int nnz = (un > 0);
  struct comm c;
  comm_split(ci, nnz, ci->id, &c);

  uint i, j;
  if (nnz) {
    sint *dof = tcalloc(sint, un);
    int level = 1;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, un, &c, 0, gs_pairwise, 0);

      int bin = (c.id >= (c.np + 1) / 2);
      for (i = 0; i < un; i++)
        dof[i] = bin;

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      if (bin == 1) {
        for (i = 0; i < un; i++)
          dof[i] = 0;
      }

      gs(dof, gs_int, gs_add, 0, gsh, bfr);

      for (i = 0; i < un; i++) {
        if (dof[i] > 0 && u2c[i] >= 0 && !nid[u2c[i]])
          nid[u2c[i]] = level;
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      level++;
    }
    free(dof);
  }

  // Calculate unqiue local and interface nodes based on compress ids.
  // Finding unique local ids is easy. To find unique interface ids, we
  // will have to sort in parallel and then manually find the unique ids.
  struct dof_t {
    ulong id, nid;
    uint p, p0, idx;
  };

  struct array arr;
  array_init(struct dof_t, &arr, crs->cn);

  uint ln = 0;
  struct dof_t t = {.id = 0, .nid = 0, .p = 0, .p0 = ci->id, .idx = 0};
  for (i = 0; i < crs->cn; i++) {
    if (!nid[i])
      ln++;
    else
      t.id = uid[i], t.idx = i, array_cat(struct dof_t, &arr, &t, 1);
  }
  crs->n[0] = ln;

  slong cnt[1] = {ln}, out[2][1], wrk[2][1];
  comm_scan(out, ci, gs_long, gs_add, cnt, 1, wrk);
  crs->s[0] = out[0][0] + 1, crs->ng[0] = out[1][0];
  // printf("id = %u n[0] = %u s[0] = %llu ng[0] = %llu\n", ci->id, crs->n[0],
  //        crs->s[0], crs->ng[0]);

  for (i = 0, ln = 0; i < crs->cn; i++) {
    if (!nid[i])
      nid[i] = crs->s[0] + ln, ln++;
  }
  assert(crs->n[0] == ln);

  // parallel_sort and set nid and send back to p0
  parallel_sort(struct dof_t, &arr, id, gs_long, 0, 0, ci, bfr);

  uint in = 0;
  if (arr.n > 0) {
    struct dof_t *pa = (struct dof_t *)arr.ptr;
    for (i = in = 1; i < arr.n; i++)
      in += (pa[i].id != pa[i - 1].id);
  }

  cnt[0] = in;
  comm_scan(out, ci, gs_long, gs_add, cnt, 1, wrk);
  crs->ng[1] = out[1][0];
  slong s = crs->ng[0] + out[0][0] + 1;
  // printf("id = %u n[1] = %u s[1] = %llu ng[1] = %llu\n", ci->id, in, s,
  //        crs->ng[1]);

  if (in) {
    struct dof_t *pa = (struct dof_t *)arr.ptr;
    i = 0;
    while (i < arr.n) {
      for (j = i + 1; j < arr.n && pa[j].id == pa[i].id; j++)
        ;
      for (; i < j; i++)
        pa[i].nid = s;
      s++;
    }
  }

  struct crystal cr;
  crystal_init(&cr, ci);
  sarray_transfer(struct dof_t, &arr, p0, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct dof_t, arr.ptr, arr.n, id, 1, bfr);
  struct dof_t *pa = (struct dof_t *)arr.ptr;
  for (i = 0; i < arr.n; i++)
    nid[pa[i].idx] = pa[i].nid;

  array_free(&arr);
  comm_free(&c);
}

// n  = ncr * nelt
// nz = ncr * ncr * nelt
struct coarse *crs_parrsb_setup(uint n, const ulong *id, uint nz, const uint *Ai,
                               const uint *Aj, const scalar *A,
                               unsigned null_space, unsigned type,
                               const struct comm *c) {
  struct coarse *crs = tcalloc(struct coarse, 1);
  // crs->un is the user vector size.
  crs->null_space = null_space, crs->type = type, crs->un = n;

  // Setup the buffer and duplicate the communicator.
  buffer_init(&crs->bfr, 1024);
  comm_dup(&crs->c, c);

  // Let's renumber the ids just in case its the schur solver. Schwarz solver
  // doesn't need re-numbering but we are going to go ahead and do it.
  slong *tid = tcalloc(slong, crs->un);
  for (uint i = 0; i < n; i++)
    tid[i] = id[i];

  // Find the mapping from user ids to unique ids (compressed ids) local to the
  // processor. Compressed vector size is `crs->cn`.
  ulong *uid = tcalloc(ulong, crs->un);
  crs->u2c = tcalloc(sint, crs->un);
  crs->cn = unique_ids(crs->u2c, uid, crs->un, tid, &crs->bfr);
  // printf("id = %u un = %u cn = %u\n", c->id, crs->un, crs->cn);

  // Now renumber unique ids based on whether they are internal or on interface.
  slong *nid = tcalloc(slong, crs->cn);
  number_dofs(nid, crs, tid, uid);
  free(tid), free(uid);

  // Now let's setup the coarse system. Create `struct mij` entries and pass
  // them into schur setup. Which processor owns the dof? All the local dofs
  // are owned by those specific preocessors -- interface dofs are owned in
  // a load balanced manner.
  uint nr = crs->ng[1] / c->np;
  ulong s0 = nr * c->np;
  uint nrem = crs->ng[1] - s0;
  uint p0 = c->np - nrem;
  s0 = p0 * nr;

  struct array mijs;
  array_init(struct mij, &mijs, n);

  struct mij m = {.r = 0, .c = 0, .idx = 0, .p = 0, .v = 0};
  for (uint k = 0; k < nz; k++) {
    sint i = crs->u2c[Ai[k]], j = crs->u2c[Aj[k]];
    if (i < 0 || j < 0 || A[k] == 0)
      continue;
    m.r = nid[i], m.c = nid[j], m.v = A[k], m.p = c->id;
    if (m.r > crs->ng[0]) {
      if (m.r - crs->ng[0] <= s0)
        m.p = (m.r - crs->ng[0] - 1) / nr;
      else
        m.p = p0 + (m.r - crs->ng[0] - s0 - 1) / (nr + 1);
    }
    array_cat(struct mij, &mijs, &m, 1);
  }

  // Now let's assemble the matrix by sending load balancing the interface rows.
  // Assembled size is `an`.
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct mij, &mijs, p, 1, &cr);

  nid = trealloc(slong, nid, crs->cn + crs->n[0] + nr + 1);
  for (uint i = 0; i < crs->cn; i++)
    nid[i] = -nid[i];

  crs->an = 0;
  if (mijs.n > 0) {
    sarray_sort_2(struct mij, mijs.ptr, mijs.n, r, 1, c, 1, &crs->bfr);
    struct mij *pm = (struct mij *)mijs.ptr;
    uint i = 0, j;
    while (i < mijs.n) {
      for (j = i + 1; j < mijs.n && pm[j].r == pm[i].r; j++)
        ;
      nid[crs->cn + crs->an] = pm[i].r, crs->an++, i = j;
    }
  }
  crs->n[1] = crs->an - crs->n[0];
  crs->s[1] = nid[crs->cn + crs->n[0]];
  crs->c2a = gs_setup(nid, crs->cn + crs->an, c, 0, gs_pairwise, 0);
  // printf("id = %u an = %u n[1] = %u nr = %u\n", c->id, crs->an, crs->n[1],
  // nr);

  switch (type) {
  case 0:
    schur_setup(crs, &mijs, &cr, &crs->bfr);
    break;
  default:
    break;
  }

  array_free(&mijs), crystal_free(&cr);

  return crs;
}

void crs_parrsb_solve(scalar *x, struct coarse *crs, scalar *b, scalar tol) {
  metric_init();

  scalar *rhs = tcalloc(scalar, crs->cn + crs->an);
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      rhs[crs->u2c[i]] += b[i];
  }
  gs(rhs, gs_double, gs_add, 1, crs->c2a, &crs->bfr);

  switch (crs->type) {
  case 0:
    schur_solve(rhs + crs->cn, crs, rhs + crs->cn, tol, &crs->bfr);
    break;
  default:
    break;
  }

  gs(rhs, gs_double, gs_add, 0, crs->c2a, &crs->bfr);
  for (uint i = 0; i < crs->un; i++) {
    if (crs->u2c[i] >= 0)
      x[i] = rhs[crs->u2c[i]];
  }
  free(rhs);

  metric_push_level();
  // metric_crs_print(&crs->c, 1);
}

void crs_parrsb_free(struct coarse *crs) {
  if (crs != NULL) {
    switch (crs->type) {
    case 0:
      schur_free(crs);
      break;
    default:
      break;
    }
    if (crs->u2c)
      free(crs->u2c);
    gs_free(crs->c2a);
    comm_free(&crs->c), buffer_free(&crs->bfr);
    free(crs), crs = NULL;
  }
}

#undef MIN
