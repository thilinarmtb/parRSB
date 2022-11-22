#include "mat.h"
#include "multigrid.h"
#include "parrsb-impl.h"

struct laplacian {
  int type, nv;
  uint nel;
  void *data;
};

//------------------------------------------------------------------------------
// Laplacian - as a `struct par_mat` in CSR mat
//
struct csr_laplacian {
  struct par_mat *M;
  struct gs_data *gsh;
  scalar *buf;
};

static void find_nbrs_rsb(struct array *arr, const struct array *elems,
                          unsigned nv, const struct comm *c, struct crystal *cr,
                          buffer *bfr) {
  size_t nelt = elems->n, npts = nelt * nv;

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong eid = out[0][0] + 1;

  struct array vertices;
  array_init(struct nbr, &vertices, npts);

  struct nbr v;
  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  for (uint i = 0; i < nelt; i++) {
    v.r = eid++;
    for (unsigned j = 0; j < nv; j++) {
      v.c = pe[i].vertices[j], v.proc = v.c % c->np;
      array_cat(struct nbr, &vertices, &v, 1);
    }
  }

  sarray_transfer(struct nbr, &vertices, proc, 1, cr);
  sarray_sort(struct nbr, vertices.ptr, vertices.n, c, 1, bfr);

  array_init(struct nbr, arr, vertices.n * 10);

  uint s = 0, e;
  struct nbr *pv = (struct nbr *)vertices.ptr;
  while (s < vertices.n) {
    e = s + 1;
    while (e < vertices.n && pv[s].c == pv[e].c)
      e++;
    for (uint i = s; i < e; i++) {
      v = pv[i];
      for (unsigned j = s; j < e; j++) {
        v.c = pv[j].r;
        array_cat(struct nbr, arr, &v, 1);
      }
    }
    s = e;
  }
  array_free(&vertices);

  sarray_transfer(struct nbr, arr, proc, 1, cr);
}

static int par_csr_init(struct laplacian *l, const struct array *elems,
                        unsigned nv, const struct comm *c, buffer *bfr) {
  l->data = tcalloc(struct csr_laplacian, 1);
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;

  struct array nbrs, eij;

  struct crystal cr;
  crystal_init(&cr, c);
  find_nbrs_rsb(&nbrs, elems, nv, c, &cr, bfr);
  crystal_free(&cr);

  compress_nbrs(&eij, &nbrs, bfr);
  struct par_mat *M = L->M = tcalloc(struct par_mat, 1);
  par_csr_setup(M, &eij, 1, bfr);
  array_free(&nbrs), array_free(&eij);

  L->gsh = setup_Q(M, c, bfr);

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  L->buf = tcalloc(scalar, nnz);

  return 0;
}

static int par_csr(scalar *v, const struct laplacian *l, scalar *u,
                   buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  if (L) {
    par_mat_vec(v, u, L->M, L->gsh, L->buf, bfr);
    return 0;
  }
  return 1;
}

static int par_csr_free(struct laplacian *l) {
  if (l && l->data) {
    struct csr_laplacian *L = (struct csr_laplacian *)l->data;
    par_mat_free(L->M), gs_free(L->gsh), tfree(L->buf);
    tfree(L);
    return 0;
  }
  return 1;
}

//------------------------------------------------------------------------------
// Laplacian - GS
//
struct gs_laplacian {
  scalar *diag, *u;
  struct gs_data *gsh;
};

static int gs_weighted_init(struct laplacian *l, const struct array *elems,
                            unsigned nv, struct comm *c, buffer *buf) {
  l->data = tcalloc(struct gs_laplacian, 1);
  struct gs_laplacian *gl = (struct gs_laplacian *)l->data;

  size_t nelt = elems->n, npts = nv * nelt;
  slong *vertices = tcalloc(slong, npts);

  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  for (uint i = 0; i < nelt; i++) {
    for (unsigned j = 0; j < nv; j++)
      vertices[i * nv + j] = pe[i].vertices[j];
  }

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_pairwise, 0);

  gl->u = tcalloc(scalar, npts);
  for (uint i = 0; i < nelt; i++) {
    for (unsigned j = 0; j < nv; j++)
      gl->u[nv * i + j] = 1.0;
  }

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(scalar, nelt);
  for (uint i = 0; i < nelt; i++) {
    gl->diag[i] = 0.0;
    for (unsigned j = 0; j < nv; j++)
      gl->diag[i] += gl->u[nv * i + j];
  }

  tfree(vertices);

  return 0;
}

static int gs_weighted(scalar *v, struct laplacian *l, scalar *u, buffer *buf) {
  struct gs_laplacian *gl = (struct gs_laplacian *)l->data;

  for (uint i = 0; i < l->nel; i++)
    for (unsigned j = 0; j < l->nv; j++)
      gl->u[l->nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  for (uint i = 0; i < l->nel; i++) {
    v[i] = gl->diag[i] * u[i];
    for (unsigned j = 0; j < l->nv; j++)
      v[i] -= gl->u[l->nv * i + j];
  }

  return 0;
}

static int gs_weighted_free(struct laplacian *l) {
  if (l && l->data) {
    struct gs_laplacian *gl = (struct gs_laplacian *)l->data;
    tfree(gl->u), tfree(gl->diag), gs_free(gl->gsh);
    tfree(gl);
    return 0;
  }
  return 1;
}

//------------------------------------------------------------------------------
// Laplacian
//
struct laplacian *laplacian_init(struct array *elems, unsigned nv,
                                 unsigned type, struct comm *c, buffer *buf) {
  struct laplacian *l = tcalloc(struct laplacian, 1);
  l->type = type, l->nv = nv, l->nel = elems->n;

  if (type & CSR)
    par_csr_init(l, elems, nv, c, buf);
  else if (type & GS)
    gs_weighted_init(l, elems, nv, c, buf);
  return l;
}

int laplacian(scalar *v, struct laplacian *l, scalar *u, buffer *buf) {
  if (l->type & CSR)
    par_csr(v, l, u, buf);
  else if (l->type & GS)
    gs_weighted(v, l, u, buf);
  else
    return 1;
  return 0;
}

void laplacian_free(struct laplacian *l) {
  if (l) {
    if (l->type & CSR)
      par_csr_free(l);
    else if (l->type & GS)
      gs_weighted_free(l);
    free(l);
  }
}
