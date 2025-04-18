#include "metrics.h"
#include "parrsb_impl.h"

#define GS 1
#define CSR 2

struct laplacian {
  int type, nv;
  uint nel;
  void *data;
};

/*
 * Laplacian - CSR based implementation.
 */
struct csr_laplacian {
  struct par_mat *M;
  struct gs_data *gsh;
  scalar *buf;
};

static int par_csr_init(laplacian l, const struct rsb_element *elems,
                        const int nv, const struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);

  struct array nbrs, eij;
  find_nbrs_rsb(&nbrs, elems, l->nel, nv, c, &cr, bfr);
  compress_nbrs(&eij, &nbrs, bfr);

  struct csr_laplacian *L = l->data = tcalloc(struct csr_laplacian, 1);
  struct par_mat *M = L->M = par_csr_setup_ext(&eij, 1, bfr);
  L->gsh = setup_Q(L->M, c, bfr);

  uint nnz = M->rn > 0 ? M->adj_off[M->rn] + M->rn : 0;
  L->buf = tcalloc(scalar, nnz);

  crystal_free(&cr);

  array_free(&nbrs);
  array_free(&eij);

  return 0;
}

static int par_csr(scalar *v, const laplacian l, scalar *u, buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  if (L != NULL) {
    mat_vec_csr(v, u, L->M, L->gsh, L->buf, bfr);
    return 0;
  }
  return 1;
}

static int par_csr_free(laplacian l) {
  if (l->data != NULL) {
    struct csr_laplacian *L = (struct csr_laplacian *)l->data;
    par_mat_free(L->M), gs_free(L->gsh), free(L->buf);
    free(L);
  }
  return 0;
}

/*
 * Laplacian - GS based implementation.
 */
struct gs_laplacian {
  scalar *diag, *u;
  struct gs_data *gsh;
};

static int gs_weighted_init(laplacian l, const struct rsb_element *elems,
                            const unsigned nv, const struct comm *c,
                            buffer *buf) {
  uint lelt = l->nel;
  uint npts = nv * lelt;
  slong *vertices = tcalloc(slong, npts);
  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++) vertices[i * nv + j] = elems[i].vertices[j];

  struct gs_laplacian *gl = l->data = tcalloc(struct gs_laplacian, 1);
  gl->u = tcalloc(scalar, npts);
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++) gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(scalar, lelt);
  for (i = 0; i < lelt; i++) {
    gl->diag[i] = 0.0;
    for (j = 0; j < nv; j++) gl->diag[i] += gl->u[nv * i + j];
  }

  if (vertices != NULL) free(vertices);

  return 0;
}

static int gs_weighted(scalar *v, laplacian l, scalar *u, buffer *bfr) {
  uint lelt = l->nel;
  unsigned nv = l->nv;
  struct gs_laplacian *gl = l->data;

  uint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++) gl->u[nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, bfr);

  for (i = 0; i < lelt; i++) {
    v[i] = gl->diag[i] * u[i];
    for (j = 0; j < nv; j++) v[i] -= gl->u[nv * i + j];
  }

  return 0;
}

static int gs_weighted_free(laplacian l) {
  struct gs_laplacian *gl = l->data;
  if (gl->u != NULL) free(gl->u);
  if (gl->diag != NULL) free(gl->diag);
  gs_free(gl->gsh);
  free(l->data);
  return 0;
}

/*
 * Laplacian - user API.
 */
int laplacian_init(laplacian *l_, const struct array *elements,
                   const element_info ei, const struct comm *c, buffer *buf) {
  metric_tic(c, RSB_LAPLACIAN_SETUP);

  laplacian l = *l_ = tcalloc(struct laplacian, 1);
  l->nv = ei->nv;
  l->type = (l->nv > 0) ? GS : CSR;
  l->nel = elements->n;

  const struct rsb_element *pe = (const struct rsb_element *)elements->ptr;
  switch (l->type) {
  case CSR: par_csr_init(l, pe, l->nv, c, buf); break;
  case GS: gs_weighted_init(l, pe, l->nv, c, buf); break;
  default: return 1; break;
  }

  metric_toc(c, RSB_LAPLACIAN_SETUP);
  return 0;
}

uint laplacian_get_size(laplacian wl) { return wl->nel; }

int laplacian_op(scalar *v, laplacian l, scalar *u, buffer *buf) {
  switch (l->type) {
  case CSR: par_csr(v, l, u, buf); break;
  case GS: gs_weighted(v, l, u, buf); break;
  default: return 1; break;
  }

  return 0;
}

int laplacian_free(laplacian *l_) {
  if (!l_ || !(*l_)) return 1;

  laplacian l = *l_;
  switch (l->type) {
  case GS: gs_weighted_free(l); break;
  case CSR: par_csr_free(l); break;
  default: return 1; break;
  }

  free(l), l = 0;
  return 0;
}

#undef CSR
#undef GS
