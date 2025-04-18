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
  // struct par_mat *M;
  struct gs_data *gsh;
  scalar *buf;
};

static int par_csr_init(laplacian l, struct array *nlist, const struct comm *c,
                        buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  crystal_free(&cr);

  return 0;
}

static int par_csr(scalar *v, const laplacian l, scalar *u, buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  // mat_vec_csr(v, u, L->M, L->gsh, L->buf, bfr);

  return 0;
}

static int par_csr_free(laplacian l) {

  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  gs_free(L->gsh), free(L->buf), free(L);

  return 0;
}

/*
 * Laplacian - GS based implementation.
 */
struct gs_laplacian {
  scalar *diag, *u;
  struct gs_data *gsh;
};

static int gs_weighted_init(laplacian l, struct array *elist,
                            const struct comm *c, buffer *buf) {
  uint lelt = l->nel;
  uint nv = l->nv;
  uint npts = nv * lelt;

  const struct rsb_element *pe = (const struct rsb_element *)elist->ptr;
  slong *vertices = tcalloc(slong, npts);
  for (uint i = 0; i < lelt; i++)
    for (uint j = 0; j < nv; j++) vertices[i * nv + j] = pe[i].vertices[j];

  struct gs_laplacian *gl = l->data = tcalloc(struct gs_laplacian, 1);
  gl->u = tcalloc(scalar, npts);
  for (uint i = 0; i < lelt; i++)
    for (uint j = 0; j < nv; j++) gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(scalar, lelt);
  for (uint i = 0; i < lelt; i++) {
    gl->diag[i] = 0.0;
    for (uint j = 0; j < nv; j++) gl->diag[i] += gl->u[nv * i + j];
  }

  free(vertices);

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
  struct gs_laplacian *gl = (struct gs_laplacian *)l->data;
  if (gl->u != NULL) free(gl->u);
  if (gl->diag != NULL) free(gl->diag);
  gs_free(gl->gsh);
  free(l->data);
  return 0;
}

/*
 * Laplacian - user API.
 */
int laplacian_init(laplacian *l_, const struct array *arr,
                   const element_info ei, const struct comm *c, buffer *buf) {
  metric_tic(c, RSB_LAPLACIAN_SETUP);

  laplacian l = *l_ = tcalloc(struct laplacian, 1);
  l->nv = ei->nv;
  l->type = (l->nv > 0) ? GS : CSR;
  l->nel = arr->n;

  switch (l->type) {
  case CSR: par_csr_init(l, arr, c, buf); break;
  case GS: gs_weighted_init(l, arr, c, buf); break;
  default: return 1; break;
  }

  metric_toc(c, RSB_LAPLACIAN_SETUP);
  return 0;
}

uint laplacian_get_size(laplacian l) {
  if (!l) return 0;
  return l->nel;
}

int laplacian_op(scalar *v, laplacian l, scalar *u, buffer *buf) {
  if (!l || !l->data) return 1;

  switch (l->type) {
  case CSR: par_csr(v, l, u, buf); break;
  case GS: gs_weighted(v, l, u, buf); break;
  default: return 1; break;
  }

  return 0;
}

int laplacian_free(laplacian *l_) {
  if (!l_ || !(*l_) || !((*l_)->data)) return 1;

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
