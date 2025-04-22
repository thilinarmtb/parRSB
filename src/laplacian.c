#include "metrics.h"
#include "parrsb_impl.h"

#define GS 1
#define CSR 2

struct laplacian {
  int type, nv;
  uint n;
  void *data;
};

/*
 * Laplacian - CSR based implementation.
 */

struct csr_laplacian {
  uint n;
  uint *di, *off;
  scalar *v;
  struct gs_data *gsh;
};

static int csr_init(laplacian l, const struct array *nlist,
                    const struct comm *c, buffer *bfr) {
  sarray_sort_2(struct graph_element, nlist->ptr, nlist->n, u, 1, v, 1, bfr);

  const graph_element pe = (const graph_element)nlist->ptr;
  uint rn = (nlist->n > 0);
  for (uint n = 1; n < nlist->n; n++)
    if (pe[n].u != pe[n - 1].u) rn++;

  struct csr_laplacian *L = l->data = tcalloc(struct csr_laplacian, 1);
  L->n = rn;
  L->di = tcalloc(uint, L->n);
  L->off = tcalloc(uint, L->n + 1);
  L->v = tcalloc(scalar, nlist->n + L->n);

  buffer_reserve(bfr, (nlist->n + L->n) * sizeof(slong));

  slong *ids = (slong *)bfr->ptr;
  uint n = 0;
  rn = 0;
  while (n < nlist->n) {
    uint i = n;
    while (i < n && pe[n].u == pe[i].u && (pe[i].v < pe[i].u))
      L->v[i] = -1, ids[i] = -pe[i].v, i++;

    uint d = L->di[rn] = i;

    while (i < n && pe[n].u == pe[i].u) L->v[i + 1] = -1, ids[i + 1] = -pe[i].v;

    ids[d] = pe[n].u, L->v[d] = i - n, L->off[++rn] = n = i;
  }

  L->gsh = gs_setup(ids, nlist->n + L->n, c, 0, gs_crystal_router, 0);

  return 0;
}

static int csr_op(scalar *v, const laplacian l, scalar *u, buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  return 0;
}

static int csr_free(laplacian l) {
  if (!l) return 1;

  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  if (!L) return 1;

  free(L->di), free(L->off), free(L->v), gs_free(L->gsh);
  free(L), l->data = 0;

  return 0;
}

/*
 * Laplacian - GS based implementation.
 */
struct gs_laplacian {
  scalar *diag, *u;
  struct gs_data *gsh;
};

static int gs_weighted_init(laplacian l, const struct array *elist,
                            const struct comm *c, buffer *buf) {
  uint ne = l->n;
  uint nv = l->nv;
  uint npts = nv * ne;

  const struct rsb_element *pe = (const struct rsb_element *)elist->ptr;
  slong *vertices = tcalloc(slong, npts);
  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) vertices[i * nv + j] = pe[i].vertices[j];

  struct gs_laplacian *gl = l->data = tcalloc(struct gs_laplacian, 1);
  gl->u = tcalloc(scalar, npts);
  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) gl->u[nv * i + j] = 1.0;

  gl->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);
  gs(gl->u, gs_double, gs_add, 0, gl->gsh, buf);

  gl->diag = tcalloc(scalar, ne);
  for (uint i = 0; i < ne; i++) {
    gl->diag[i] = 0.0;
    for (uint j = 0; j < nv; j++) gl->diag[i] += gl->u[nv * i + j];
  }

  free(vertices);

  return 0;
}

static int gs_weighted_op(scalar *v, laplacian l, scalar *u, buffer *bfr) {
  uint ne = l->n;
  unsigned nv = l->nv;
  struct gs_laplacian *gl = l->data;

  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) gl->u[nv * i + j] = u[i];

  gs(gl->u, gs_double, gs_add, 0, gl->gsh, bfr);

  for (uint i = 0; i < ne; i++) {
    v[i] = gl->diag[i] * u[i];
    for (uint j = 0; j < nv; j++) v[i] -= gl->u[nv * i + j];
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
  laplacian l = *l_ = tcalloc(struct laplacian, 1);
  l->nv = ei->nv;
  l->type = (l->nv > 0) ? GS : CSR;
  l->n = arr->n;

  switch (l->type) {
  case CSR: csr_init(l, arr, c, buf); break;
  case GS: gs_weighted_init(l, arr, c, buf); break;
  default: break;
  }

  return 0;
}

uint laplacian_get_size(laplacian l) {
  if (!l) return 0;
  return l->n;
}

int laplacian_op(scalar *v, laplacian l, scalar *u, buffer *buf) {
  if (!l || !l->data) return 1;

  switch (l->type) {
  case CSR: csr_op(v, l, u, buf); break;
  case GS: gs_weighted_op(v, l, u, buf); break;
  default: break;
  }

  return 0;
}

int laplacian_free(laplacian *l_) {
  if (!l_ || !(*l_) || !((*l_)->data)) return 1;

  laplacian l = *l_;
  switch (l->type) {
  case GS: gs_weighted_free(l); break;
  case CSR: csr_free(l); break;
  default: break;
  }

  free(l), l = 0;
  return 0;
}

#undef CSR
#undef GS
