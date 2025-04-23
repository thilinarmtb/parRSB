#include "metrics.h"
#include "parrsb_impl.h"

#define GS 1
#define CSR 2

struct laplacian {
  int type, nv;
  void *data;
};

/*
 * Laplacian - CSR based implementation.
 */

struct csr_laplacian {
  uint n;
  uint *di, *off;
  scalar *v, *wrk;
  struct gs_data *gsh;
};

static void filter_entries(struct array *nlist, const struct comm *c) {
  if (nlist->n == 0) return;

  // Find the minimum and maximum row ids.
  slong max = LONG_MIN, min = LONG_MAX;
  const graph_element pe = (const graph_element)nlist->ptr;
  for (uint i = 0; i < nlist->n; i++) {
    if (pe[i].u > (ulong)max) max = pe[i].u;
    if (pe[i].u < (ulong)min) min = pe[i].u;
  }

  slong wrk;
  comm_allreduce(c, gs_long, gs_max, &max, 1, &wrk);
  comm_allreduce(c, gs_long, gs_min, &min, 1, &wrk);

  // Filter out all the entires `< min` and `> max`.
  uint n = 0;
  for (uint i = 0; i < nlist->n; i++)
    if (pe[i].v >= (ulong)min || pe[i].v <= (ulong)max) pe[n++] = pe[i];
  nlist->n = n;
}

static void csr_init(laplacian l, struct array *nlist, const struct comm *c,
                     buffer *bfr) {
  filter_entries(nlist, c);

  const graph_element pe = (const graph_element)nlist->ptr;
  uint rn = (nlist->n > 0);
  for (uint n = 1; n < nlist->n; n++)
    if (pe[n].u != pe[n - 1].u) rn++;

  struct csr_laplacian *L = l->data = tcalloc(struct csr_laplacian, 1);
  L->n = rn;
  L->di = tcalloc(uint, L->n);
  L->off = tcalloc(uint, L->n + 1);

  uint nnz = nlist->n + L->n;
  L->v = tcalloc(scalar, nnz);
  L->wrk = tcalloc(scalar, nnz);

  buffer_reserve(bfr, nnz * sizeof(slong));

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

  L->gsh = gs_setup(ids, nnz, c, 0, gs_crystal_router, 0);
}

static uint csr_size(laplacian l) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  return L->n;
}

static void csr_op(scalar *v, const laplacian l, scalar *u, buffer *bfr) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;

  for (uint i = 0; i < L->n; i++) L->wrk[L->di[i]] = u[i];

  gs(L->wrk, gs_scalar, gs_add, 0, L->gsh, bfr);

  for (uint i = 0; i < L->n; i++) {
    scalar s = 0;
    for (uint j = L->off[i], je = L->off[i + 1]; j < je; j++)
      s += L->wrk[j] * L->v[j];
    v[i] = s;
  }
}

static void csr_free(laplacian l) {
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;

  gs_free(L->gsh);
  free(L->di), free(L->off), free(L->v), free(L->wrk);
  free(L), l->data = 0;
}

/*
 * Laplacian - GS based implementation.
 */
struct gs_laplacian {
  uint n;
  scalar *diag, *wrk;
  struct gs_data *gsh;
};

static void gs_weighted_init(laplacian l, const struct array *elist,
                             const struct comm *c, buffer *bfr) {
  uint ne = elist->n;
  uint nv = l->nv;
  uint npts = nv * ne;

  buffer_reserve(bfr, npts * sizeof(slong));

  const struct rsb_element *pe = (const struct rsb_element *)elist->ptr;
  slong *vertices = (slong *)bfr->ptr;
  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) vertices[i * nv + j] = pe[i].vertices[j];

  struct gs_laplacian *L = l->data = tcalloc(struct gs_laplacian, 1);
  L->n = ne;
  L->gsh = gs_setup(vertices, npts, c, 0, gs_crystal_router, 0);

  L->wrk = tcalloc(scalar, npts);
  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) L->wrk[nv * i + j] = 1.0;
  gs(L->wrk, gs_scalar, gs_add, 0, L->gsh, bfr);

  L->diag = tcalloc(scalar, ne);
  for (uint i = 0; i < ne; i++) {
    L->diag[i] = 0.0;
    for (uint j = 0; j < nv; j++) L->diag[i] += L->wrk[nv * i + j];
  }
}

static uint gs_weighted_size(laplacian l) {
  struct gs_laplacian *L = (struct gs_laplacian *)l->data;
  return L->n;
}

static void gs_weighted_op(scalar *v, laplacian l, scalar *u, buffer *bfr) {
  unsigned nv = l->nv;

  struct gs_laplacian *L = l->data;
  uint ne = L->n;
  for (uint i = 0; i < ne; i++)
    for (uint j = 0; j < nv; j++) L->wrk[nv * i + j] = u[i];

  gs(L->wrk, gs_scalar, gs_add, 0, L->gsh, bfr);

  for (uint i = 0; i < ne; i++) {
    v[i] = L->diag[i] * u[i];
    for (uint j = 0; j < nv; j++) v[i] -= L->wrk[nv * i + j];
  }
}

static void gs_weighted_free(laplacian l) {
  struct gs_laplacian *L = (struct gs_laplacian *)l->data;
  gs_free(L->gsh);
  free(L->wrk), free(L->diag);
  free(L), l->data = 0;
}

/*
 * Laplacian - user API.
 */
int laplacian_init(laplacian *l_, struct array *arr, const element_info ei,
                   const struct comm *c, buffer *buf) {
  laplacian l = *l_ = tcalloc(struct laplacian, 1);
  l->nv = ei->nv;
  l->type = (l->nv > 0) ? GS : CSR;

  switch (l->type) {
  case CSR: csr_init(l, arr, c, buf); break;
  case GS: gs_weighted_init(l, arr, c, buf); break;
  default: break;
  }

  return 0;
}

uint laplacian_get_size(laplacian l) {
  if (!l || !l->data) return 0;

  switch (l->type) {
  case CSR: return csr_size(l); break;
  case GS: return gs_weighted_size(l); break;
  default: break;
  }

  return 0;
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

  free(l), *l_ = 0;
  return 0;
}

#undef CSR
#undef GS
