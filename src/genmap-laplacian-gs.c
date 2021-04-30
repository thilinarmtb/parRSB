#include <genmap-impl.h>
#include <genmap-partition.h>

#define SWAP(v, i, j)                                                          \
  do {                                                                         \
    slong tmp = v[i];                                                          \
    v[i] = v[j];                                                               \
    v[j] = tmp;                                                                \
  } while (0)

#define SORT4(v)                                                               \
  do {                                                                         \
    if (v[0] > v[1])                                                           \
      SWAP(v, 0, 1);                                                           \
    if (v[0] > v[2])                                                           \
      SWAP(v, 0, 2);                                                           \
    if (v[0] > v[3])                                                           \
      SWAP(v, 0, 3);                                                           \
    if (v[1] > v[2])                                                           \
      SWAP(v, 1, 2);                                                           \
    if (v[1] > v[3])                                                           \
      SWAP(v, 1, 3);                                                           \
    if (v[2] > v[3])                                                           \
      SWAP(v, 2, 3);                                                           \
  } while (0)

void genmap_number_faces_and_edges(genmap_handle h, struct comm *c) {
  int nv = genmap_get_nvertices(h);
  int ndim = (nv == 8) ? 3 : 2;
  int nf = 2 * ndim;                      // number of faces
  int nfv = 2 * (ndim - 1);               // number of face vertices
  int ne = (ndim == 3) ? nv + nf - 2 : 0; // number of edges

  sint lelt = genmap_get_nel(h);
  uint nlocal = nv + nf + ne;

  struct rsb_element *elements = genmap_get_elements(h);

  // Copy the vertices first and find the max_id along with it
  slong max_id = 0, buf;
  sint e, i;
  for (e = 0; e < lelt; e++) {
    for (i = 0; i < nv; i++)
      if (elements[e].vertices[i] > max_id)
        max_id = elements[e].vertices[i];
  }

  comm_allreduce(c, gs_long, gs_max, &max_id, 1, &buf);

  slong v[4];
  int j;
  for (e = 0; e < lelt; e++) {
    for (i = 0; i < ne; i++) {
      // We are going to number each edge uniquely.
      // If edge_i = (v1, v2), v_min = GENMAP_MIN(v1, v2), v_max =
      // GENMAP_MAX(v1, v2) and max_id = maximum global id unique id of edge_i =
      // max_id + v_min + (v_max - 1) * max_id
      v[0] = elements[e].vertices[edges3D[i][0]];
      v[1] = elements[e].vertices[edges3D[i][1]];
      elements[e].vertices[nv + i] = max_id + GENMAP_MIN(v[0], v[1]) +
                                     (GENMAP_MAX(v[0], v[1]) - 1) * max_id;
    }

    for (i = 0; i < nf; i++) {
      for (j = 0; j < nfv; j++)
        v[j] = elements[e].vertices[faces3D[i][j] - 1];

      SORT4(v);

      // Two edges or three vertices uniquely define a face
      elements[e].vertices[nv + ne + i] = max_id + max_id * max_id + v[0] +
                                          max_id * (v[1] - 1) +
                                          max_id * max_id * (v[2] - 1);
    }
  }
}

// FIXME: Only works for 3D as of now
static int genmap_unweighted_laplacian_gs_init(genmap_handle h,
                                               struct comm *c) {
  int nv = genmap_get_nvertices(h);
  int ndim = (nv == 8) ? 3 : 2;
  int nf = 2 * ndim;                      // number of faces
  int ne = (ndim == 3) ? nv + nf - 2 : 0; // number of edges

  sint lelt = genmap_get_nel(h);
  uint nlocal = nv + nf + ne;
  uint size = lelt * nlocal;

  // Store vertices, edges and faces
  GenmapLong *vef;
  GenmapMalloc(size, &vef);

  struct rsb_element *elements = genmap_get_elements(h);
  int e, i;
  for (e = 0; e < lelt; e++) {
    for (i = 0; i < nlocal; i++)
      vef[e * nlocal + i] = elements[e].vertices[i];
  }

  if (h->gs != NULL)
    gs_free(h->gs);

  h->gs = gs_setup(vef, size, c, 0, gs_crystal_router, 0);

  GenmapScalar *u;
  GenmapMalloc(size, &u);
  for (i = 0; i < size; i++)
    u[i] = 1.0;

  gs(u, gs_double, gs_add, 0, h->gs, &h->buf);

  GenmapRealloc(lelt, &h->diagonal);
  for (e = 0; e < lelt; e++) {
    h->diagonal[e] = -2.0;
    for (i = 0; i < nv; i++)
      h->diagonal[e] += u[e * nlocal + i];
    for (i = 0; i < ne; i++)
      h->diagonal[e] -= u[e * nlocal + nv + i];
    for (i = 0; i < nf; i++)
      h->diagonal[e] += u[e * nlocal + nv + ne + i];
  }

  GenmapFree(u);
  GenmapFree(vef);

  return 0;
}

static int genmap_unweighted_laplacian_gs(genmap_handle h, GenmapScalar *u,
                                          GenmapScalar *v) {
  int nv = genmap_get_nvertices(h);
  int ndim = (nv == 8) ? 3 : 2;
  int nf = 2 * ndim;                      // number of faces
  int nfv = 2 * (ndim - 1);               // number of face vertices
  int ne = (ndim == 3) ? nv + nf - 2 : 0; // number of edges
  uint nlocal = nv + nf + ne;

  sint lelt = genmap_get_nel(h);
  uint size = lelt * nlocal;

  GenmapScalar *ucv;
  GenmapMalloc(size, &ucv);

  sint i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nlocal; j++)
      ucv[i * nlocal + j] = u[i];

  gs(ucv, gs_double, gs_add, 0, h->gs, &h->buf);

  for (i = 0; i < lelt; i++) {
    v[i] = (h->diagonal[i] + 2.0) * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= ucv[nlocal * i + j];
    for (j = 0; j < ne; j++)
      v[i] += ucv[nlocal * i + nv + j];
    for (j = 0; j < nf; j++)
      v[i] -= ucv[nlocal * i + nv + ne + j];
  }

  GenmapFree(ucv);

  return 0;
}

static int genmap_weighted_laplacian_gs_init(genmap_handle h, struct comm *c) {
  GenmapInt lelt = genmap_get_nel(h);
  GenmapInt nv = genmap_get_nvertices(h);

  GenmapRealloc(lelt, &h->diagonal);
  GenmapUInt numPoints = (GenmapUInt)nv * lelt;

  GenmapLong *vertices;
  GenmapMalloc(numPoints, &vertices);

  struct rsb_element *elements = genmap_get_elements(h);
  GenmapInt i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++)
      vertices[i * nv + j] = elements[i].vertices[j];
  }

  if (h->gs != NULL)
    gs_free(h->gs);

  h->gs = gs_setup(vertices, numPoints, c, 0, gs_crystal_router, 0);

  GenmapScalar *u;
  GenmapMalloc(numPoints, &u);

  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      u[nv * i + j] = 1.;

  gs(u, gs_double, gs_add, 0, h->gs, &h->buf);

  for (i = 0; i < lelt; i++) {
    h->diagonal[i] = -nv;
    for (j = 0; j < nv; j++)
      h->diagonal[i] += u[nv * i + j];
  }

  GenmapFree(u);
  GenmapFree(vertices);

  return 0;
}

static int genmap_weighted_laplacian_gs(genmap_handle h, GenmapScalar *u,
                                        GenmapScalar *v) {
  GenmapInt lelt = genmap_get_nel(h);
  GenmapInt nv = genmap_get_nvertices(h);

  GenmapScalar *ucv;
  GenmapMalloc((size_t)(nv * lelt), &ucv);

  GenmapInt i, j;
  for (i = 0; i < lelt; i++)
    for (j = 0; j < nv; j++)
      ucv[nv * i + j] = u[i];

  gs(ucv, gs_double, gs_add, 0, h->gs, &h->buf);

  for (i = 0; i < lelt; i++) {
    v[i] = (h->diagonal[i] + nv) * u[i];
    for (j = 0; j < nv; j++)
      v[i] -= ucv[nv * i + j];
  }

  GenmapFree(ucv);

  return 0;
}

int genmap_laplacian_gs_init(genmap_handle h, struct comm *c) {
  metric_tic(c, LAPLACIANSETUP);
  if (h->options->rsb_laplacian_weighted == 0)
    genmap_unweighted_laplacian_gs_init(h, c);
  else
    genmap_weighted_laplacian_gs_init(h, c);
  metric_toc(c, LAPLACIANSETUP);
}

int genmap_laplacian_gs(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  if (h->options->rsb_laplacian_weighted == 0)
    genmap_unweighted_laplacian_gs(h, u, v);
  else
    genmap_weighted_laplacian_gs(h, u, v);
}
