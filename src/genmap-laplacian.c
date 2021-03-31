#include <genmap-impl.h>

struct unique_id {
  ulong id;
  int proc;
  int shared;
};

typedef struct {
  GenmapULong sequenceId;
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

static void genmap_find_neighbors(struct array *nbrs, genmap_handle h,
                                  struct comm *cc) {
  metric_tic(cc, FIRSTHALF);
  sint lelt = genmap_get_nel(h);
  sint nv = genmap_get_nvertices(h);

  genmap_comm_scan(h, cc);
  ulong elem_id = genmap_get_local_start_index(h) + 1;
  ulong sequenceId = elem_id * nv;

  size_t size = lelt * nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  struct rsb_element *elems = genmap_get_elements(h);
  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.sequenceId = sequenceId,
                    .nNeighbors = 0,
                    .elementId = elem_id,
                    .vertexId = elems[i].vertices[j],
                    .workProc = cc->id};
      array_cat(vertex, &vertices, &vrt, 1);
      sequenceId++;
    }
    elem_id++;
  }
  assert(vertices.n == lelt * nv);

  // Compress vertices and figure out shared ids
  sarray_sort(vertex, vertices.ptr, vertices.n, vertexId, 1, &h->buf);

  struct array unique;
  array_init(struct unique_id, &unique, lelt);

  ulong vertex_id = 0; // vertexId > 0
  vertex *vp = vertices.ptr;
  for (i = 0; i < vertices.n; i++) {
    if (vp[i].vertexId > vertex_id) {
      vertex_id = vp[i].vertexId;
      struct unique_id id = {
          .proc = vertex_id % cc->np, .id = vertex_id, .shared = 0};
      array_cat(struct unique_id, &unique, &id, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, cc);

  sarray_transfer(struct unique_id, &unique, proc, 1, &cr);
  sarray_sort(struct unique_id, unique.ptr, unique.n, id, 1, &h->buf);

  struct unique_id *up = unique.ptr;
  sint s = 0, e;
  while (s < unique.n) {
    e = s + 1;
    while (e < unique.n && up[e].id == up[s].id)
      e++;
    if (e > s + 1)
      while (s < e) {
        up[s].shared = 1;
        s++;
      }
    s = e;
  }

  sarray_transfer(struct unique_id, &unique, proc, 1, &cr);
  sarray_sort(struct unique_id, unique.ptr, unique.n, id, 1, &h->buf);

  struct array shared;
  array_init(vertex, &shared, lelt);
  struct array local;
  array_init(vertex, &local, lelt);

  // Now send only the shared ids to workProc
  up = unique.ptr;
  int workProc;
  s = 0;
  for (i = 0; i < unique.n; i++) {
    vertex_id = up[i].id;
    if (up[i].shared == 1) {
      workProc = vertex_id % cc->np;
      while (s < vertices.n && vp[s].vertexId == vertex_id) {
        vp[s].workProc = workProc;
        array_cat(vertex, &shared, &vp[s], 1);
        s++;
      }
    } else {
      workProc = cc->id;
      while (s < vertices.n && vp[s].vertexId == vertex_id) {
        vp[s].workProc = workProc;
        array_cat(vertex, &local, &vp[s], 1);
        s++;
      }
    }
  }
  array_free(&unique);
  array_free(&vertices);

  sarray_transfer(vertex, &shared, workProc, 1, &cr);

  array_init(vertex, &vertices, 2 * lelt);
  array_cat(vertex, &vertices, local.ptr, local.n);
  array_cat(vertex, &vertices, shared.ptr, shared.n);
  array_free(&local);
  array_free(&shared);

  metric_toc(cc, FIRSTHALF);

  metric_tic(cc, SECONDHALF);

  sarray_sort(vertex, vertices.ptr, vertices.n, vertexId, 1, &h->buf);
  size = vertices.n;
  vp = vertices.ptr;

  struct array csr;
  array_init(csr_entry, &csr, 10);

  // FIXME: Assumes quads or hexes
  s = 0;
  csr_entry t;
  while (s < size) {
    e = s + 1;
    while (e < size && vp[s].vertexId == vp[e].vertexId)
      e++;
    int n_neighbors = GENMAP_MIN(e, size) - s;

    for (i = s; i < GENMAP_MIN(e, size); i++) {
      t.r = vp[i].elementId;
      t.proc = vp[i].workProc;
      for (j = 0; j < n_neighbors; j++) {
        t.c = vp[s + j].elementId;
        array_cat(csr_entry, &csr, &t, 1);
      }
    }
    s = e;
  }
  array_free(&vertices);

  sarray_transfer(csr_entry, &csr, proc, 1, &cr);
  sarray_sort_2(csr_entry, csr.ptr, csr.n, r, 1, c, 1, &h->buf);

  metric_toc(cc, SECONDHALF);

  array_init(entry, nbrs, lelt);

  if (csr.n > 0) {

    entry ee = {0, 0, 0, 0, 0, 0.0}, ep = {0, 0, 0, 0, 0.0};
    csr_entry *ap = csr.ptr;
    ep.r = ap[0].r;
    ep.c = ap[0].c;

    array_cat(entry, nbrs, &ep, 1);

    for (i = 1; i < csr.n; i++) {
      ee.r = ap[i].r;
      ee.c = ap[i].c;
      if (ee.r != ep.r || ee.c != ep.c) {
        array_cat(entry, nbrs, &ee, 1);
        ep = ee;
      }
    }

    sarray_sort_2(entry, nbrs->ptr, nbrs->n, r, 1, c, 1, &h->buf);
  }

  array_free(&csr);
  crystal_free(&cr);
}

int GenmapInitLaplacian(genmap_handle h, struct comm *c) {
  struct array entries;

  metric_tic(c, FINDNBRS);
  genmap_find_neighbors(&entries, h, c);
  metric_toc(c, FINDNBRS);

  metric_tic(c, CSRMATSETUP);
  csr_mat_setup(&entries, c, &h->M);
  metric_toc(c, CSRMATSETUP);

  array_free(&entries);

  metric_tic(c, CSRTOPSETUP);
  h->gsh = get_csr_top(h->M, c);
  metric_toc(c, CSRTOPSETUP);

  GenmapRealloc(h->M->row_off[h->M->rn], &h->b);

#if 0
  int nnz = h->M->row_off[h->M->rn];
  double fro[2] = {0.0, 0.0}, buf[2];
  for (int i = 0; i < nnz; i++) {
    fro[0] += h->M->v[i];
    fro[1] += h->M->v[i] * h->M->v[i];
  }
  comm_allreduce(c, gs_double, gs_add, &fro, 2, &buf);
  if (c->id == 0)
    printf("nrom(G,'1')=%g\nnorm(G,'fro')=%g\n", fro[0], fro[1]);
#endif

  return 0;
}

int GenmapLaplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  csr_mat_gather(h->M, h->gsh, u, h->b, &h->buf);
  csr_mat_apply(v, h->M, h->b);

  return 0;
}

#undef SWAP
#undef SORT3
