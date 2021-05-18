#include <genmap-impl.h>
#include <genmap-partition.h>

struct unique_id {
  ulong id;
  int proc;
  int shared;
};

typedef struct {
  int nNeighbors;
  GenmapULong elementId;
  GenmapULong vertexId;
  uint workProc;
} vertex;

static void genmap_find_neighbors(struct array *nbrs, ulong *elem_id,
                                  genmap_handle h, struct comm *cc) {
  metric_tic(cc, FIRSTHALF);
  sint lelt = genmap_get_nel(h);
  int nv = genmap_get_nvertices(h);
  size_t size = lelt * nv;

  slong out[2][1], buf[2][1];
  slong in = lelt;
  comm_scan(out, cc, gs_long, gs_add, &in, 1, buf);
  ulong eid = out[0][0] + 1;

  struct array vertices;
  array_init(vertex, &vertices, size);

  struct rsb_element *elems = genmap_get_elements(h);
  sint i, j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < nv; j++) {
      vertex vrt = {.nNeighbors = 0,
                    .elementId = (elem_id == NULL) ? eid : elem_id[i],
                    .vertexId = elems[i].vertices[j],
                    .workProc = cc->id};
      array_cat(vertex, &vertices, &vrt, 1);
    }
    eid++;
  }
  assert(vertices.n == lelt * nv);

  sarray_sort(vertex, vertices.ptr, vertices.n, vertexId, 1, &h->buf);

  /* Compress vertices and figure out shared ids */
  struct array unique;
  array_init(struct unique_id, &unique, lelt);

  ulong vertex_id = 0; /* vertexId > 0 */
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

  /* Split the vertex ids to shared and local */
  struct array shared;
  array_init(vertex, &shared, lelt);
  struct array local;
  array_init(vertex, &local, lelt);

  /* Now send only the shared ids to workProc */
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

  /* Concatenate local ids and recieved shared ids */
  array_init(vertex, &vertices, 2 * lelt);
  array_cat(vertex, &vertices, local.ptr, local.n);
  array_cat(vertex, &vertices, shared.ptr, shared.n);
  array_free(&local);
  array_free(&shared);

  metric_toc(cc, FIRSTHALF);

  /* Generate CSR entries */
  metric_tic(cc, SECONDHALF);

  sarray_sort(vertex, vertices.ptr, vertices.n, vertexId, 1, &h->buf);
  size = vertices.n;
  vp = vertices.ptr;

  array_init(csr_entry, nbrs, 10);

  s = 0;
  csr_entry t;
  t.v = -1.0;
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
        array_cat(csr_entry, nbrs, &t, 1);
      }
    }
    s = e;
  }

  sarray_transfer(csr_entry, nbrs, proc, 1, &cr);

  array_free(&vertices);
  crystal_free(&cr);

  metric_toc(cc, SECONDHALF);
}

static int genmap_unweighted_laplacian_csr_init(struct array *csr,
                                                csr_entry *ptr, uint pn,
                                                genmap_handle h,
                                                struct comm *c) {
  uint nr;
  sint diag;
  uint e = 0, e1, s;
  while (e < pn) {
    s = e;
    nr = 0;
    diag = -1;
    while (e < pn && ptr[s].r == ptr[e].r) {
      for (e1 = e; e < pn && ptr[e1].r == ptr[e].r && ptr[e1].c == ptr[e].c;
           e++)
        ;
      array_cat(csr_entry, csr, &ptr[e1], 1);
      nr++;
      if (ptr[e1].r == ptr[e1].c)
        diag = csr->n - 1;
    }
    assert(diag >= 0);
    ((csr_entry *)csr->ptr)[diag].v = nr - 1.0;
  }

  return 0;
}

static int genmap_weighted_laplacian_csr_init(struct array *csr, csr_entry *ptr,
                                              uint pn, genmap_handle h,
                                              struct comm *c) {
  GenmapScalar v, v1;
  sint diag;
  uint e = 0, e1, s;
  while (e < pn) {
    s = e;
    v = 0.0;
    diag = -1;
    while (e < pn && ptr[s].r == ptr[e].r) {
      v1 = 0.0;
      for (e1 = e; e < pn && ptr[e1].r == ptr[e].r && ptr[e1].c == ptr[e].c;
           e++)
        v1 += ptr[e].v;
      ptr[e1].v = v1;
      array_cat(csr_entry, csr, &ptr[e1], 1);
      v += v1;
      if (ptr[e1].r == ptr[e1].c)
        diag = csr->n - 1;
    }
    assert(diag >= 0);
    ((csr_entry *)csr->ptr)[diag].v = -h->nv - v;
  }

  return 0;
}

int genmap_laplacian_csr_init(struct csr_mat_ **M, ulong *elem_id,
                              genmap_handle h, struct comm *c) {
  metric_tic(c, LAPLACIANSETUP);

  metric_tic(c, FINDNBRS);
  struct array entries;
  genmap_find_neighbors(&entries, elem_id, h, c);
  metric_toc(c, FINDNBRS);

  sarray_sort_2(csr_entry, entries.ptr, entries.n, r, 1, c, 1, &h->buf);

  struct array csr;
  array_init(csr_entry, &csr, 10);

  if (h->options->rsb_laplacian_weighted == 0)
    genmap_unweighted_laplacian_csr_init(&csr, entries.ptr, entries.n, h, c);
  else
    genmap_weighted_laplacian_csr_init(&csr, entries.ptr, entries.n, h, c);

  sarray_sort_2(csr_entry, csr.ptr, csr.n, r, 1, c, 1, &h->buf);

  metric_tic(c, CSRMATSETUP);
  csr_mat_setup(M, &csr, c, &h->buf);
  metric_toc(c, CSRMATSETUP);

  array_free(&entries);
  array_free(&csr);

  metric_toc(c, LAPLACIANSETUP);
}

int genmap_laplacian_csr(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  csr_mat_apply(v, h->M, u, &h->buf);
  return 0;
}

#undef SWAP
#undef SORT3
