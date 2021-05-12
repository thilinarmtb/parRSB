#include <genmap-impl.h>
#include <genmap-partition.h>

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

  sarray_sort(vertex, vertices.ptr, vertices.n, vertexId, 1, &h->buf);

  // Compress vertices and figure out shared ids
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

  // Split the vertex ids to shared and local
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

  // Concatenate local ids and recieved shared ids
  array_init(vertex, &vertices, 2 * lelt);
  array_cat(vertex, &vertices, local.ptr, local.n);
  array_cat(vertex, &vertices, shared.ptr, shared.n);
  array_free(&local);
  array_free(&shared);

  metric_toc(cc, FIRSTHALF);

  // Generate CSR entries
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
  array_free(&vertices);

  sarray_transfer(csr_entry, nbrs, proc, 1, &cr);

  metric_toc(cc, SECONDHALF);

  crystal_free(&cr);
}

static int genmap_unweighted_laplacian_csr_init(genmap_handle h,
                                                struct comm *c) {
  metric_tic(c, FINDNBRS);
  struct array entries;
  genmap_find_neighbors(&entries, h, c);
  metric_toc(c, FINDNBRS);

  sarray_sort_2(csr_entry, entries.ptr, entries.n, r, 1, c, 1, &h->buf);

  struct array csr;
  array_init(csr_entry, &csr, 10);

  csr_entry *ptr = entries.ptr;
  sint diag;
  uint e = 0, e1, nr, s;
  while (e < entries.n) {
    s = e;
    nr = 0;
    diag = -1;
    while (e < entries.n && ptr[s].r == ptr[e].r) {
      for (e1 = e; e < entries.n && ptr[e1].c == ptr[e].c; e++)
        ;
      array_cat(csr_entry, &csr, &ptr[e1], 1);
      nr++;
      if (ptr[e1].r == ptr[e1].c)
        diag = csr.n - 1;
    }
    assert(diag >= 0);
    ((csr_entry *)csr.ptr)[diag].v = nr - 1.0;
  }

  array_free(&entries);

  sarray_sort_2(csr_entry, csr.ptr, csr.n, r, 1, c, 1, &h->buf);

  metric_tic(c, CSRMATSETUP);
  csr_mat_setup(h->M, &csr, c, &h->buf);
  array_free(&csr);
  metric_toc(c, CSRMATSETUP);

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

static int genmap_weighted_laplacian_csr_init(genmap_handle h, struct comm *c) {
  metric_tic(c, FINDNBRS);
  struct array entries;
  genmap_find_neighbors(&entries, h, c);
  metric_toc(c, FINDNBRS);

  sarray_sort_2(csr_entry, entries.ptr, entries.n, r, 1, c, 1, &h->buf);

  struct array csr;
  array_init(csr_entry, &csr, 10);

  int nv = genmap_get_nvertices(h);

  csr_entry *ptr = entries.ptr;
  sint diag;
  GenmapScalar v, v1;
  uint e = 0, e1, s;
  while (e < entries.n) {
    s = e;
    v = 0.0;
    diag = -1;
    while (e < entries.n && ptr[s].r == ptr[e].r) {
      v1 = 0.0;
      for (e1 = e;
           e < entries.n && ptr[e1].r == ptr[e].r && ptr[e1].c == ptr[e].c; e++)
        v1 += ptr[e].v;
      ptr[e1].v = v1;
      array_cat(csr_entry, &csr, &ptr[e1], 1);
      v += v1;
      if (ptr[e1].r == ptr[e1].c)
        diag = csr.n - 1;
    }
    assert(diag >= 0);
    ((csr_entry *)csr.ptr)[diag].v = -nv - v;
  }

  array_free(&entries);

  sarray_sort_2(csr_entry, csr.ptr, csr.n, r, 1, c, 1, &h->buf);

  metric_tic(c, CSRMATSETUP);
  csr_mat_setup(h->M, &csr, c, &h->buf);
  array_free(&csr);
  metric_toc(c, CSRMATSETUP);

  return 0;
}

int genmap_laplacian_csr_init(genmap_handle h, struct comm *c) {
  metric_tic(c, LAPLACIANSETUP);

  if (h->M != NULL)
    csr_mat_free(h->M);

  h->M = tmalloc(struct csr_mat_, 1);

  if (h->options->rsb_laplacian_weighted == 0)
    genmap_unweighted_laplacian_csr_init(h, c);
  else
    genmap_weighted_laplacian_csr_init(h, c);

  metric_toc(c, LAPLACIANSETUP);
}

int genmap_laplacian_csr(genmap_handle h, GenmapScalar *u, GenmapScalar *v) {
  csr_mat_apply(v, h->M, u, &h->buf);
  return 0;
}

#undef SWAP
#undef SORT3
