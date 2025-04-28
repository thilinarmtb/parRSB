#include "con-impl.h"

int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};
int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS] = {
    {1, 2, 4}, {0, 3, 5}, {0, 3, 6}, {1, 2, 7},
    {0, 5, 6}, {1, 4, 7}, {2, 4, 7}, {3, 5, 6}};

/*
 * Mesh struct.
 */
struct mesh_t *mesh_init(uint nelt, unsigned nv, unsigned ndim, unsigned nnbrs,
                         double *coord, long long *pinfo, uint npinfo,
                         const struct comm *c) {
  struct mesh_t *m = tcalloc(struct mesh_t, 1);
  m->nelt = nelt;
  m->nv = nv;
  m->ndim = ndim;
  m->nnbrs = nnbrs;

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  ulong start = out[0][0];
  m->nelgt = out[1][0];

  array_init(struct point_t, &m->elements, nelt * nv);
  struct point_t p = {.origin = c->id};
  for (uint i = 0; i < nelt; i++) {
    for (uint k = 0; k < nv; k++) {
      uint j = PRE_TO_SYM_VERTEX[k];
      for (uint l = 0; l < ndim; l++)
        p.x[l] = coord[i * nv * ndim + j * ndim + l];
      p.elementId = start + i, p.sequenceId = nv * (start + i) + k;
      array_cat(struct point_t, &m->elements, &p, 1);
    }
  }

  array_init(struct boundary_t, &m->boundary, npinfo);
  struct boundary_t b;
  for (uint i = 0; i < npinfo; i++) {
    b.elementId = pinfo[4 * i + 0] - 1;
    b.faceId = PRE_TO_SYM_FACE[pinfo[4 * i + 1] - 1];
    b.bc[0] = pinfo[4 * i + 2] - 1;
    b.bc[1] = PRE_TO_SYM_FACE[pinfo[4 * i + 3] - 1];
    array_cat(struct boundary_t, &m->boundary, &b, 1);
  }

  return m;
}

int mesh_free(struct mesh_t *m) {
  array_free(&m->elements), array_free(&m->boundary), free(m), m = 0;
  return 0;
}

/*
 * Find the minimum distance between a vertex and its neighbors.
 */
static inline double distance_2d(struct point_t *a, struct point_t *b) {
  return diff_sqr(a->x[0], b->x[0]) + diff_sqr(a->x[1], b->x[1]);
}

static inline double distance_3d(struct point_t *a, struct point_t *b) {
  return distance_2d(a, b) + diff_sqr(a->x[2], b->x[2]);
}

int find_min_neighbor_distance(Mesh mesh) {
  struct point_t *p = (struct point_t *)mesh->elements.ptr;
  uint ndim = mesh->ndim;
  uint nv = mesh->nv;

  if (ndim < 2 || ndim > 3) return 1;

  uint i, j, k, neighbor;
  if (ndim == 3) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          scalar d = distance_3d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  }

  if (ndim == 2) {
    for (i = 0; i < mesh->elements.n; i += nv) {
      for (j = 0; j < nv; j++) {
        p[i + j].dx = SCALAR_MAX;
        for (k = 0; k < mesh->nnbrs; k++) {
          neighbor = NEIGHBOR_MAP[j][k];
          scalar d = distance_2d(&p[i + j], &p[i + neighbor]);
          p[i + j].dx = MIN(p[i + j].dx, d);
        }
      }
    }
  }

  return 0;
}

/*
 * Global numbering.
 */
int set_global_id(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = (struct point_t *)mesh->elements.ptr;

  sint bin = (nPoints > 0);
  struct comm nonZeroRanks;
  comm_split(c, bin, c->id, &nonZeroRanks);

  if (bin == 1) {
    slong count = 0;
    for (uint i = 0; i < nPoints; i++)
      if (points[i].ifSegment) count++;

    slong in = count, out[2][1], buf[2][1];
    comm_scan(out, &nonZeroRanks, gs_long, gs_add, &in, 1, buf);
    slong start = out[0][0];

    count = -1;
    for (uint i = 0; i < nPoints; i++) {
      if (points[i].ifSegment) count++;
      assert(start + count >= 0);
      points[i].globalId = start + count;
    }
  }

  comm_free(&nonZeroRanks);

  return 0;
}

int send_back(Mesh mesh, struct comm *c, buffer *bfr) {
  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct point_t, &mesh->elements, origin, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct point_t, mesh->elements.ptr, mesh->elements.n, sequenceId,
              1, bfr);

  return 0;
}

static int transfer_boundary_faces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  struct array *boundary = &mesh->boundary;
  BoundaryFace ptr = (struct boundary_t *)boundary->ptr;
  int nFaces = boundary->n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = ptr[i].elementId;
    if (eid < N)
      ptr[i].proc = eid / nelt;
    else
      ptr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct boundary_t, boundary, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

/*
 * Sanity checks.
 */
typedef struct {
  ulong sequenceId;
  int nnbrs;
  ulong elementId;
  ulong vertexId;
  uint workProc;
} vertex;

typedef struct {
  uint *offsets;
  ulong *elements;
  ulong *globalIds;
  uint size;
} VToEMap;

typedef struct {
  ulong id;
} LongID;

typedef struct {
  ulong procId;
} ProcID;

static VToEMap *getVToEMap(Mesh m, struct comm *c, buffer *bfr) {
  uint nelt = m->nelt;
  uint nv = m->nv;

  slong out[2][1], buf[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  ulong elemId = out[0][0];
  ulong sequenceId = elemId * nv;

  size_t size = nelt * (size_t)nv;
  struct array vertices;
  array_init(vertex, &vertices, size);

  // Create (globalId, elementId) pairs and send them to globalId % np
  Point ptr = m->elements.ptr;
  uint i, j;
  for (i = 0; i < nelt; i++) {
    for (j = 0; j < nv; j++) {
      ulong globalId = ptr[i * nv + j].globalId + 1;
      vertex t = {.elementId = elemId,
                  .sequenceId = sequenceId,
                  .vertexId = globalId,
                  .workProc = globalId % c->np};
      array_cat(vertex, &vertices, &t, 1);
      sequenceId++;
    }

    elemId++;
  }

  sarray_sort_2(vertex, vertices.ptr, vertices.n, vertexId, 1, elementId, 1,
                bfr);

  struct array vtcsCmpct;
  array_init(vertex, &vtcsCmpct, 10);
  vertex *vPtr = (vertex *)vertices.ptr;

  if (vertices.n > 0) {
    vertex prev = vPtr[0];
    array_cat(vertex, &vtcsCmpct, &prev, 1);

    for (i = 1; i < vertices.n; i++) {
      if ((vPtr[i].elementId != prev.elementId) ||
          (vPtr[i].vertexId != prev.vertexId)) {
        prev = vPtr[i];
        array_cat(vertex, &vtcsCmpct, &prev, 1);
      }
    }
  }
  array_free(&vertices);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(vertex, &vtcsCmpct, workProc, 1, &cr);

  // Find all the elements which share globalId and send the union
  // back to all the processors which has globalId
  // FIXME: Assumes quads or hexes
  vPtr = vtcsCmpct.ptr;
  sarray_sort_2(vertex, vPtr, vtcsCmpct.n, vertexId, 1, workProc, 0, bfr);

  struct array a;
  array_init(vertex, &a, 10);
  struct array procs;
  array_init(ProcID, &procs, 10);

  vPtr = vtcsCmpct.ptr;
  uint s = 0, e;
  vertex t;
  ProcID p;
  while (s < vtcsCmpct.n) {
    procs.n = 0;

    p.procId = vPtr[s].workProc;
    array_cat(ProcID, &procs, &p, 1);
    for (e = s + 1; e < vtcsCmpct.n && vPtr[s].vertexId == vPtr[e].vertexId;
         e++) {
      if (vPtr[e].workProc != p.procId) {
        p.procId = vPtr[e].workProc;
        array_cat(ProcID, &procs, &p, 1);
      }
    }

    ProcID *pPtr = procs.ptr;
    e = MIN(e, vtcsCmpct.n);
    for (i = 0; i < procs.n; i++) {
      t.workProc = pPtr[i].procId;
      for (j = s; j < e; j++) {
        t.vertexId = vPtr[j].vertexId;
        t.sequenceId = vPtr[j].sequenceId;
        t.elementId = vPtr[j].elementId;
        array_cat(vertex, &a, &t, 1);
      }
    }
    s = e;
  }
  array_free(&vtcsCmpct);
  array_free(&procs);

  sarray_transfer(vertex, &a, workProc, 1, &cr);
  sarray_sort_2(vertex, a.ptr, a.n, vertexId, 1, elementId, 1, bfr);
  crystal_free(&cr);

  // create the map
  if (a.n == 0) return NULL;

  VToEMap *map = calloc(1, sizeof(VToEMap));
  map->elements = calloc(a.n, sizeof(ulong));

  uint nGIds = 1, prev = 0;
  vertex *aPtr = (vertex *)a.ptr;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId) nGIds++;
    prev = i;
  }

  map->size = nGIds;
  map->globalIds = calloc(nGIds, sizeof(ulong));
  map->offsets = calloc(nGIds + 1, sizeof(ulong));

  map->elements[0] = aPtr[0].elementId;
  map->globalIds[0] = aPtr[0].vertexId;
  map->offsets[0] = 0;

  prev = 0;
  uint nOffsets = 0;
  for (i = 1; i < a.n; i++) {
    if (aPtr[i].vertexId != aPtr[prev].vertexId) {
      nOffsets++;
      map->globalIds[nOffsets] = aPtr[i].vertexId;
      map->offsets[nOffsets] = prev = i;
    }
    map->elements[i] = aPtr[i].elementId;
  }
  map->offsets[++nOffsets] = a.n;
  assert(nOffsets == nGIds);

  array_free(&a);

  return map;
}

// key must be present in globalIds
static uint getPosition(VToEMap *map, ulong key) {
  ulong *globalIds = map->globalIds;

  uint begin = 0;
  uint end = map->size;
  uint mid = 0;
  while (begin < end) {
    mid = (begin + end) / 2;

    if (key == globalIds[mid])
      return mid;
    else if (key < globalIds[mid])
      end = mid;
    else
      begin = mid;
  };

  if (globalIds[mid] != key) return UINT_MAX;
  return mid;
}

static void freeVToEMap(VToEMap *map) {
  free(map->globalIds);
  free(map->offsets);
  free(map->elements);
  free(map);
}

static int face_check(Mesh mesh, struct comm *c, buffer *bfr) {
  VToEMap *map = getVToEMap(mesh, c, bfr);

  uint nelt = mesh->nelt;
  uint ndim = mesh->ndim;

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (ndim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  Point ptr = mesh->elements.ptr;
  uint nf = (ndim == 3) ? 6 : 4;
  uint nfv = (ndim == 3) ? 4 : 2;
  uint nv = (ndim == 3) ? 8 : 4;

  struct array shared;
  array_init(LongID, &shared, 200);

  int err = 0;

  uint i, j, k, l;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nf && err == 0; j++) {
      shared.n = 0;

      for (k = 0; k < nfv; k++) {
        ulong globalId = ptr[i * nv + faces[j][k] - 1].globalId + 1;
        uint indx = getPosition(map, globalId);
        assert(indx < UINT_MAX);
        LongID elemId;
        for (l = map->offsets[indx]; l < map->offsets[indx + 1]; l++) {
          elemId.id = map->elements[l];
          array_cat(LongID, &shared, &elemId, 1);
        }
      }

      sarray_sort(LongID, shared.ptr, shared.n, id, 1, bfr);

      ulong prev = 0;
      int ncount = 1;
      LongID *sptr = shared.ptr;
      for (l = 1; l < shared.n; l++) {
        if (sptr[l].id != sptr[prev].id) {
          if (ncount == 3) {
            err = 1;
            break;
          }
          prev = l;
          ncount = 1;
        } else {
          ncount++;
        }
      }

      if (ncount == 3) {
        err = 1;
        break;
      }
    }
  }

  array_free(&shared);
  freeVToEMap(map);

  return err;
}

static int element_check(Mesh mesh, struct comm *c, buffer *bfr) {
  uint nelt = mesh->nelt;
  uint ndim = mesh->ndim;
  uint nv = (ndim == 3) ? 8 : 4;

  LongID globalIds[8];
  Point ptr = mesh->elements.ptr;
  uint i, j;
  int err = 0;
  for (i = 0; i < nelt && err == 0; i++) {
    for (j = 0; j < nv; j++) globalIds[j].id = ptr[i * nv + j].globalId + 1;

    sarray_sort(LongID, globalIds, nv, id, 1, bfr);

    for (j = 0; j < nv - 1; j++) {
      if (globalIds[j].id == globalIds[j + 1].id) {
        err = 1;
        break;
      }
    }
  }

  return err;
}

/*
 * Input:
 *   nelt: Number of elements
 *   ndim: The dimension of the problem
 *   nv: Number of vertices in an element
 *     nv = 8 if ndim == 3 (Hex)
 *     nv = 4 if ndim = 2 (Quad)
 *   pinfo: Periodic BC information
 *   npinfo: Number of periodic BCs
 *   tol: Tolerance similar to gencon
 *   coord [nelt, nv, ndim]: Coordinates of elements vertices in preprocessor
 *     ordering.
 * Output:
 *   vtx[nelt, nv]: Global numbering of vertices of elements
 */
int parrsb_conn_mesh(long long *vtx, double *coord, uint nelt, unsigned ndim,
                     long long *pinfo, int npinfo, double tol, MPI_Comm comm,
                     int verbose) {
  struct comm c;
  comm_init(&c, comm);

  parrsb_print(&c, verbose, "Running parCon ...");

  double duration[8] = {0};
  const char *name[8] = {
      "transfer_boundary_faces    ", "find_min_neighbor_distance ",
      "find_unique_vertices       ", "set_global_id              ",
      "element_check              ", "face_check                 ",
      "match_periodic_faces       ", "copy_output                "};

  buffer bfr;
  buffer_init(&bfr, 1024);

  parrsb_barrier(&c);
  double tall = comm_time(), t;

  unsigned nv = (ndim == 3) ? 8 : 4;
  unsigned nnbrs = ndim;
  Mesh mesh = mesh_init(nelt, nv, ndim, nnbrs, coord, pinfo, npinfo, &c);

  parrsb_print(&c, verbose - 1, "\t%s ...", name[0]);
  parrsb_barrier(&c), t = comm_time();
  transfer_boundary_faces(mesh, &c);
  duration[0] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[1]);
  parrsb_barrier(&c), t = comm_time();
  find_min_neighbor_distance(mesh);
  duration[1] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[2]);
  parrsb_barrier(&c), t = comm_time();
  find_unique_vertices(mesh, &c, tol, verbose - 1, &bfr);
  duration[2] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[3]);
  parrsb_barrier(&c), t = comm_time();
  set_global_id(mesh, &c);
  send_back(mesh, &c, &bfr);
  duration[3] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[4]);
  parrsb_barrier(&c), t = comm_time();
  sint err = 0;
  con_chk_err(element_check(mesh, &c, &bfr), err, &c);
  if (err) {
    if (c.id == 0) fprintf(stderr, "parCon: element check failed !\n");
    goto cleanup;
  }
  duration[4] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[5]);
  parrsb_barrier(&c), t = comm_time();
  con_chk_err(face_check(mesh, &c, &bfr), err, &c);
  if (err) {
    if (c.id == 0) fprintf(stderr, "parCon: face check failed !\n");
    goto cleanup;
  }
  duration[5] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[6]);
  parrsb_barrier(&c), t = comm_time();
  match_periodic_faces(mesh, &c, verbose - 1, &bfr);
  duration[6] = comm_time() - t;

  parrsb_print(&c, verbose - 1, "\t%s ...", name[7]);
  parrsb_barrier(&c), t = comm_time();
  Point ptr = mesh->elements.ptr;
  for (uint i = 0; i < nelt; i++) {
    for (uint j = 0; j < mesh->nv; j++)
      vtx[i * mesh->nv + j] = ptr[i * mesh->nv + j].globalId + 1;
  }
  duration[7] = comm_time() - t;

  // Report timing info and finish
  double gmin[8], gmax[8], buf[8];
  for (unsigned i = 0; i < 8; i++) gmax[i] = gmin[i] = duration[i];
  comm_allreduce(&c, gs_double, gs_min, gmin, 8, buf);
  comm_allreduce(&c, gs_double, gs_max, gmax, 8, buf);

  for (unsigned i = 0; i < 7; i++) {
    parrsb_print(&c, verbose - 1, "%s: %e %e (min max)", name[i], gmin[i],
                 gmax[i]);
  }

  parrsb_barrier(&c);
  tall = comm_time() - tall;
  parrsb_print(&c, verbose, "parCon (tol = %e) finished in %g s", tol, tall);

cleanup:
  buffer_free(&bfr), mesh_free(mesh), comm_free(&c);

  return err;
}

/*
 * Fortran interface.
 */
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nelt, int *ndim,
                       long long *pinfo, int *npinfo, double *tol,
                       MPI_Fint *fcomm, int *err, int *verbose) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*fcomm);
  *err = parrsb_conn_mesh(vtx, coord, *nelt, *ndim, pinfo, *npinfo, *tol, c,
                          *verbose);
}
