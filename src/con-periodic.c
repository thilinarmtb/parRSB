#include "con-impl.h"

#include <math.h>

//==============================================================================
// Handle periodic BCs
//
int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{1, 5, 7, 3}, {2, 4, 8, 6},
                                                   {1, 2, 6, 5}, {3, 7, 8, 4},
                                                   {1, 3, 4, 2}, {5, 6, 8, 7}};

int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{3, 1, 0, 0}, {2, 4, 0, 0},
                                                   {1, 2, 0, 0}, {4, 3, 0, 0},
                                                   {0, 0, 0, 0}, {0, 0, 0, 0}};

#define distance2D(a, b) (diff_sqr(a.x[0], b.x[0]) + diff_sqr(a.x[1], b.x[1]))
#define distance3D(a, b) (distance2D(a, b) + diff_sqr(a.x[2], b.x[2]))

struct mpair_t {
  uint proc;
  ulong orig, min;
};

static int compress_periodic_vertices(Mesh mesh, struct comm *c, buffer *bfr) {
  parallel_sort(struct point_t, &mesh->elements, globalId, gs_long, 0, 0, c,
                bfr);

  Point points = mesh->elements.ptr;
  uint npoints = mesh->elements.n;

  uint i, nunique = 0;
  if (npoints > 0) {
    ulong current = points[0].globalId;
    points[0].globalId = nunique;
    for (i = 1; i < npoints; i++)
      if (points[i].globalId == current)
        points[i].globalId = nunique;
      else {
        current = points[i].globalId, ++nunique;
        points[i].globalId = nunique;
      }
  }

  slong out[2][1], buf[2][1], in[1];
  if (npoints > 0)
    in[0] = nunique + 1;
  else
    in[0] = 0;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  for (i = 0; i < npoints; i++) points[i].globalId += start;

  return 0;
}

static ulong findMinBelowI(ulong min, uint I, struct array *arr) {
  struct mpair_t *ptr = (struct mpair_t *)arr->ptr;
  for (uint i = 0; i < I; i++)
    if (ptr[i].orig == min) return ptr[i].min;
  return min;
}

static int renumber_periodic_vertices(Mesh mesh, struct comm *c,
                                      struct array *matched, buffer *bfr) {
  uint size1 = mesh->elements.n, size2 = matched->n;
  slong *mids = tcalloc(slong, size1 + 2 * size2),
        *mnew = tcalloc(slong, size1 + 2 * size2),
        *mcur = tcalloc(slong, size1);

  struct point_t *pe = (struct point_t *)mesh->elements.ptr;
  for (uint i = 0; i < size1; i++) mids[i] = pe[i].globalId;
  struct mpair_t *pm = (struct mpair_t *)matched->ptr;
  for (uint i = 0; i < size2; i++) mids[size1 + i] = pm[i].orig;
  struct gs_data *gsh = gs_setup(mids, size1 + size2, c, 0, gs_pairwise, 0);

  for (uint i = 0; i < size1; i++) mnew[i] = pe[i].globalId;
  for (uint i = 0; i < size2; i++) mnew[size1 + i] = pm[i].min;
  gs(mnew, gs_long, gs_min, 0, gsh, bfr);

  sint changed, wrk;
  do {
    for (uint i = 0; i < size1; i++) mcur[i] = mnew[i];
    for (uint i = 0; i < size2; i++) mids[size1 + size2 + i] = -mnew[size1 + i];
    struct gs_data *gsh1 =
        gs_setup(mids, size1 + 2 * size2, c, 0, gs_pairwise, 0);

    gs(mnew, gs_long, gs_min, 0, gsh1, bfr);
    gs_free(gsh1);

    for (uint i = 0; i < size2; i++) mnew[size1 + i] = mnew[size1 + size2 + i];
    gs(mnew, gs_long, gs_min, 0, gsh, bfr);

    changed = 0;
    for (uint i = 0; i < size1; i++) changed += (mnew[i] != mcur[i]);
    comm_allreduce(c, gs_int, gs_max, &changed, 1, &wrk);
  } while (changed);

  for (uint i = 0; i < size1; i++) pe[i].globalId = mcur[i];

  gs_free(gsh);
  free(mids), free(mnew), free(mcur);

  return 0;
}

static int findConnectedPeriodicPairs(Mesh mesh, BoundaryFace f_,
                                      BoundaryFace g_, struct array *matched) {
  struct boundary_t f = *f_, g = *g_;

  int nvf = mesh->nv / 2;
  int ndim = mesh->ndim;

  int i, j;
  scalar fMax = 0.0, gMax = 0.0;

  for (i = 0; i < ndim; i++) {
    scalar meanF = 0.0, meanG = 0.0;

    for (j = 0; j < nvf; j++) {
      fMax = MAX(fMax, fabs(f.face[j].x[i]));
      gMax = MAX(gMax, fabs(g.face[j].x[i]));
      meanF += f.face[j].x[i];
      meanG += g.face[j].x[i];
    }

    for (j = 0; j < nvf; j++) {
      f.face[j].x[i] -= (meanF / nvf);
      g.face[j].x[i] -= (meanG / nvf);
    }
  }

  int shift = 0, k;
  scalar d2Min = 1.e20, d2;
  for (i = 0; i < nvf; i++) {
    d2 = 0.0;
    for (j = 0; j < nvf; j++) {
      k = (j + i) % nvf;
      k = nvf - 1 - k;
      if (ndim == 3)
        d2 += distance3D(f.face[j], g.face[k]);
      else if (ndim == 2)
        d2 += distance2D(f.face[j], g.face[k]);
    }
    if (d2 < d2Min) {
      d2Min = d2;
      shift = i;
    }
  }
  d2Min = sqrt(d2Min);

  scalar fgMax = MAX(fMax, gMax);
  scalar tol = (1e-3) * fgMax;
  if (d2Min > tol) {
    fprintf(stderr,
            "Faces did not match: (d2Min,tol,face1,face2): "
            "%lf %lf %lld %lld\n",
            d2Min, tol, f.faceId, g.faceId);
    exit(1);
  }

  struct mpair_t m;
  for (i = 0; i < nvf; i++) {
    k = (i + shift) % nvf;
    k = nvf - 1 - k;
    m.min = MIN(f.face[i].globalId, g.face[k].globalId);
    m.orig = MAX(f.face[i].globalId, g.face[k].globalId);
    array_cat(struct mpair_t, matched, &m, 1);
  }

  return 0;
}

static int find_connected_periodic_faces(Mesh mesh, struct array *matched) {
  sint bSize = mesh->boundary.n;
  BoundaryFace ptr = mesh->boundary.ptr;
  sint i, j;

  for (i = 0; i < bSize - 1; i++) {
    for (j = i + 1; j < bSize; j++)
      if ((ulong)ptr[j].bc[0] == ptr[i].elementId &&
          (ulong)ptr[j].bc[1] == ptr[i].faceId) {
        findConnectedPeriodicPairs(mesh, &ptr[i], &ptr[j], matched);
      }
  }
  return 0;
}

static int gather_matching_periodic_faces(Mesh mesh, struct comm *c) {
  uint size = c->np;

  BoundaryFace bPtr = mesh->boundary.ptr;
  int nFaces = mesh->boundary.n;

  slong nelgt = mesh->nelgt;
  sint nelt = nelgt / size;
  sint nrem = nelgt - nelt * size;
  slong N = (size - nrem) * nelt;

  sint i;
  slong eid;
  for (i = 0; i < nFaces; i++) {
    eid = MAX((ulong)bPtr[i].bc[0], bPtr[i].elementId);
    if (eid < N)
      bPtr[i].proc = eid / nelt;
    else
      bPtr[i].proc = (eid - N) / (nelt + 1) + size - nrem;
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct boundary_t, &mesh->boundary, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

static int set_periodic_face_coords(Mesh mesh, struct comm *c, buffer *buf) {
  BoundaryFace bPtr = mesh->boundary.ptr;
  sint bSize = mesh->boundary.n;
  if (bSize == 0) return 0;

  Point ePtr = mesh->elements.ptr;
  sint eSize = mesh->elements.n;
  if (eSize == 0) return 0;

  /* Need boundary array to be sorted by elementId */
  sarray_sort(struct boundary_t, bPtr, bSize, elementId, 1, buf);

  /* Need element array to be sorted by sequenceId */
  sarray_sort(struct point_t, ePtr, eSize, sequenceId, 1, buf);

  int faces[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
  if (mesh->ndim == 3)
    memcpy(faces, faces3D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));
  else
    memcpy(faces, faces2D, GC_MAX_FACES * GC_MAX_FACE_VERTICES * sizeof(int));

  sint i = 0, k = 0;
  int nv = mesh->nv, nvf = mesh->nv / 2, j;
  while (i < bSize) {
    while (k < eSize && ePtr[k].elementId < bPtr[i].elementId) k += nv;
    // copy vertices to boundary face
    if (k < eSize && ePtr[k].elementId == bPtr[i].elementId) {
      int faceId = bPtr[i].faceId;
      for (j = 0; j < nvf; j++)
        bPtr[i].face[j] = ePtr[k + faces[faceId][j] - 1];
    }
    i++;
  }

  return 0;
}

int match_periodic_faces(Mesh mesh, struct comm *c, int verbose, buffer *bfr) {
  const char *functions[6] = {
      "set_periodic_face_coords      ", "gather_matching_periodic_faces",
      "find_connected_periodic_faces ", "renumber_periodic_vertices    ",
      "compress_periodic_vertices    ", "send_back                     "};

  parrsb_print(c, verbose, "\t\t%s ...", functions[0]);
  set_periodic_face_coords(mesh, c, bfr);

  parrsb_print(c, verbose, "\t\t%s ...", functions[1]);
  gather_matching_periodic_faces(mesh, c);

  struct array matched;
  array_init(struct mpair_t, &matched, 10);
  matched.n = 0;

  parrsb_print(c, verbose, "\t\t%s ...", functions[2]);
  find_connected_periodic_faces(mesh, &matched);

  parrsb_print(c, verbose, "\t\t%s ...", functions[3]);
  renumber_periodic_vertices(mesh, c, &matched, bfr);
  array_free(&matched);

  parrsb_print(c, verbose, "\t\t%s ...", functions[4]);
  compress_periodic_vertices(mesh, c, bfr);

  parrsb_print(c, verbose, "\t\t%s ...", functions[5]);
  send_back(mesh, c, bfr);

  return 0;
}

static scalar transform_face(scalar R[4][3], const scalar face[4][3],
                             const scalar *const coord, sint nv, sint ndim) {
  return DBL_MAX;
}

static sint calculate_R_and_t(scalar R[4][3], long long *const gid, uint nf,
                              const sint *const bid, sint nv, sint ndim,
                              const scalar *const coord, struct comm *c) {
  sint root = INT_MAX;
  uint index = UINT_MAX;
  for (uint i = 0; i < nf; i++) {
    if (bid[i] == 1) continue;
    root = c->id, index = i;
    break;
  }

  char wrk[MAX(sizeof(scalar), sizeof(sint))];
  comm_allreduce(c, gs_int, gs_min, &root, 1, &wrk);

  scalar face[4][3];
  if (c->id != root) goto bcast_face;
  for (uint i = 0; i < nv * ndim; i++)
    face[0][i] = coord[index * nv * ndim + i];
bcast_face:
  comm_bcast(c, face, nv * ndim * sizeof(scalar), root);

  scalar error = DBL_MAX;
  index = UINT_MAX;
  for (uint i = 0; i < nf; i++) {
    if (bid[i] == 0) continue;
    scalar error_i = transform_face(R, face, &coord[i * nv * ndim], nv, ndim);
    if (error_i < error) {
      error = error_i;
      index = i;
    }
  }

  scalar global_error = error;
  comm_allreduce(c, gs_scalar, gs_min, &global_error, 1, &wrk);

  if (global_error > 1e-14) return 1;

  if (fabs(global_error - error) < 1e-13)
    root = c->id;
  else
    root = INT_MAX;
  comm_allreduce(c, gs_int, gs_min, &root, 1, &wrk);

  if (c->id != root) goto bcast_R_and_t;
  transform_face(R, face, &coord[index * nv * ndim], nv, ndim);
bcast_R_and_t:
  comm_bcast(c, R, 4 * 3 * sizeof(scalar), root);

  return 0;
}

static void transform_points(scalar *coord_, sint nf, const sint *const bid,
                             sint nv, sint ndim, const scalar *const coord,
                             const scalar R[4][3]) {
  const size_t nc = (size_t)nf * nv;
  const scalar t[3] = {R[3][0], R[3][1], R[3][2]};

  if (ndim == 2) goto dim2;

  for (uint i = 0; i < nc; i++) {
    const scalar *coord_i = &coord[i * ndim];
    scalar *coord_o = &coord_[i * ndim];

    if (bid[i] == 0) {
      coord_o[0] = coord_i[0];
      coord_o[1] = coord_i[1];
      coord_o[2] = coord_i[2];
    } else if (bid[i] == 1) {
      coord_o[0] = R[0][0] * coord_i[0] + R[0][1] * coord_i[1] +
                   R[0][2] * coord_i[2] + t[0];
      coord_o[1] = R[1][0] * coord_i[0] + R[1][1] * coord_i[1] +
                   R[1][2] * coord_i[2] + t[1];
      coord_o[2] = R[2][0] * coord_i[0] + R[2][1] * coord_i[1] +
                   R[2][2] * coord_i[2] + t[2];
    }
  }

  return;
dim2:
  for (uint i = 0; i < nc; i++) {
    const scalar *coord_i = &coord[i * ndim];
    scalar *coord_o = &coord_[i * ndim];

    if (bid[i] == 0) {
      coord_o[0] = coord_i[0];
      coord_o[1] = coord_i[1];
    } else if (bid[i] == 1) {
      coord_o[0] = R[0][0] * coord_i[0] + R[0][1] * coord_i[1];
      coord_o[1] = R[1][0] * coord_i[0] + R[1][1] * coord_i[1];
    }
  }
}

static int number_points(long long *const gid, sint nf, sint nv, sint ndim,
                         scalar *coord, scalar tol, struct comm *c,
                         buffer *bfr) {
  unsigned nnbrs = (nv == 4) ? 2 : 1;
  struct mesh_t *mesh = mesh_init(nf, nv, ndim, nnbrs, coord, 0, 0, c);

  find_min_neighbor_distance(mesh);

  find_unique_vertices(mesh, c, tol, 0, bfr);
  set_global_id(mesh, c);
  send_back(mesh, c, bfr);

  con_chk_err(element_check(mesh, c, bfr), "element check failed.", c);
  con_chk_err(face_check(mesh, c, bfr), "face_check failed.", c);

  mesh_free(mesh);
  return 0;
}

static void update_global_ids(long long *const gid, sint nf, sint nv,
                              const long long *const new_gid,
                              const struct comm *const c, buffer *bfr) {
  const size_t size = (size_t)nf * nv;
  struct gs_data *gsh = gs_setup(new_gid, size, c, 0, gs_pairwise, 0);

  gs(gid, gs_long, gs_min, 0, gsh, bfr);

  gs_free(gsh);
}

int match_periodic_faces_automatically(long long *const gid, uint nf,
                                       const sint *const bid, sint nv,
                                       sint ndim, const scalar *const coord,
                                       scalar tol, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  sint err = 0;

  // Match one periodic face from bid = 0 with all the faces of bid = 1 and
  // calculate the rotation matrix and translation vector in case there is
  // a match.
  scalar R[4][3];
  err |= calculate_R_and_t(R, gid, nf, bid, nv, ndim, coord, &c);

  // Transform the points in B using R and t.
  scalar *transformed_coord = tcalloc(scalar, nf * nv * ndim);
  transform_points(transformed_coord, nf, bid, nv, ndim, coord, R);

  // Globally number points:
  buffer bfr;
  buffer_init(&bfr, 1024);

  long long *new_gid = tcalloc(long long, nf *nv);
  err |= number_points(new_gid, nf, nv, ndim, transformed_coord, tol, &c, &bfr);

  update_global_ids(gid, nf, nv, new_gid, &c, &bfr);

  sint wrk;
  comm_allreduce(&c, gs_int, gs_max, &err, 1, &wrk);

  free(transformed_coord), free(new_gid);
  buffer_free(&bfr), comm_free(&c);

  return err;
}

#undef distance2D
#undef distance3D
