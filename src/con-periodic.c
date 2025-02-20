#include "con-impl.h"

#include <math.h>
#include <time.h>

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

/*
 * Lanczos algorithm for symmetric eigenvalue problems.
 */
static inline scalar dot(const scalar *a, const scalar *b, sint n) {
  scalar sum = 0.0;
  for (sint i = 0; i < n; i++) sum += a[i] * b[i];
  return sum;
}

static inline void scale(scalar *a, const scalar *b, scalar c, sint n) {
  for (sint i = 0; i < n; i++) a[i] = b[i] * c;
}

static inline void add2s1(scalar *a, scalar *b, scalar c, sint n) {
  for (sint i = 0; i < n; i++) a[i] = c * a[i] + b[i];
}

static inline void add2s2(scalar *a, scalar *b, scalar c, sint n) {
  for (sint i = 0; i < n; i++) a[i] = a[i] + c * b[i];
}

static inline void ax(scalar *y, const scalar A[3][3], const scalar *x,
                      sint n) {
  for (sint i = 0; i < n; i++) {
    scalar sum = 0.0;
    for (sint j = 0; j < n; j++) sum += A[i][j] * x[j];
    y[i] = sum;
  }
}

static sint lanczos(scalar *diag, scalar *upper, scalar rr[3][4], sint ndim,
                    const scalar A[3][3], scalar tol) {
  scalar r[3], p[3], w[3];
  for (sint i = 0; i < ndim; i++) {
    r[i] = rand() / (scalar)RAND_MAX;
    p[i] = w[i] = 0.0;
  }

  scalar rtr = dot(r, r, ndim);
  scalar rnorm = sqrt(rtr);
  scalar rtol = rnorm * tol;
  scalar rni = 1.0 / rnorm;

  scalar rr_[4][3];
  scale(&rr_[0][0], r, rni, ndim);

  scalar rtz1 = 1, rtz2;
  scalar pap = 0, pap_old;
  scalar alpha, beta;

  sint iter;
  for (iter = 0; iter < ndim; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0) beta = 0.0;

    add2s1(p, r, beta, ndim);
    ax(w, A, p, ndim);

    pap_old = pap, pap = dot(w, p, ndim);

    alpha = rtz1 / pap;
    add2s2(r, w, -alpha, ndim);

    rtr = dot(r, r, ndim);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;
    scale(&rr_[iter + 1][0], r, rni, ndim);

    if (iter == 0) {
      diag[iter] = pap / rtz1;
    } else {
      diag[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      iter++;
      break;
    }
  }

  // Transpose rr
  for (int i = 0; i < ndim; i++)
    for (int j = 0; j < iter; j++) rr[i][j] = rr_[j][i];

  return iter;
}

/*
 * Power iteration.
 */
static scalar norm2(const scalar *x, sint n) { return sqrt(dot(x, x, n)); }

static scalar normi(const scalar *x, sint n) {
  scalar norm = 0;
  for (sint i = 0; i < n; i++)
    if (norm < fabs(x[i])) norm = fabs(x[i]);
  return norm;
}

static void ortho(scalar *x, const scalar *y, sint n) {
  scalar dotxy = dot(x, y, n);
  scalar normy = sqrt(dot(y, y, n));
  scalar y_[3];
  scale(y_, y, 1.0 / normy, n);
  add2s2(x, y_, -dotxy / normy, n);
}

static void copy(scalar *a, const scalar *b, sint n) {
  for (sint i = 0; i < n; i++) a[i] = b[i];
}

static int power_iteration(scalar *evecs, scalar *evals, sint ndim,
                           const scalar A[3][3], scalar tol) {
  for (int i = 0; i < ndim; i++) {
    scalar x[3];
    for (sint i = 0; i < 3; i++) x[i] = rand() / (scalar)RAND_MAX;

    scalar norm_1 = norm2(x, ndim);
    scale(x, x, 1.0 / norm_1, ndim);
    scalar norm_0;
    do {
      ax(evecs + i * ndim, A, x, ndim);
      for (sint j = 0; j < i; j++)
        ortho(evecs + i * ndim, evecs + j * ndim, ndim);

      norm_0 = norm_1;
      norm_1 = norm2(evecs + i * ndim, ndim);
      scale(x, evecs + i * ndim, 1.0 / norm_1, ndim);
    } while (fabs(norm_1 - norm_0) > tol);

    evals[i] = norm_1;
    copy(evecs + i * ndim, x, ndim);
  }

  return 0;
}

static void svd(scalar U[3][3], scalar S[3][3], scalar V[3][3], sint ndim,
                scalar A[3][3], scalar tol) {
  // Find the eigenvectors of A^T A. These are the columns of V.
  scalar ATA[3][3];
  for (int i = 0; i < 9; i++) ATA[0][i] = 0.0;

  for (sint i = 0; i < ndim; i++) {
    for (sint j = 0; j < ndim; j++)
      for (sint k = 0; k < ndim; k++) ATA[i][j] += A[k][i] * A[k][j];
  }

  scalar evecs[9], evals[3];
  power_iteration(evecs, evals, ndim, ATA, tol);

  // V = evecs
  for (sint i = 0; i < 3; i++)
    for (sint j = 0; j < ndim; j++) V[j][i] = evecs[i * ndim + j];

  // Singular values are the square roots of the eigenvalues of A^T A.
  for (int i = 0; i < 9; i++) S[0][i] = 0.0;
  for (int i = 0; i < 3; i++) S[i][i] = evals[i] = sqrt(evals[i]);

  // Calculate U.
  for (sint i = 0; i < ndim; i++) {
    for (sint j = 0; j < ndim; j++) {
      U[j][i] = 0.0;
      for (sint k = 0; k < ndim; k++) U[j][i] += A[j][k] * V[k][i];
      U[j][i] /= evals[i];
    }
  }
}

// find the translation vector `t` and the rotation matrix `R` that maps the
// face `face0` to the face `face1`.
static scalar transform_face(scalar R[3][3], scalar t[3],
                             const scalar face1[4][3], const scalar face0[4][3],
                             sint nv, sint ndim, scalar tol) {
  // Find the centroids of the faces.
  scalar t0[3], t1[3];
  for (int i = 0; i < 3; i++) t0[i] = t1[i] = 0.0;
  for (sint i = 0; i < nv; i++)
    for (sint j = 0; j < ndim; j++) t1[j] += face1[i][j], t0[j] += face0[i][j];
  for (sint i = 0; i < ndim; i++) t1[i] /= nv, t0[i] /= nv;

  // Next we find the rotation matrix R. To do so, we form the (face0^T x
  // face1) matrix (ndim x ndim).
  scalar C[3][3];
  for (sint i = 0; i < ndim; i++) {
    for (sint j = 0; j < ndim; j++) {
      C[i][j] = 0.0;
      for (sint k = 0; k < nv; k++)
        C[i][j] += (face0[k][i] - t0[i]) * (face1[k][j] - t1[j]);
    }
  }

  // Compute the SVD of the matrix C = USV^T.
  scalar U[3][3], S[3][3], V[3][3];
  svd(U, S, V, ndim, C, tol);

  // Calculate the rotation matrix R = VU^T.
  for (sint i = 0; i < ndim; i++) {
    for (sint j = 0; j < ndim; j++) {
      R[i][j] = 0.0;
      for (sint k = 0; k < ndim; k++) R[i][j] += V[i][k] * U[j][k];
    }
  }

  // Calculate the translation vector t (fix this).
  for (sint i = 0; i < ndim; i++) {
    t[i] = t1[i];
    for (sint j = 0; j < ndim; j++) t[i] -= R[i][j] * t0[j];
  }

  // Transform the face0 using R and t.
  scalar face1_[4][3];
  for (sint i = 0; i < nv; i++) {
    for (sint j = 0; j < ndim; j++) {
      face1_[i][j] = 0.0;
      for (sint k = 0; k < ndim; k++) face1_[i][j] += R[j][k] * face0[i][k];
      face1_[i][j] += t[j];
    }
  }

  // Calculate the error.
  scalar err = 0;
  for (sint i = 0; i < nv; i++)
    for (sint j = 0; j < ndim; j++) err += diff_sqr(face1[i][j], face1_[i][j]);

  return sqrt(err / (nv * ndim));
}

static sint calculate_R_and_t(scalar R[3][3], scalar t[3], slong *const gid,
                              uint nf, const sint *const bid, sint nv,
                              sint ndim, const scalar *const coord, scalar tol,
                              struct comm *c) {
  srand(time(0));

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
  for (sint i = 0; i < nv; i++)
    for (sint j = 0; j < ndim; j++) face[i][j] = coord[(index + i) * ndim + j];
bcast_face:
  comm_bcast(c, face, nv * ndim * sizeof(scalar), root);

  scalar error = DBL_MAX;
  index = UINT_MAX;
  for (uint i = 0; i < nf; i++) {
    if (bid[i] == 0) continue;
    scalar error_i = transform_face(
        R, t, face, (const scalar(*)[3])(&coord[i * nv * ndim]), nv, ndim, tol);
    if (error_i < error) {
      error = error_i;
      index = i;
    }
  }

  scalar global_error = error;
  comm_allreduce(c, gs_scalar, gs_min, &global_error, 1, &wrk);

  if (global_error > tol) return 1;

  const scalar eps = 1e-15;
  if (fabs(global_error - error) < eps)
    root = c->id;
  else
    root = INT_MAX;
  comm_allreduce(c, gs_int, gs_min, &root, 1, &wrk);

  if (c->id != root) goto bcast_R_and_t;
  transform_face(R, t, face, (const scalar(*)[3])(&coord[index * nv * ndim]),
                 nv, ndim, tol);
bcast_R_and_t:
  comm_bcast(c, R, 3 * 3 * sizeof(scalar), root);
  comm_bcast(c, t, 3 * sizeof(scalar), root);

  return 0;
}

// `coordo` should be zero initialized.
static void transform_points(scalar *coordo, sint nf, const sint *const bid,
                             sint nv, sint ndim, const scalar *const coordi,
                             const scalar R[4][3]) {
  const scalar t[3] = {R[3][0], R[3][1], R[3][2]};

  for (size_t i = 0; i < nf; i++) {
    const size_t offset = i * (size_t)nv * (size_t)ndim;
    const scalar *coordi_i = &coordi[offset];
    scalar *coordo_i = &coordo[offset];

    if (bid[i] == 0) {
      for (int j = 0; j < nv; j++) {
        for (int k = 0; k < ndim; k++)
          coordo_i[j * nv + k] = coordi_i[j * nv + k];
      }
    } else if (bid[i] == 1) {
      for (int j = 0; j < nv; j++) {
        for (int k = 0; k < ndim; k++) {
          for (int l = 0; l < ndim; l++)
            coordo_i[j * nv + l] += R[k][l] * coordi_i[j * nv + l] + t[l];
        }
      }
    }
  }
}

static int number_points(slong *const gid, sint nf, sint nv, sint ndim,
                         scalar *coord, scalar tol, struct comm *c,
                         buffer *bfr) {
  unsigned nnbrs = (nv == 4) ? 2 : 1;
  struct mesh_t *mesh = mesh_init(nf, nv, ndim, nnbrs, coord, 0, 0, c);

  find_min_neighbor_distance(mesh);

  find_unique_vertices(mesh, c, tol, 0, bfr);
  set_global_id(mesh, c);
  send_back(mesh, c, bfr);

#define cleanup_before_return()                                                \
  { mesh_free(mesh); }

  con_chk_err(element_check(mesh, c, bfr), "element check failed.", c);
  con_chk_err(face_check(mesh, c, bfr), "face_check failed.", c);

  cleanup_before_return();
#undef cleanup_before_return
  return 0;
}

static void update_global_ids(slong *const gid, sint nf, sint nv,
                              const slong *const new_gid,
                              const struct comm *const c, buffer *bfr) {
  const size_t size = (size_t)nf * nv;
  struct gs_data *gsh = gs_setup(new_gid, size, c, 0, gs_pairwise, 0);

  gs(gid, gs_long, gs_min, 0, gsh, bfr);

  gs_free(gsh);
}

int automatic_periodic_face_match(slong *const gid, uint nf,
                                  const sint *const bid, sint nv, sint ndim,
                                  const scalar *const coord, scalar tol,
                                  MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  const size_t ngids = (size_t)nf * (size_t)nv;
  slong *new_gid = tcalloc(slong, ngids);

  const size_t ncoords = ngids * (size_t)ndim;
  scalar *transformed_coord = tcalloc(scalar, ncoords);

  buffer bfr;
  buffer_init(&bfr, 1024);

#define cleanup_before_return()                                                \
  {                                                                            \
    free(transformed_coord), free(new_gid);                                    \
    buffer_free(&bfr), comm_free(&c);                                          \
  }

  // Match one periodic face from bid = 0 with all the faces of bid = 1 and
  // calculate the rotation matrix and translation vector in case there is
  // a match.
  scalar R[3][3], t[3];
  con_chk_err(calculate_R_and_t(R, t, gid, nf, bid, nv, ndim, coord, tol, &c),
              "calculate_R_and_t failed.", &c);

  // Transform the points in bid = 1 set using R and t.
  transform_points(transformed_coord, nf, bid, nv, ndim, coord, R);

  // Globally number points:
  con_chk_err(
      number_points(new_gid, nf, nv, ndim, transformed_coord, tol, &c, &bfr),
      "number_points failed.", &c);

  // Update the global ids of the original points.
  update_global_ids(gid, nf, nv, new_gid, &c, &bfr);

  cleanup_before_return();
#undef cleanup_before_return
  return 0;
}

static int test_power_iteration_00(const scalar tol) {
  scalar A[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  scalar evecs[9], evals[3];
  power_iteration(evecs, evals, 3, A, tol);

  sint err = 0;
  err |= (fabs((evals[0] - 1.0) / 1.0) > tol);
  err |= (fabs((evals[1] - 1.0) / 1.0) > tol);
  err |= (fabs((evals[2] - 1.0) / 1.0) > tol);
  err |= (fabs(dot(evecs, evecs, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 3, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 6, evecs + 6, 3) - 1) > tol);
  err |= (fabs(dot(evecs, evecs + 3, 3)) > tol);
  err |= (fabs(dot(evecs, evecs + 6, 3)) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 6, 3)) > tol);

  return err;
}

static int test_power_iteration_01(const scalar tol) {
  scalar A[3][3] = {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}};

  scalar evecs[9], evals[3];
  power_iteration(evecs, evals, 3, A, tol);

  sint err = 0;
  err |= (fabs((evals[0] - 3.0) / 3.0) > tol);
  err |= (fabs((evals[1] - 2.0) / 2.0) > tol);
  err |= (fabs((evals[2] - 1.0) / 1.0) > tol);
  err |= (fabs(dot(evecs, evecs, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 3, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 6, evecs + 6, 3) - 1) > tol);
  err |= (fabs(dot(evecs, evecs + 3, 3)) > tol);
  err |= (fabs(dot(evecs, evecs + 6, 3)) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 6, 3)) > tol);

  return err;
}

static int test_power_iteration_02(const scalar tol) {
  scalar A[3][3] = {{5, -10, -5}, {2, 14, 2}, {-4, -8, 6}};
  scalar evecs[9], evals[3];
  power_iteration(evecs, evals, 3, A, tol);

  sint err = 0;
  err |= (fabs((evals[0] - 10.0) / 10.0) > tol);
  err |= (fabs((evals[1] - 10.0) / 10.0) > tol);
  err |= (fabs((evals[2] - 5.0) / 5.0) > tol);
  err |= (fabs(dot(evecs, evecs, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 3, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 6, evecs + 6, 3) - 1) > tol);
  err |= (fabs(dot(evecs, evecs + 3, 3)) > tol);
  err |= (fabs(dot(evecs, evecs + 6, 3)) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 6, 3)) > tol);

  return err;
}

static int test_power_iteration_03(const scalar tol) {
  scalar A[3][3] = {
      {1.0, 0.0, 0.0}, {0.0, 5.0e-01, 8.660254e-01}, {0.0, 0.0, 0.0}};

  scalar evals[3], evecs[9];
  power_iteration(evecs, evals, 3, A, tol);

  sint err = 0;
  err |= (fabs((evals[0] - 1.0)) > tol);
  err |= (fabs((evals[1] - 0.5)) > tol);
  err |= (fabs((evals[2] - 0.0)) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 3, 3) - 1) > tol);
  err |= (fabs(dot(evecs + 6, evecs + 6, 3) - 1) > tol);
  err |= (fabs(dot(evecs, evecs + 3, 3)) > tol);
  err |= (fabs(dot(evecs, evecs + 3, 3)) > tol);
  err |= (fabs(dot(evecs, evecs + 6, 3)) > tol);
  err |= (fabs(dot(evecs + 3, evecs + 6, 3)) > tol);

  return err;
}

static int test_transform_face_00(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++)
    for (sint j = 0; j < 3; j++)
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;

  scalar R[3][3], t[3];
  scalar error = transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {0.0, 0.0, 0.0};

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 9; i++) R_err[i] = R[0][i] - R_expected[0][i];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > 1e-8);
  err |= (normi(R_err, 9) > 1e-8);
  return err;
}

static int test_transform_face_01(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;
      face1[i][j] += 10.0;
    }
  }

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {10.0, 10.0, 10.0};

  scalar R[3][3], t[3];
  scalar error = transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 9; i++) R_err[i] = R[0][i] - R_expected[0][i];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > 1e-8);
  err |= (normi(R_err, 9) > 1e-8);
  return err;
}

static int test_transform_face_02(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;
      face1[i][j] += j;
    }
  }

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {0.0, 1.0, 2.0};

  scalar R[3][3], t[3];
  scalar error = transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 9; i++) R_err[i] = R[0][i] - R_expected[0][i];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > 1e-8);
  err |= (normi(R_err, 9) > 1e-8);
  return err;
}

static int test_transform_face_03(const scalar tol) {
  scalar face0[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}};

  // Rotate by 60 degrees around the x-axis and translate by (1, 2, 3).
  scalar theta = M_PI / 3.0;
  scalar R_orig[3][3] = {{1.0, 0.0, 0.0},
                         {0.0, cos(theta), -sin(theta)},
                         {0.0, sin(theta), cos(theta)}};
  scalar t_orig[3] = {1.0, 2.0, 3.0};

  scalar face1[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = 0;
      for (sint k = 0; k < 3; k++) face1[i][j] += R_orig[j][k] * face0[i][k];
      face1[i][j] += t_orig[j];
    }
  }

  scalar R[3][3], t[3];
  scalar error = transform_face(R, t, face1, face0, 4, 3, tol);

  sint err = 0;

  scalar errs[12];
  for (sint i = 0; i < 12; i++) errs[i] = R[0][i] - R_orig[0][i];
  err |= (normi(errs, 12) > 1e-8);

  for (sint i = 0; i < 3; i++) errs[i] = t[i] - t_orig[i];
  err |= (normi(errs, 3) > 1e-8);

  return err;
}

#define chk_test(call, count, c)                                               \
  {                                                                            \
    sint err = (call);                                                         \
    sint wrk;                                                                  \
    comm_allreduce(&(c), gs_int, gs_max, &err, 1, &wrk);                       \
    if (err) {                                                                 \
      if ((c).id == 0) fprintf(stderr, #call " failed.\n");                    \
      (count)++;                                                               \
    }                                                                          \
  }

int test_automatic_periodic_face_match(slong *const gid, uint nf,
                                       const sint *const bid, sint nv,
                                       sint ndim, const scalar *const coord,
                                       scalar tol, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  sint errs = 0;
  chk_test(test_power_iteration_00(tol), errs, c);
  chk_test(test_power_iteration_01(tol), errs, c);
  chk_test(test_power_iteration_02(tol), errs, c);
  chk_test(test_power_iteration_03(tol), errs, c);

  chk_test(test_transform_face_00(tol), errs, c);
  chk_test(test_transform_face_01(tol), errs, c);
  chk_test(test_transform_face_02(tol), errs, c);
  chk_test(test_transform_face_03(tol), errs, c);

  comm_free(&c);
  return errs;
}

#undef distance2D
#undef distance3D
