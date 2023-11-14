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

  for (i = 0; i < npoints; i++)
    points[i].globalId += start;

  return 0;
}

static int renumber_periodic_vertices(Mesh mesh, struct comm *c,
                                      struct array *matches, buffer *bfr) {
  uint size1 = mesh->elements.n, size2 = matches->n;
  slong *mids = tcalloc(slong, size1 + 2 * size2),
        *mnew = tcalloc(slong, size1 + 2 * size2),
        *mcur = tcalloc(slong, size1);

  struct point_t *pe = (struct point_t *)mesh->elements.ptr;
  for (uint i = 0; i < size1; i++)
    mids[i] = pe[i].globalId;
  struct mpair_t *pm = (struct mpair_t *)matches->ptr;
  for (uint i = 0; i < size2; i++)
    mids[size1 + i] = pm[i].orig;
  struct gs_data *gsh = gs_setup(mids, size1 + size2, c, 0, gs_pairwise, 0);

  for (uint i = 0; i < size1; i++)
    mnew[i] = pe[i].globalId;
  for (uint i = 0; i < size2; i++)
    mnew[size1 + i] = pm[i].min;
  gs(mnew, gs_long, gs_min, 0, gsh, bfr);

  sint changed, wrk;
  do {
    for (uint i = 0; i < size1; i++)
      mcur[i] = mnew[i];
    for (uint i = 0; i < size2; i++)
      mids[size1 + size2 + i] = -mnew[size1 + i];
    struct gs_data *gsh1 =
        gs_setup(mids, size1 + 2 * size2, c, 0, gs_pairwise, 0);

    gs(mnew, gs_long, gs_min, 0, gsh1, bfr);
    gs_free(gsh1);

    for (uint i = 0; i < size2; i++)
      mnew[size1 + i] = mnew[size1 + size2 + i];
    gs(mnew, gs_long, gs_min, 0, gsh, bfr);

    changed = 0;
    for (uint i = 0; i < size1; i++)
      changed += (mnew[i] != mcur[i]);
    comm_allreduce(c, gs_int, gs_max, &changed, 1, &wrk);
  } while (changed);

  for (uint i = 0; i < size1; i++)
    pe[i].globalId = mcur[i];

  gs_free(gsh);
  free(mids), free(mnew), free(mcur);

  return 0;
}

static int match_connected_periodic_pair(Mesh mesh, BoundaryFace f_,
                                         BoundaryFace g_,
                                         struct array *matches) {
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
    array_cat(struct mpair_t, matches, &m, 1);
  }

  return 0;
}

static sint find_connected_periodic_face(Mesh mesh, const sint i,
                                         const sint *matched) {
  sint size = mesh->boundary.n;
  BoundaryFace ptr = mesh->boundary.ptr;

  ulong elementId = ptr[i].bc[0];
  ulong faceId = ptr[i].bc[1];

  uint start = i + 1, end = size;
  assert(start > 0 && "start must be larger than 0");
  assert(start < end && "start must be less than size");

  while (start < end) {
    uint mid = (start + end) / 2;

    // Found a match !
    if (elementId == ptr[mid].elementId) {
      uint j = mid;
      do {
        if (faceId == ptr[j].faceId)
          return j;
        j--;
      } while (j > 0 && elementId == ptr[j].elementId);

      j = mid + 1;
      while (j < size && elementId == ptr[j].elementId) {
        if (faceId == ptr[j].faceId)
          return j;
        j++;
      }
      // Ooops!
      return UINT_MAX;
    }

    if (elementId < ptr[mid].elementId) {
      end = mid;
      continue;
    }

    // If elementId > ptr[mid].elementId:
    start = mid;
  }

  return UINT_MAX;
}

static int find_connected_periodic_faces(Mesh mesh, struct array *matches,
                                         buffer *bfr) {
  sint size = mesh->boundary.n;
  BoundaryFace ptr = mesh->boundary.ptr;

  sarray_sort(struct boundary_t, ptr, size, elementId, 1, bfr);

  sint *matched = tcalloc(sint, size);

  for (sint i = 0; i < size - 1; i++) {
    if (matched[i])
      continue;
    uint j = find_connected_periodic_face(mesh, i, matched);
    assert(j < UINT_MAX && "No matching face found !");
    matched[i] = matched[j] = 1;
    match_connected_periodic_pair(mesh, &ptr[i], &ptr[j], matches);
  }
  // Sanity check:
  for (uint i = 0; i < size; i++)
    assert(matched[i]);

  free(matched);

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
  if (bSize == 0)
    return 0;

  Point ePtr = mesh->elements.ptr;
  sint eSize = mesh->elements.n;
  if (eSize == 0)
    return 0;

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
    while (k < eSize && ePtr[k].elementId < bPtr[i].elementId)
      k += nv;
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

  struct array matches;
  array_init(struct mpair_t, &matches, 1024);

  parrsb_print(c, verbose, "\t\t%s ...", functions[2]);
  find_connected_periodic_faces(mesh, &matches, bfr);

  parrsb_print(c, verbose, "\t\t%s ...", functions[3]);
  renumber_periodic_vertices(mesh, c, &matches, bfr);
  array_free(&matches);

  parrsb_print(c, verbose, "\t\t%s ...", functions[4]);
  compress_periodic_vertices(mesh, c, bfr);

  parrsb_print(c, verbose, "\t\t%s ...", functions[5]);
  send_back(mesh, c, bfr);

  return 0;
}

#undef distance2D
#undef distance3D
