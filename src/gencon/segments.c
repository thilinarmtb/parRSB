#include <math.h>
#include <stdlib.h>

#include <gencon-impl.h>
#include <sort.h>

static void initSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  /* Initialize globalId and ifSegment */
  slong out[2][1], buf[2][1], in[1];
  in[0] = nPoints;
  comm_scan(out, c, gs_long, gs_add, in, 1, buf);
  slong start = out[0][0];

  uint i;
  for (i = 0; i < nPoints; i++) {
    points[i].ifSegment = 0;
    points[i].globalId = start + i;
  }

  /* First rank with nPoints > 0 set ifSegment = 1 */
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sortLocalSegments(Mesh mesh, int dim, buffer *bfr) {
  sint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint s = 0, e;
  while (s < nPoints) {
    // find the length of the segment
    for (e = s + 1; e < nPoints && points[e].ifSegment == 0; e++)
      ;

    // sort start to end based on dim
    switch (dim) {
    case 0:
      sarray_sort_2(struct Point_private, &points[s], e - s, x[0], 3,
                    sequenceId, 1, bfr);
      break;
    case 1:
      sarray_sort_2(struct Point_private, &points[s], e - s, x[1], 3,
                    sequenceId, 1, bfr);
      break;
    case 2:
      sarray_sort_2(struct Point_private, &points[s], e - s, x[2], 3,
                    sequenceId, 1, bfr);
      break;
    default:
      break;
    }

    sint i, sum = 0;
    for (i = s; i < e; i++) {
      sum += points[i].ifSegment;
      points[i].ifSegment = 0;
    }

    if (sum > 0)
      points[s].ifSegment = 1;

    s = e;
  }

  return 0;
}

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  sortLocalSegments(mesh, dim, bfr);
}

static int findLocalSegments(Mesh mesh, struct comm *c, int i,
                             GenmapScalar tolSquared) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;
  int nDim = mesh->nDim;

  sint j;
  for (j = 1; j < npts; j++) {
    GenmapScalar d = sqrDiff(pts[j].x[i], pts[j - 1].x[i]);

    GenmapScalar dx = min(pts[j].dx, pts[j - 1].dx) * tolSquared;

    if (d > dx)
      pts[j].ifSegment = 1;
  }

  sint rank = c->id;
  sint size = c->np;

  struct Point_private lastp = pts[npts - 1];
  lastp.proc = (rank + 1) % size;

  struct array arr;
  array_init(struct Point_private, &arr, 1);
  array_cat(struct Point_private, &arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &arr, proc, 1, &cr);
  crystal_free(&cr);

  uint n = arr.n;
  assert(n == 1);
  lastp = ((struct Point_private *)arr.ptr)[0];

  if (rank > 0) {
    GenmapScalar d = sqrDiff(lastp.x[i], pts->x[i]);
    GenmapScalar dx = min(lastp.dx, pts->dx) * tolSquared;
    if (d > dx)
      pts->ifSegment = 1;
  }

  array_free(&arr);

  return 0;
}

static int findSegments(Mesh mesh, struct comm *c, int i, GenmapScalar tolSquared) {
  findLocalSegments(mesh, c, i, tolSquared);
  return 0;
}

static int mergeSegments(Mesh mesh, struct comm *c, buffer *bfr) {
  uint npoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  int n, ifseg = 0;
  uint sendn = 0;
  for (n = 0; n < npoints; n++)
    if (points[n].ifSegment == 1) {
      ifseg = 1;
      sendn = n;
      break;
    }

  sint out[2][1], buf[2][1], in[1];
  in[0] = ifseg * (c->id + 1);
  comm_scan(out, c, gs_int, gs_max, in, 1, buf);
  sint rank = out[0][0] - 1;

  if (rank == -1) {
    if (c->id > 0)
      rank = c->id - 1;
    else if (c->id == 0)
      rank = 0;
  }

  // If rank > 0, send i = 0,... n-1 where points[i].ifSegment == 0 to
  // rank with previous ifSegment == 1
  for (n = 0; n < sendn; n++)
    points[n].proc = rank;
  for (; n < npoints; n++)
    points[n].proc = c->id;

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, &mesh->elements, proc, 0, &cr);
  crystal_free(&cr);

  sarray_sort(struct Point_private, mesh->elements.ptr, mesh->elements.n,
              globalId, 1, bfr);

  return 0;
}

slong countSegments(Mesh mesh, int verbose, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++)
    if (points[i].ifSegment > 0) {
      count++;
      if (verbose > 0)
        printf("nid = %d point = %.10e %.10e %.10e\n", c->id, points[i].x[0],
               points[i].x[1], points[i].x[2]);
    }

  slong in, buf[2][1];
  in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);
  return in;
}

#define sort_by_coord(mesh, c, xa, xb, xc, bfr)                                \
  do {                                                                         \
    parallel_sort(struct Point_private, &(mesh->elements), x[xa], gs_scalar,   \
                  bin_sort, 0, c, bfr);                                        \
    uint nPoints = mesh->elements.n;                                           \
    Point points = mesh->elements.ptr;                                         \
                                                                               \
    int nDim = mesh->nDim;                                                     \
    if (nDim == 13)                                                            \
      sarray_sort_3(struct Point_private, points, nPoints, x[xa], 3, x[xb], 3, \
                    x[xc], 3, bfr);                                            \
    else if (nDim == 12)                                                       \
      sarray_sort_2(struct Point_private, points, nPoints, x[xa], 3, x[xb], 3, \
                    bfr);                                                      \
    sarray_sort_2(struct Point_private, points, nPoints, x[xa], 3, sequenceId, \
                  1, bfr);                                                     \
  } while (0)

int findUniqueVertices(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose,
                       buffer *bfr) {
  struct comm comm_seg;
  comm_dup(&comm_seg, c);

  initSegments(mesh, &comm_seg);

  GenmapScalar tolSquared = tol * tol;
  int nDim = mesh->nDim;
  int merge = 1;

  int t, d;
  for (t = 0; t < nDim; t++) {
    for (d = 0; d < nDim; d++) {
      sortSegments(mesh, &comm_seg, d, bfr);
      findSegments(mesh, &comm_seg, d, tolSquared);

      slong n_seg = countSegments(mesh, t == 0 && d == 0, &comm_seg);
      if (comm_seg.id == 0)
        printf("\tlocglob: %d %d %lld\n", t + 1, d + 1, n_seg);

      if (merge > 0) {
        mergeSegments(mesh, &comm_seg, bfr);
        merge = 0;
      }

      comm_free(&comm_seg);
      //genmap_comm_split(c, bin, c->id, &comm_seg);
    }
  }

  comm_free(&comm_seg);

  return 0;
}
