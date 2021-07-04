#include <math.h>
#include <stdlib.h>

#include <gencon-impl.h>
#include <sort.h>

static void initSegment(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  uint i;
  for (i = 0; i < nPoints; i++) {
    points[i].ifSegment = 0;
    points[i].globalId = 1;
  }

  /* First rank with nPoints > 0 set ifSegment = 1 */
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  sint buf[2];
  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sendLastPoint(struct array *arr, Mesh mesh, struct comm *c) {
  Point pts = mesh->elements.ptr;
  sint npts = mesh->elements.n;

  struct Point_private lastp = pts[npts - 1];
  lastp.proc = (c->id + 1) % c->np;

  array_init(struct Point_private, arr, 1);
  array_cat(struct Point_private, arr, &lastp, 1);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct Point_private, arr, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

static int sortSegments(Mesh mesh, struct comm *c, int dim, buffer *bfr) {
  if (c->np > 1) {
    /* Parallel sort -- we haven't localized the problem yet */
    switch (dim) {
    case 0:
      parallel_sort(struct Point_private, &mesh->elements, x[0], gs_scalar,
                    bin_sort, 1, c, bfr);
      break;
    case 1:
      parallel_sort(struct Point_private, &mesh->elements, x[1], gs_scalar,
                    bin_sort, 1, c, bfr);
      break;
    case 2:
      parallel_sort(struct Point_private, &mesh->elements, x[2], gs_scalar,
                    bin_sort, 1, c, bfr);
      break;
    default:
      break;
    }

    initSegment(mesh, c);
  } else {
    /* Local sort: Segments are local */
  }

  return 0;
}

static int findSegments(Mesh mesh, struct comm *c, int i,
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

  struct array arr;
  sendLastPoint(&arr, mesh, c);

  if (c->id > 0) {
    struct Point_private *lastp = arr.ptr;
    GenmapScalar d = sqrDiff(lastp->x[i], pts->x[i]);
    GenmapScalar dx = min(lastp->dx, pts->dx) * tolSquared;
    if (d > dx)
      pts->ifSegment = 1;
  }

  array_free(&arr);

  return 0;
}

#if 0
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

int numberSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  uint count = 0;
  uint i;
  for (i = 0; i < nPoints; i++)
    if (points[i].ifSegment > 0)
      count++;

  slong buf[2][1], out[2][1];
  slong in = count;
  comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
  slong start = out[0][0];

  for (i = 0; i < nPoints; i++) {
    if (points[i].ifSegment > 0)
      start++;
    points[i].globalId = start;
  }

  return 0;
}
#endif

slong countSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++) {
    if (points[i].ifSegment > 0) {
      count++;
    }
  }

  slong buf[2][1];
  slong in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);

  return in;
}

static int rearrangeSegments(Mesh mesh, struct comm *c) {
  struct comm seg;
  comm_init(&seg, c);

  while (seg.np > 1 || countSegments(mesh, &seg) > 1) {
    uint nPoints = mesh->elements.n;
    Point points = mesh->elements.ptr;

    /* comm_scan */
    slong out[2][1], buf[2][1], in[1];
    in[0] = nPoints;
    comm_scan(out, c, gs_long, gs_add, in, 1, buf);
    slong start = out[0][0];
    slong nelg = out[1][0];

    double min = DBL_MAX;
    int inc_proc;

    uint i, index;
    for (i = 0; i < nPoints; i++) {
      if (points[i].ifSegment > 0) {
        double frac0 = fabs((start + i)/nelg - (c->id + 0.0)/c->np);
        double frac1 = fabs((start + i)/nelg - (c->id + 1.0)/c->np);
        if (frac0 < min) {
          inc_proc = 0;
          min = frac0;
          index = i;
        }
        if (frac1 < min) {
          inc_proc = 1;
          min = frac1;
          index = i;
        }
      }
    }

    double dbuf[2];
    double ming = min;
    comm_allreduce(c, gs_double, gs_min, &ming, 1, dbuf);

    sint rankg = -1;
    if (fabs(ming - min) < 1e-15)
      rankg = c->id + inc_proc;
    comm_allreduce(c, gs_int, gs_max, &rankg, 1, buf);

    int bin = 1;
    if (c->id < rankg)
      bin = 0;

    /* Transfer and split */
  }

  comm_free(c);
  comm_init(c, &seg);
  comm_free(&seg);
}

int findUniqueVertices(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose,
                       buffer *bfr) {
  GenmapScalar tolSquared = tol * tol;
  int nDim = mesh->nDim;

  struct comm seg;
  comm_init(&seg, c);

  int t, d;
  for (t = 0; t < nDim; t++) {
    for (d = 0; d < nDim; d++) {
      sortSegments(mesh, &seg, d, bfr);
      findSegments(mesh, &seg, d, tolSquared);

      slong n_seg = countSegments(mesh, c);
      if (c->id == 0)
        printf("\tlocglob: %d %d %lld\n", t + 1, d + 1, n_seg);

      rearrangeSegments(mesh, &seg);
    }
  }

  comm_free(&seg);

  return 0;
}
