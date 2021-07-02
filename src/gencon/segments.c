#include <math.h>
#include <stdlib.h>

#include <gencon-impl.h>
#include <sort.h>

static void initSegments(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  uint i;
  for (i = 0; i < nPoints; i++)
    points[i].ifSegment = 0;

  /* First rank with nPoints > 0 set ifSegment = 1 */
  sint rank = c->id;
  if (nPoints == 0)
    rank = c->np;

  sint buf[2];
  comm_allreduce(c, gs_int, gs_min, &rank, 1, buf);

  if (c->id == rank)
    points[0].ifSegment = 1;
}

static int sendLastElement(struct array *arr, Mesh mesh, struct comm *c) {
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
  /* Parallel sort: first by globalId then by x[dim] */
  switch (dim) {
    case 0:
      parallel_sort_2(struct Point_private, &mesh->elements, globalId, 1, x[0],
                      gs_scalar, bin_sort, 1, c, bfr);
      break;
    case 1:
      parallel_sort_2(struct Point_private, &mesh->elements, globalId, 1, x[1],
                      gs_scalar, bin_sort, 1, c, bfr);
      break;
    case 2:
      parallel_sort_2(struct Point_private, &mesh->elements, globalId, 1, x[2],
                      gs_scalar, bin_sort, 1, c, bfr);
      break;
    default:
      break;
  }

  /* Identify the starting points of segments after sort */
  struct Point_private *pnt = mesh->elements.ptr;
  pnt->ifSegment = 1;

  struct array arr;
  sendLastElement(&arr, mesh, c);

  if (c->id > 0) {
    struct Point_private *pnt1 = arr.ptr;
    if (pnt->globalId == pnt1->globalId)
      pnt->ifSegment = 0;
  }

  uint npts = mesh->elements.n;
  ulong prev = pnt[0].globalId;
  uint i;
  for (i = 1; i < npts; i++) {
    if (pnt[i].globalId != prev) {
      pnt[i].ifSegment = 1;
      prev = pnt[i].globalId;
    } else
      pnt[i].ifSegment = 0;
  }

  array_free(&arr);

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
  sendLastElement(&arr, mesh, c);

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
#endif

slong countSegments(Mesh mesh, int verbose, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;

  sint count = 0, i;
  for (i = 0; i < nPoints; i++) {
    if (points[i].ifSegment > 0) {
      count++;
      if (verbose > 0)
        printf("nid = %d point = %.10e %.10e %.10e\n", c->id, points[i].x[0],
               points[i].x[1], points[i].x[2]);
    }
  }

  slong buf[2][1];
  slong in = count;
  comm_allreduce(c, gs_long, gs_add, &in, 1, buf);

  return in;
}

static int findBin(Mesh mesh, struct comm *c) {
  uint nPoints = mesh->elements.n;
  Point points = mesh->elements.ptr;


  slong buf[2][1], out[2][1];
  sint in = (nPoints > 0) ? points[0].ifSegment : 0;
  comm_scan(out, c, gs_int, gs_add, &in, 1, buf);

  return out[0][0] + in;
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

int findUniqueVertices(Mesh mesh, struct comm *c, GenmapScalar tol, int verbose,
                       buffer *bfr) {
  initSegments(mesh, c);

  GenmapScalar tolSquared = tol * tol;
  int nDim = mesh->nDim;
  int merge = 1;

  int t, d;
  for (t = 0; t < nDim; t++) {
    for (d = 0; d < nDim; d++) {
      numberSegments(mesh, c);

      sortSegments(mesh, c, d, bfr);

      findSegments(mesh, c, d, tolSquared);

      slong n_seg = countSegments(mesh, t == 0 && d == 0, c);
      if (c->id == 0)
        printf("\tlocglob: %d %d %lld\n", t + 1, d + 1, n_seg);
    }
  }

  return 0;
}
