#ifndef _GENMAP_GENCON_IMPL_H_
#define _GENMAP_GENCON_IMPL_H_

#include <gencon.h>
#include <genmap-impl.h>

/* Boundary condition types */
#define GC_PERIODIC "P  "

/* Header lengths */
#define GC_RE2_HEADER_LEN 80
#define GC_CO2_HEADER_LEN 132

extern int NEIGHBOR_MAP[GC_MAX_VERTICES][GC_MAX_NEIGHBORS];
extern int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES];
extern int PRE_TO_SYM_FACE[GC_MAX_FACES];

#define sqrDiff(x, y) (((x) - (y)) * ((x) - (y)))
#define distance2D(a, b)                                                       \
  (sqrDiff((a).x[0], (b).x[0]) + sqrDiff((a).x[1], (b).x[1]))
#define distance3D(a, b) (distance2D(a, b) + sqrDiff((a).x[2], (b).x[2]))

#define READ_T(coords, buf, T, nVertex)                                        \
  do {                                                                         \
    memcpy((coords), buf, sizeof(T) * nVertex);                                \
  } while (0)

#define WRITE_INT(dest, val)                                                   \
  do {                                                                         \
    memcpy(dest, &(val), sizeof(int));                                         \
  } while (0)

// TODO: Use rsb_element here
struct Point_private {
  GenmapScalar dx;
  GenmapScalar x[GC_MAX_DIM];
  uint proc;
  uint origin;
  int ifSegment;
  ulong sequenceId;
  ulong elementId;
  ulong globalId;
};

struct Face_private {
  struct Point_private vertex[GC_MAX_FACE_VERTICES];
};

struct Element_private {
  struct Point_private vertex[GC_MAX_VERTICES];
};

struct Boundary_private {
  ulong elementId;
  ulong faceId;
  struct Face_private face;
  uint proc;
  long bc[2];
};

struct Mesh_private {
  ulong nelgt;
  ulong nelgv;
  uint nelt;
  int nDim;
  int nVertex;
  int nNeighbors;
  struct array elements;
  struct array boundary;
};

int transferBoundaryFaces(Mesh mesh, struct comm *c);

#endif
