#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef GENMAP_DEBUG
#include <stdio.h>
#endif

#include <genmap-multigrid-precon.h>
#include <genmap.h>

#define GENMAP_RCB_ELEMENT 0
#define GENMAP_RSB_ELEMENT 1

#define GENMAP_ALIGN 32

#define GENMAP_SP_TOL 1e-05
#define GENMAP_DP_TOL 1e-12
#define GENMAP_TOL GENMAP_DP_TOL

#define GENMAP_MIN(a, b) ((a) < (b) ? (a) : (b))
#define GENMAP_MAX(a, b) ((a) > (b) ? (a) : (b))

#define MAXDIM 3 /* Maximum dimension of the mesh */
#define MAXNV 8  /* Maximum number of vertices per element */

/*
 Preprocessor Corner notation:      Symmetric Corner notation:

         4+-----+3    ^ s                    3+-----+4    ^ s
         /     /|     |                      /     /|     |
        /     / |     |                     /     / |     |
      8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
       |     | /     /                     |     | /     /
       |     |/     /                      |     |/     /
      5+-----+6    t                      5+-----+6    t



                   i) Preprocessor notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    3   / |     |
                             4--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     1

                  ii) Symmetric notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    4   / |     |
                             1--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     3

   EFACE(IFACE)  - Given face number IFACE in symmetric notation,
                   returns preprocessor notation face number.

   EFACE1(IFACE) - Given face number IFACE in preprocessor notation,
                   returns symmetric notation face number.

The following variables all take the symmetric notation of IFACE
as arguments:

   ICFACE(i,IFACE) - Gives the 4 vertices which reside on face IFACE
                     as depicted below, e.g. ICFACE(i,2)=2,4,6,8.

                      3+-----+4    ^ Y
                      /  2  /|     |
   Edge 1 extends    /     / |     |
     from vertex   7+-----+8 +2    +----> X
     1 to 2.        |  4  | /     /
                    |     |/     /
                   5+-----+6    Z
                       3

   IEDGFC(i,IFACE) - Gives the 4 edges which border the face IFACE
                     Edge numbering is as follows:
                        Edge = 1,2,3,4     run in +r direction
                        Edge = 5,6,7,8     run in +s direction
                        Edge = 9,10,11,12  run in +t direction

                     Ordering of each edge is such that a monotonically
                     increasing sequence of vertices is associated with
                     the start point of a corresponding set of
                     monotonically increasing edge numbers, e.g.,

   ICEDG(i,IEDGE)  - Gives 3 variables for determining the stride along
                     a given edge, IEDGE;  i=1 gives the starting vertex
                                           i=2 gives the stopping vertex
                                           i=3 gives the stride size.

*/

/* Upper bounds for faces, vertices, and edges */
#define GC_MAX_VERTICES 8
#define GC_MAX_EDGES 12
#define GC_MAX_FACES 6

/* Upper bounds for faces */
#define GC_MAX_FACE_VERTICES 4
#define GC_MAX_EDGE_VERTICES 2

/* Upper bound for number of dimensions */
#define GC_MAX_DIM 3
#define GC_MAX_NEIGHBORS 3

/* 1 - indexed */
extern int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
extern int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];

/* 0 - indexed */
extern int edges3D[GC_MAX_EDGES][GC_MAX_EDGE_VERTICES];

/* rcb_element is used for rcb and rib */
struct rcb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
};

/* rsb_element should be a superset of rcb_element */
struct rsb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[MAXDIM];
  GenmapScalar fiedler;
  GenmapLong vertices[8];
  GenmapInt part;
  GenmapULong globalId0;
};

int rcb(struct comm *ci, struct array *elements, int ndim, buffer *bfr);
int rib(struct comm *ci, struct array *elements, int ndim, buffer *bfr);

struct genmap_handle_private {
  genmap_comm global;
  genmap_comm local;

  GenmapLong nel;
  GenmapLong Nnodes;
  GenmapLong start;
  int nv;

  struct array *elements;

  /* Weighted Laplacian */
  struct gs_data *gsw;
  GenmapScalar *weights;
  buffer buf;

  /* Un-weighted Laplacian */
  struct gs_data *gs;
  GenmapScalar *diagonal;

  struct gs_data *gsh;
  csr_mat M;
  GenmapScalar *b;

  parRSB_options *options;
};

struct genmap_vector_private {
  GenmapInt size;
  GenmapScalar *data;
};

int GenmapMallocArray(size_t n, size_t unit, void *p);
int GenmapCallocArray(size_t n, size_t unit, void *p);
int GenmapReallocArray(size_t n, size_t unit, void *p);
int GenmapFree(void *p);

#define GenmapMalloc(n, p) GenmapMallocArray((n), sizeof(**(p)), p)
#define GenmapCalloc(n, p) GenmapCallocArray((n), sizeof(**(p)), p)
#define GenmapRealloc(n, p) GenmapReallocArray((n), sizeof(**(p)), p)

/* Genmap Metrics */
typedef enum {
  RCB,
  WEIGHTEDLAPLACIANSETUP,
  FIEDLER,
  NFIEDLER,
  FIEDLERSORT,
  BISECTANDREPAIR,
  LANCZOS,
  NLANCZOS,
  WEIGHTEDLAPLACIAN,
  TQLI,
  LAPLACIANSETUP,
  LAPLACIANGSSETUP,
  FINDNBRS,
  FIRSTHALF,
  SECONDHALF,
  CSRMATSETUP,
  CSRTOPSETUP,
  PRECONDSETUP,
  RQI,
  NRQI,
  PROJECT,
  NPROJECT,
  GRAMMIAN,
  LAPLACIAN,
  VCYCLE,
  END
} metric;

void metric_init();
void metric_acc(metric m, double count);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_print(struct comm *c);
void metric_finalize();

/* Laplacian and gencon */
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

/* Components */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv);

void split_and_repair_partitions(genmap_handle h, struct comm *lc, int level);
/* Matrix inverse */
void matrix_inverse(int N, double *A);

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, slong start,
                      struct comm *c);
int GenmapVectorDump(const char *fname, GenmapScalar *y, uint size,
                     struct comm *c);
int GenmapCentroidDump(const char *fname, genmap_handle h, sint g_id,
                       struct comm *c);
int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c);

#endif
