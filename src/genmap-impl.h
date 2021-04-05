#ifndef _GENMAP_IMPL_H_
#define _GENMAP_IMPL_H_

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <genmap.h>

#include <genmap-elements.h>
#include <genmap-metrics.h>
#include <genmap-multigrid.h>
#include <genmap-partition.h>

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

/* Dump data */
int GenmapFiedlerDump(const char *fname, genmap_handle h, slong start,
                      struct comm *c);

int GenmapVectorDump(const char *fname, GenmapScalar *y, uint size,
                     struct comm *c);

int GenmapCentroidDump(const char *fname, genmap_handle h, sint gid,
                       struct comm *c);

int GenmapElementIdDump(const char *fname, genmap_handle h, struct comm *c);

#endif
