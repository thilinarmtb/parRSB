#include <stdio.h>
#include <stdlib.h>

#include <genmap-impl.h>

/* 1 - indexed */
int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{1, 5, 7, 3}, {2, 4, 8, 6},
                                                   {1, 2, 6, 5}, {3, 7, 8, 4},
                                                   {1, 3, 4, 2}, {5, 6, 8, 7}};

int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES] = {{3, 1, 0, 0}, {2, 4, 0, 0},
                                                   {1, 2, 0, 0}, {4, 3, 0, 0},
                                                   {0, 0, 0, 0}, {0, 0, 0, 0}};

/* 0 - indexed */
int edges3D[GC_MAX_EDGES][GC_MAX_EDGE_VERTICES] = {
    {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 4}, {2, 6},
    {3, 7}, {1, 5}, {0, 2}, {1, 3}, {4, 6}, {5, 7}};

int genmap_init(genmap_handle *h_, comm_ext ce, parRSB_options *options) {
  GenmapMalloc(1, h_);
  genmap_handle h = *h_;

  GenmapMalloc(1, &h->global);
  comm_init(h->global, ce);

  GenmapMalloc(1, &h->local);
  comm_init(h->local, ce);

  /* GS based Laplacian */
  h->gs = NULL;
  h->diagonal = NULL;

  /* CSR based Laplacian */
  h->M = NULL;

  buffer_init(&h->buf, 1024);

  h->options = options;

  h->start = 0;
  h->nel = 0;

  return 0;
}

int genmap_finalize(genmap_handle h) {
  buffer_free(&h->buf);

  /* GS based Laplacian */
  if (h->diagonal != NULL)
    GenmapFree(h->diagonal);
  if (h->gs)
    gs_free(h->gs);

  /* CSR based Laplacian */
  if (h->M)
    csr_mat_free(h->M);

  if (h->global != NULL) {
    comm_free(h->global);
    GenmapFree(h->global);
  }

  if (h->local != NULL) {
    comm_free(h->local);
    GenmapFree(h->local);
  }

  GenmapFree(h);

  return 0;
}

void *genmap_get_elements(genmap_handle h) {
  return (struct rsb_element *)h->elements->ptr;
}
void genmap_set_elements(genmap_handle h, struct array *elements) {
  h->elements = elements;
}

int genmap_get_nvertices(genmap_handle h) { return h->nv; }
void genmap_set_nvertices(genmap_handle h, int nv) { h->nv = nv; }

GenmapInt genmap_get_nel(genmap_handle h) { return h->elements->n; }
GenmapULong genmap_get_partition_nel(genmap_handle h) { return h->nel; }

int GenmapMallocArray(size_t n, size_t unit, void *p) {
  int ierr = posix_memalign((void **)p, GENMAP_ALIGN, n * unit);
  if (ierr)
    printf("GenmapMallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  return ierr;
}

int GenmapCallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = calloc(n, unit);
  if (n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapCallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapReallocArray(size_t n, size_t unit, void *p) {
  int ierr = 0;
  *(void **)p = realloc(*(void **)p, n * unit);
  if (n && unit && !*(void **)p) {
    ierr = 1;
    printf("GenmapReallocArray Failed: %s:%d\n", __FILE__, __LINE__);
  }
  return ierr;
}

int GenmapFree(void *p) {
  free(p);
  p = NULL;
  return 0;
}
