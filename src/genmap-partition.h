#ifndef _GENMAP_PARTITION_H_
#define _GENMAP_PARTITION_H_

#include <genmap-elements.h>
#include <genmap.h>

/* RCB/RIB */
int rcb(struct comm *ci, struct array *elements, int ndim, buffer *bfr);
int rib(struct comm *ci, struct array *elements, int ndim, buffer *bfr);

/* Eigen */
int GenmapTQLI(genmap_handle h, genmap_vector diag, genmap_vector upper,
               genmap_vector **eVec, genmap_vector *eVal);

int genmap_inverse_power(double *y, int N, double *A, int verbose);

int genmap_power(double *y, int N, double *A, int verbose);

/* Laplacian */
void genmap_number_faces_and_edges(genmap_handle h, struct comm *c);

int genmap_laplacian_csr_init(genmap_handle h, struct comm *c);
int genmap_laplacian_csr(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

int genmap_laplacian_gs_init(genmap_handle h, struct comm *c);
int genmap_laplacian_gs(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

int genmap_laplacian_init(genmap_handle h, struct comm *c);
int genmap_laplacian(genmap_handle h, GenmapScalar *u, GenmapScalar *v);

/* Lanczos */
int GenmapLanczosLegendary(genmap_handle h, struct comm *c, genmap_vector f,
                           GenmapInt niter, genmap_vector **rr,
                           genmap_vector diag, genmap_vector upper);

int GenmapLanczos(genmap_handle h, struct comm *c, genmap_vector init,
                  GenmapInt iter, genmap_vector **q, genmap_vector alpha,
                  genmap_vector beta);

/* Fiedler */
int GenmapFiedlerLanczos(genmap_handle h, struct comm *c, int maxIter,
                         int global);

int GenmapFiedlerRQI(genmap_handle h, struct comm *c, int maxIter, int global);

/* Components and repair */
sint get_components(sint *component, struct rsb_element *elements,
                    struct comm *c, buffer *buf, uint nelt, uint nv);

void split_and_repair_partitions(genmap_handle h, struct comm *lc, int level);

/* Matrix inverse */
void matrix_inverse(int N, double *A);

#endif
