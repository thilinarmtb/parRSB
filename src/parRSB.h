#ifndef _PARRSB_H_
#define _PARRSB_H_

#include <mpi.h>

typedef struct {
  /* General options */
  int global_partitioner; // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int local_partitioner;  // -1 - None, 0 - RSB, 1 - RCB, 2 - RIB (Default: -1)
  int debug_level;        // 0, 1, 2, .. etc (Default: 0)
  int print_timing_info;  // 0 or 1 (Default: 0)

  /* RSB specific */
  int rsb_algo;                     // 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_prepartition;             // 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_grammian;                 // 0 or 1 (Default: 1)
  int rsb_laplacian_implementation; // 1 - gather-scatter, 2 - CSR, 3 - both
                                    // (Default: 1)
  int rsb_laplacian_weighted;       // 0 - unweighted, 1 - weighted (Default: 1)
} parRSB_options;

extern parRSB_options parrsb_default_options;

/*
 * part = [nel], out,
 * seq = [nel], out,
 * vtx = [nel x nv], in,
 * coord = [nel x nv x ndim], in,
 * nel = in,
 * nv = in,
 * options = in */
int parRSB_partMesh(int *part, int *seq, long long *vtx, double *coord, int nel,
                    int nv, parRSB_options *options, MPI_Comm comm);

/* vtx = [nelt, nv], out
 * coord = [nelt, nv, ndim], in (vertices are orders in preprocessor ordering),
 * nel = in,
 * ndim = in,
 * periodicInfo = [nPeriodicFaces x 4], in,
 * tol = in,
 * comm = in,
 * verbose = in
 */
int parRSB_findConnectivity(long long *vtx, double *coord, int nel, int nDim,
                            long long *periodicInfo, int nPeriodicFaces,
                            double tol, MPI_Comm comm, int verbose);

/* Misc */

int parrsb_read_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                     char *name, MPI_Comm comm, int read);

int parrsb_distribute_elements(unsigned int *nelt, long long **vl,
                               double **coord, int *part, int nv,
                               MPI_Comm comm);

void parrsb_part_stat(long long *vtx, int nel, int nv, MPI_Comm comm);

struct csr_mat_ *parrsb_numbering(unsigned int *nelt, unsigned int *nlevels,
                                  unsigned int **level_off, MPI_Comm **comms,
                                  long long *vl, double *coord, int nv,
                                  MPI_Comm comm);

void parrsb_ilu0(unsigned int nlevels, unsigned int *level_off, MPI_Comm *comms,
                 struct csr_mat_ *M, MPI_Comm world);

void parrsb_lu_solve(double *x, struct csr_mat_ *M, double *b,
                     unsigned int nlevels, unsigned int *level_off,
                     MPI_Comm *comms, MPI_Comm world);
#endif
