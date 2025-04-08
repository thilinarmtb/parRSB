#ifndef _PARRSB_H_
#define _PARRSB_H_

#include "gslib.h"

#if !defined(GSLIB_USE_MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GSLIB_USE_GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * parRSB options:
 */
typedef struct parrsb_options *parrsb_options;

int parrsb_options_get_default(parrsb_options *options);

int parrsb_options_copy(parrsb_options *dest, const parrsb_options src);

void parrsb_options_print(const parrsb_options options);

int parrsb_options_free(parrsb_options *options);

/*
 * parRSB partitioning algorithms.
 */
int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int *const tag,
                     const int nel, const int nv, const parrsb_options options,
                     MPI_Comm comm);

int parrsb_part_graph(int *part, size_t num_nodes, long long *nodes,
                      size_t *offsets, long long *neighbors,
                      const parrsb_options options, MPI_Comm comm);

/*
 * Connectivity calculation algorithms.
 */
int parrsb_conn_mesh(long long *vtx, double *coord, uint nel, unsigned nDim,
                     long long *periodicInfo, int nPeriodicFaces, double tol,
                     MPI_Comm comm);

#define fparrsb_conn_mesh                                                      \
  GS_FORTRAN_UNPREFIXED(fparrsb_conn_mesh, FPARRSB_CONN_MESH)
void fparrsb_conn_mesh(long long *vtx, double *coord, int *nel, int *nDim,
                       long long *periodicInfo, int *nPeriodicFaces,
                       double *tol, MPI_Fint *fcomm, int *err);

/*
 * I/O routines.
 */
int parrsb_read_mesh(unsigned *nel, unsigned *nv, long long **vl,
                     double **coord, unsigned *nbcs, long long **bcs,
                     char *name, MPI_Comm comm, int read);

int parrsb_dump_con(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

int parrsb_dump_map(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

/*
 * Auxiliary functions.
 */
typedef struct {
  char *mesh;  // Mesh name, required.
  double tol;  // gencon tolerance, default: 0.2
  int test;    // run tests, default: 0
  int dump;    // dump the connectivity or map file, default: 1
  int nactive; // # of active MPI ranks, default: INT_MAX
  int verbose; // Verbosity, default: 0
} parrsb_cmd_line_opts;

parrsb_cmd_line_opts *parrsb_parse_cmd_opts(int argc, char *argv[]);

void parrsb_cmd_opts_free(parrsb_cmd_line_opts *opts);

int parrsb_dist_mesh(unsigned *nelt, long long **vl, double **coord, int *part,
                     unsigned nv, MPI_Comm comm);

int parrsb_setup_mesh(unsigned *nelt, unsigned *nv, long long **vl,
                      double **coord, parrsb_cmd_line_opts *opts,
                      MPI_Comm comm);

void parrsb_print_part_stat(long long *vtx, unsigned nelt, unsigned nv,
                            MPI_Comm comm);

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm comm);

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm);
#define parrsb_check_error(err, comm)                                          \
  parrsb_check_error_(err, __FILE__, __LINE__, comm)

#ifdef __cplusplus
}
#endif

#endif
