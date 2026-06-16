#ifndef _PARRSB_H_
#define _PARRSB_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <mpi.h>

/*
 * parRSB options:
 */
typedef struct parrsb_options *parrsb_options_t;

int parrsb_options_get_default(parrsb_options_t *options);

int parrsb_options_update(parrsb_options_t options);

int parrsb_options_set_partitioner(parrsb_options_t options, int partitioner);

void parrsb_options_print(const parrsb_options_t options);

int parrsb_options_free(parrsb_options_t *options);

/*
 * parRSB partitioning algorithms.
 */
int parrsb_part_mesh(int *part, const long long *const vtx,
                     const double *const xyz, const int nel, const int nv,
                     const parrsb_options_t options, MPI_Comm comm);

int parrsb_part_graph(int *part, unsigned num_nodes,
                      const long long *const nodes,
                      const unsigned *const offsets,
                      const long long *const neighbors,
                      const parrsb_options_t options, const MPI_Comm comm);

/*
 * Connectivity calculation algorithms for SEM meshes.
 */
int parrsb_conn_mesh(long long *vtx, const double *const coord, unsigned nel,
                     unsigned ndim, const long long *const pinfo, int npinfo,
                     double tol, MPI_Comm comm);

/*
 * I/O routines.
 */
int parrsb_read_mesh(unsigned *nel, unsigned *nv, double **coord,
                     unsigned *nbcs, long long **bcs, char *name,
                     MPI_Comm comm);

int parrsb_read_conn(unsigned *nel, unsigned *nv, long long **vtx, char *name,
                     MPI_Comm comm);

int parrsb_dump_con(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

int parrsb_dump_map(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm);

/*
 * Auxiliary functions.
 */
int parrsb_dist_mesh(unsigned *nelt, long long **vl, double **coord, int *part,
                     unsigned nv, MPI_Comm comm);

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif
