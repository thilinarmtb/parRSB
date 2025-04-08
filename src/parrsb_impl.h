#ifndef _PARRSB_IMPL_H_
#define _PARRSB_IMPL_H_

#define _POSIX_C_SOURCE 200809L

#include "parrsb.h"
#include "types.h"

#include <float.h>
#include <stdlib.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

struct parrsb_options {
  // General options
  int partitioner; // Partition algo: 0 - RSB, 1 - RCB, 2 - RIB (Default: 0)
  int tagged;      // Tagged partitioning: 0 - No, 1 - Yes (Default: 0)
  int levels;      // Number of levels: 1, or 2 (Default: 2)
  int find_disconnected_comps; // Find number of components: 0 - No, 1 - Yes
                               // (Default: 1)
  int repair; // Repair disconnected components: 0 - No, 1 - Yes (Default: 0)
  int verbose_level; // Verbose level: 0, 1, 2, .. etc (Default: 1)
  int profile_level; // Profile level: 0, 1, 2, .. etc (Default: 0)
  // RSB common (Lanczos and MG) options
  int rsb_algo; // RSB algo: 0 - Lanczos, 1 - MG (Default: 0)
  int rsb_pre;  // RSB pre-partition : 0 - None, 1 - RCB , 2 - RIB (Default: 1)
  int rsb_max_iter;   // Max iterations in Lanczos / MG (Default: 50)
  int rsb_max_passes; // Max Lanczos restarts / Inverse iterations (Default: 50)
  double rsb_tol;     // Tolerance for Lanczos or RQI (Default: 1e-5)
  int rsb_dump_stats; // Dump partition statistics to a text file.
  // RSB MG specific options
  int rsb_mg_grammian; // MG Grammian: 0 or 1 (Default: 0)
  int rsb_mg_factor;   // MG Coarsening factor (Default: 2, should be > 1)
};

/*
 * Set the visibility of a symbol.
 */
#define VISIBILITY(mode) __attribute__((visibility(#mode)))

/*
 * Declare a symbol as external.
 */
#if defined(__cplusplus)
#define EXTERN extern "C" VISIBILITY(default)
#else
#define EXTERN extern VISIBILITY(default)
#endif

/*
 * Declare a symbol as internal.
 */
#if defined(__cplusplus)
#define INTERN extern "C" VISIBILITY(hidden)
#else
#define INTERN extern VISIBILITY(hidden)
#endif

#define MAXDIM 3 // Maximum dimension of the mesh.
#define MAXNV 8  // Maximum number of vertices per element.

/*
 * RCB / RIB. `struct rcb_element` is used for RCB and RIB partitioning.
 */
struct rcb_element {
  uint proc, origin;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
};

INTERN int rcb(struct array *elements, size_t unit_size, int ndim,
               const struct comm *c, buffer *bfr);
INTERN int rib(struct array *elements, size_t unit_size, int ndim,
               const struct comm *c, buffer *bfr);

/*
 * RSB. `struct rsb_element` = `struct rcb_element` + `vertices`. Order is
 * important.
 */
struct rsb_element {
  uint proc, origin;
  ulong globalId;
  scalar coord[MAXDIM], fiedler;
  slong vertices[MAXNV];
};

INTERN void rsb(struct array *elements, int nv, const parrsb_options options,
                const struct comm *comms, buffer *bfr);

INTERN int fiedler(struct array *elements, int nv, const parrsb_options options,
                   const struct comm *gsc, buffer *buf);

/*
 * Laplacian.
 */
#define GS 1
#define CSR 2

typedef struct laplacian *laplacian;
INTERN int laplacian_init(laplacian *l, struct rsb_element *elems, uint nel,
                          int nv, int type, const struct comm *c, buffer *bfr);
INTERN uint laplacian_get_size(laplacian wl);
INTERN int laplacian_op(scalar *v, laplacian l, scalar *u, buffer *bfr);
INTERN int laplacian_free(laplacian *l);

/*
 * Helper routines.
 */
INTERN void parrsb_barrier(struct comm *c);

INTERN void parrsb_print(const struct comm *c, int verbose, const char *fmt,
                         ...);

INTERN void parrsb_print_stack(void);

INTERN int log2ll(long long n);

INTERN int nv_to_ndim(int nv);

/*
 * Matrix inverse (local to MPI process).
 */
INTERN int matrix_inverse(int N, scalar *A);

/*
 * General eigenvalue related.
 */
INTERN int power_iter(double *y, uint N, double *A);

INTERN int inv_power_iter(double *y, uint N, double *A);

INTERN int tqli(scalar *evectors, scalar *evalues, sint n, scalar *diagonal,
                scalar *upper, int id);
#endif
