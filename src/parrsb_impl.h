#ifndef _PARRSB_IMPL_H_
#define _PARRSB_IMPL_H_

#define _POSIX_C_SOURCE 200809L

#include <float.h>
#include <stdlib.h>

#include "gslib.h"
#include "parrsb.h"

typedef double scalar;
#define SCALAR_MAX DBL_MAX
#define SCALAR_MIN DBL_MIN
#define SCALAR_TOL 1e-12
#define gs_scalar gs_double

#if !defined(GSLIB_USE_MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GSLIB_USE_GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

/*
 * Set the visibility of a symbol.
 */
#define VISIBILITY_(mode) __attribute__((visibility(#mode)))
#define VISIBILITY(mode) VISIBILITY_(mode)

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
#if !defined(PARRSB_VISIBILITY)
#define PARRSB_VISIBILITY hidden
#endif

#if defined(__cplusplus)
#define INTERN extern "C" VISIBILITY(PARRSB_VISIBILITY)
#else
#define INTERN extern VISIBILITY(PARRSB_VISIBILITY)
#endif

/*
 * Commonly used macros.
 */
// Macro for finding min of two numbers.
#define MIN(a, b) ((a) < (b) ? (a) : (b))
// Macro for finding max of two numbers.
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*
 * Compile time constants.
 */
// Maximum dimension of the mesh.
#define MAXDIM 3
// Maximum number of vertices per element in the mesh.
#define MAXNV 8
// Maximum depth of the call stack for arena.
#define MAXDEPTH 32

/*
 * parRSB configuration options.
 */
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
 * Element info struct keeps track of align, size and number of
 * vertices and dimensions (latter two in case of a mesh) of the
 * input.
 */
struct element_info {
  size_t size, align;
  int nd, nv;
};

typedef struct element_info *element_info;

typedef enum { MESH = 0, GRAPH = 1 } element_type_t;

void element_info_init(element_info *ei, element_type_t type);

void element_info_free(element_info *ei);

element_type_t element_info_type(element_info ei);

/*
 * Base element for mesh and graph partitioning. All of the structures used for
 * partitioning should append other fields to this structure.
 */
struct base_element {
  uint proc, origin;
  ulong globalId;
  scalar fiedler;
};

typedef struct base_element *base_element;

/*
 * General graph partitioning.
 */
struct graph_element {
  uint proc, origin;
  ulong globalId;
  scalar fiedler;
  ulong u, v;
};

typedef struct graph_element *graph_element;

int graph_load_balance(struct array *nlist, uint nn, const long long *nodes,
                       const unsigned *offsets, const long long *neighbors,
                       struct crystal *cr, buffer *bfr);

void graph_restore(int *part, struct array *nlist, struct crystal *cr,
                   buffer *bfr);

/*
 * RCB / RIB. `struct rcb_element` is used for mesh partitioning with RCB and
 * RIB.
 */
struct rcb_element {
  uint proc, origin;
  ulong globalId;
  scalar fiedler;
  scalar coord[MAXDIM];
};

typedef struct rcb_element *rcb_element;

INTERN int rcb(struct array *elements, const element_info ei,
               const struct comm *c, buffer *bfr);

INTERN int rib(struct array *elements, const element_info ei,
               const struct comm *c, buffer *bfr);

/*
 * Laplacian.
 */
typedef struct laplacian *laplacian;

int laplacian_init(laplacian *l, struct array *elements, const element_info ei,
                   const struct comm *c);
uint laplacian_size(laplacian l);
int laplacian_op(scalar *v, laplacian l, scalar *u);
int laplacian_free(laplacian *l);

/*
 * RSB. `struct rsb_element` = `struct rcb_element` + `vertices`. Order is
 * important. Used for mesh partitioning with RSB.
 */
struct rsb_element {
  uint proc, origin;
  ulong globalId;
  scalar fiedler;
  scalar coord[MAXDIM];
  slong vertices[MAXNV];
};

typedef struct rsb_element *rsb_element;

INTERN void rsb(struct array *elements, const element_info ei,
                const parrsb_options_t options, const struct comm *comms,
                buffer *bfr);

int lanczos(scalar *f, laplacian l, scalar *vi, const struct comm *c,
            const parrsb_options_t opts);

INTERN int fiedler(scalar *fiedler, laplacian l, const parrsb_options_t opts,
                   const struct comm *c);

/*
 * Find the number of disconnected components in the graph.
 */
INTERN slong get_components(sint *component, const struct array *elems,
                            const element_info ei, const struct comm *c,
                            buffer *buf);

/*
 * Helper routines.
 */
INTERN void parrsb_barrier(struct comm *c);

INTERN void parrsb_info(const struct comm *c, int verbose, const char *fmt,
                        ...);

INTERN void parrsb_error(int cond, const struct comm *c, int verbose,
                         const char *fmt, ...);

INTERN void parrsb_print_stack(void);

INTERN int log2ll(long long n);

INTERN int nv_to_ndim(int nv);

INTERN int matrix_inverse(int N, scalar *A);

/*
 * General eigenvalue related routines.
 */
INTERN int power_iter(double *y, uint N, double *A);

INTERN int inv_power_iter(double *y, uint N, double *A);

INTERN int tqli(scalar *evectors, scalar *evalues, sint n, scalar *diagonal,
                scalar *upper, int id);
#endif
