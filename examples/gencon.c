/*
 * Generate connectivity from Nek5000 mesh (.re2) file.
 */
#include "parrsb_example.h"

static void check_connectivity_aux(unsigned nelt, unsigned nv,
                                   const long long *const vl1,
                                   const long long *const vl2, struct comm *c,
                                   buffer *bfr) {
  unsigned size = nelt * nv;
  long long *minp = tcalloc(long long, size);
  long long *maxp = tcalloc(long long, size);

  struct gs_data *gsh = gs_setup(vl1, size, c, 0, gs_pairwise, 0);
  for (unsigned i = 0; i < size; i++) minp[i] = maxp[i] = vl2[i];
  gs(minp, gs_long, gs_min, 0, gsh, bfr);
  gs(maxp, gs_long, gs_max, 0, gsh, bfr);
  gs_free(gsh);

  unsigned i = 0;
  for (; i < size; i++)
    if (minp[i] != maxp[i]) break;

  free(minp), free(maxp);
  parrsb_check_error(i != size, c->c);
}

static void check_connectivity(const long long *const vlp, char *name,
                               MPI_Comm comm) {
  unsigned nelt, nv;
  long long *vls = NULL;
  parrsb_check_error(parrsb_read_conn(&nelt, &nv, &vls, name, comm), comm);

  struct comm c;
  comm_init(&c, comm);

  buffer bfr;
  buffer_init(&bfr, 1024);

  check_connectivity_aux(nelt, nv, vlp, vls, &c, &bfr);
  check_connectivity_aux(nelt, nv, vls, vlp, &c, &bfr);

  unsigned i = 0, size = nelt * nv;
  if (c.np == 1) {
    for (; i < size; i++)
      if (vls[i] != vlp[i]) break;
    parrsb_check_error(i != size, comm);
  }

  comm_free(&c);
  buffer_free(&bfr);
  free(vls);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  parrsb_example_opts in = parrsb_example_opts_parse(argc, argv, comm);

  // Read the geometry from the .re2 file.
  unsigned nelt, nbcs, nv;
  double *coord = 0;
  long long *bcs = 0;
  parrsb_check_error(
      parrsb_read_mesh(&nelt, &nv, &coord, &nbcs, &bcs, in->mesh, comm), comm);

  // Find connectivity.
  unsigned ndim = (nv == 8 ? 3 : 2);
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  parrsb_check_error(
      parrsb_conn_mesh(vl, coord, nelt, ndim, bcs, nbcs, in->tol, comm), comm);

  // Turns on testing if test is on
  if (in->test) check_connectivity(vl, in->mesh, comm);

  // Free resources
  free(vl), free(coord), free(bcs);
  parrsb_example_opts_free(&in);
  MPI_Finalize();

  return 0;
}
