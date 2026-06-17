/*
 * Generate partitions from Nek5000's mesh (.re2) file.
 */
#include "parrsb_example.h"

static void partition(long long *vl, double *coord, unsigned nelt, unsigned nv,
                      parrsb_example_opts in, MPI_Comm world, MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(world, &rank);
  MPI_Comm_size(world, &size);

  // Print and store pre-partition statistics.
  if (in->verbose > 0) {
    if (rank == 0) printf("Partition statistics before RSB:\n");
    fflush(stdout);
    parrsb_print_part_stat(vl, nelt, nv, comm);
  }

  int nss[6];
  parrsb_part_stat(NULL, NULL, nss, NULL, vl, nelt, nv, comm);

  // Partition the mesh.
  int *part = (int *)calloc(nelt, sizeof(int));

  parrsb_options_t options;
  parrsb_options_get_default(&options);
  int err = parrsb_part_mesh(part, vl, coord, nelt, nv, options, comm);
  parrsb_options_free(&options);
  parrsb_check_error(err, comm);

  // Redistribute data based on identified partitions.
  parrsb_check_error(parrsb_dist_mesh(&nelt, &vl, &coord, part, nv, comm),
                     comm);

  // Store post-partition statistics.
  if (in->verbose > 0) {
    if (rank == 0) printf("Partition statistics after RSB:\n");
    fflush(stdout);
    parrsb_print_part_stat(vl, nelt, nv, comm);
  }
  parrsb_part_stat(NULL, NULL, &nss[3], NULL, vl, nelt, nv, comm);

  if (in->test && in->nactive > 1) parrsb_check_error((nss[2] < nss[5]), comm);

  free(part);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm world = MPI_COMM_WORLD;

  parrsb_example_opts in = parrsb_example_opts_parse(argc, argv, world);

  int rank, size;
  MPI_Comm_rank(world, &rank);
  MPI_Comm_size(world, &size);

  if (in->nactive > size) in->nactive = size;

  MPI_Comm comm;
  MPI_Comm_split(world, rank < in->nactive, rank, &comm);
  if (rank >= in->nactive) goto finalize;

  // Read the geometry from the .re2 file.
  unsigned nelt, nbcs, nv, err;
  double *coord = 0;
  long long *bcs = 0;
  err = parrsb_read_mesh(&nelt, &nv, &coord, &nbcs, &bcs, in->mesh, comm);
  parrsb_check_error(err, comm);

  // Find connectivity.
  long long *vl = (long long *)calloc(nelt * nv, sizeof(long long));
  unsigned ndim = (nv == 8 ? 3 : 2);
  err = parrsb_conn_mesh(vl, coord, nelt, ndim, bcs, nbcs, in->tol, comm);
  parrsb_check_error(err, comm);

  partition(vl, coord, nelt, nv, in, world, comm);

  free(vl), free(coord), free(bcs);

finalize:
  parrsb_example_opts_free(&in);
  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}
