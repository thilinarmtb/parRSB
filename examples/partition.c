/*
 * Parition a mesh with RSB using Nek5000's vertex connectivity (co2) file.
 * Geometric partitioning methods which is mesh geometry in .re2 file can
 * be used either as an alternative to RSB or to increase RSB convergence.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <parRSB.h>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc != 3) {
    if (rank == 0)
      printf("Usage: %s <#nread> <mesh file>\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  int n_read = atoi(argv[1]);
  char *mesh = argv[2];

  int color = 0;
  if (rank < n_read)
    color = 3; /* Read both .re2 (1) + and .co2 (2) */
  MPI_Comm comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &comm);

  /* Read Nek5000 mesh and print mesh statistics */
  unsigned int nelt = 0;
  int nv = 0;
  long long *vl = NULL;
  double *coord = NULL;
  int ierr = read_nek_mesh(&nelt, &nv, &vl, &coord, mesh, comm, color);

  /* Find the partition, distribute elements and print mesh statistics */
  int *part = NULL;
  if (ierr == 0) {
    parrsb_part_stat(vl, nelt, nv, MPI_COMM_WORLD);

    parRSB_options options = parrsb_default_options;
    part = (int *)calloc(nelt, sizeof(int));
    ierr |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, &options,
                            MPI_COMM_WORLD);
  }

  if (ierr == 0) {
    /* Compress coords to find centroids */
    parrsb_distribute_elements(nelt, nv, part, vl, coord, MPI_COMM_WORLD);
    parrsb_part_stat(vl, nelt, nv, MPI_COMM_WORLD);
  }

  if (part != NULL)
    free(part);
  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}
