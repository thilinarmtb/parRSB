#include <genmap-iterative.h>
#include <genmap.h>

#include <parRSB.h>

static int read_mesh_and_con(unsigned int *nelt_, int *nv_, long long **vl,
                             double **coord, char *mesh, double tol,
                             MPI_Comm comm) {
  /* Read the geometry from the .re2 file */
  unsigned int nelt = 0;
  int nv = 0;
  int ierr = parrsb_read_mesh(&nelt, &nv, NULL, coord, mesh, comm, 1);

  /* Calculate connectivity */
  if (ierr == 0) {
    *vl = (long long *)calloc(nelt * nv, sizeof(long long));
    int ndim = (nv == 8) ? 3 : 2;
    ierr |=
        parRSB_findConnectivity(*vl, *coord, nelt, ndim, NULL, 0, tol, comm, 0);
  }

  *nelt_ = nelt;
  *nv_ = nv;

  return ierr;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 2) {
    if (rank == 0)
      printf("Usage: ./%s <mesh file> [tol]\n", argv[0]);
    MPI_Finalize();
    return 1;
  }

  char *mesh = argv[1];
  double tol = (argc > 2) ? atof(argv[2]) : 0.2;

  /* Read in the mesh */
  unsigned int nelt = 0;
  int nv = 0;
  long long *vl = NULL;
  double *coord = NULL;
  read_mesh_and_con(&nelt, &nv, &vl, &coord, mesh, tol, MPI_COMM_WORLD);

  /* Generate CSR matrix with RSB ordering */
  unsigned int nlevels = 0;
  unsigned int *level_off = NULL;
  MPI_Comm *comms = NULL;
  struct csr_mat_ *M = parrsb_numbering(&nelt, &nlevels, &level_off, &comms, vl,
                                        coord, nv, MPI_COMM_WORLD);

  double *x = tcalloc(double, M->rn);
  double *b = tcalloc(double, M->rn);

  int i;
  for (i = 0; i < M->rn; i++)
    b[i] = i + 1;

  parrsb_lu_solve(x, M, b, nlevels, level_off, comms, MPI_COMM_WORLD);

  if (x != NULL)
    free(x);
  if (b != NULL)
    free(b);

  if (M != NULL)
    csr_mat_free(M);
  if (comms != NULL) {
    for (i = 0; i < nlevels; i++)
      MPI_Comm_free(&comms[i]);
    free(comms);
  }
  if (level_off != NULL)
    free(level_off);

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Finalize();

  return 0;
}
