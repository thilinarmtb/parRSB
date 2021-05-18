#include <genmap-iterative.h>
#include <genmap.h>

#include <parRSB.h>

#define write_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    T v = val;                                                                 \
    memcpy(dest, &(val), sizeof(T) * nunits);                                  \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

int dump_csr_mat(const char *fname, struct csr_mat_ *M, MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_File file;
  int err = MPI_File_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);
  if (err != 0)
    return err;

  long nelt = M->rn;
  long nnz = M->row_off[nelt];

  long nnzg;
  MPI_Allreduce(&nnz, &nnzg, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  /* Write i, j, v, rank */
  uint write_size = (2 * sizeof(long) + sizeof(double) + sizeof(int)) * nnz;
  if (rank == 0)
    write_size += sizeof(long); /* for nelgt */

  char *pbuf, *pbuf0;
  pbuf = pbuf0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0)
    write_T(pbuf0, nnzg, long, 1);

  uint i, j;
  for (i = 0; i < nelt; i++) {
    slong row = M->row_id[i];
    for (j = M->row_off[i]; j < M->row_off[i + 1]; j++) {
      write_T(pbuf0, row, long, 1);
      write_T(pbuf0, M->col[j], long, 1);
      write_T(pbuf0, M->v[j], double, 1);
      write_T(pbuf0, rank, int, 1);
    }
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, pbuf, write_size, MPI_BYTE, &st);
  if (err != 0 && rank == 0)
    fprintf(stderr, "%s:%d Error opening file %s for writing.\n", __FILE__,
            __LINE__, fname);

  err += MPI_File_close(&file);
  MPI_Barrier(comm);

  free(pbuf);

  return err;
}

#undef write_T

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
  unsigned int *map = NULL;
  long long *gid = NULL;
  unsigned int nlevels = 0;
  unsigned int *level_off = NULL;
  struct csr_mat_ *R = parrsb_numbering(&nelt, &nlevels, &level_off, vl, coord,
                                        nv, MPI_COMM_WORLD);

  /* Dump the matrix */
  dump_csr_mat("Reordered.dump", R, MPI_COMM_WORLD);

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (gid != NULL)
    free(gid);
  if (map != NULL)
    free(map);
  if (level_off != NULL)
    free(level_off);
  if (R != NULL)
    csr_mat_free(R);

  MPI_Finalize();

  return 0;
}
