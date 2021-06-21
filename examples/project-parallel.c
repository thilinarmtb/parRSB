#include <genmap-csr-mat.h>
#include <genmap-iterative.h>
#include <genmap.h>

#include <parRSB.h>

static int project_ilu(genmap_handle h, unsigned int lelt) {
  struct comm *gc = h->global;
  struct precond *d = precond_setup(3, h, gc);

  genmap_vector init;
  genmap_vector_create(&init, lelt);

  slong out[2][1], buf[2][1];
  slong in = lelt;
  comm_scan(out, gc, gs_long, gs_add, &in, 1, buf);

  uint i;
  for (i = 0; i < lelt; i++)
    //init->data[i] = out[0][0] + i + 1;
    init->data[i] = h->M->row_id[i];

  genmap_vector_ortho_one(gc, init, out[1][0]);

  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  project(h, gc, d, init, 100, y);

  precond_free(d);
  genmap_destroy_vector(y);
  genmap_destroy_vector(init);

  return 0;
}

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

  long long *vl = NULL;
  double *coord = NULL;
  unsigned int nelt = 0;
  int nv = 0;
  read_mesh_and_con(&nelt, &nv, &vl, &coord, mesh, tol, MPI_COMM_WORLD);

  genmap_handle h =
      parrsb_numbering_w_handle(&nelt, vl, coord, nv, MPI_COMM_WORLD);

  csr_mat_dump("Reordered.dump", h->M, MPI_COMM_WORLD);

  parRSB_options options = parrsb_default_options;
  options.rsb_laplacian_implementation = 2;
  h->options = &options;

  project_ilu(h, nelt);
  genmap_finalize(h);

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Finalize();

  return 0;
}
