#include <genmap-impl.h>
#include <genmap-iterative.h>
#include <genmap-precond.h>

int genmap_ilu_project(genmap_handle h) {
  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_number_faces_and_edges(h, gc);
  genmap_comm_scan(h, gc);

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;
  assert(np == 1);

  metric_tic(lc, RCB);
  rcb(lc, h->elements, ndim, &h->buf);
  metric_toc(lc, RCB);

  uint lelt = genmap_get_nel(h);
  genmap_vector init;
  genmap_vector_create(&init, lelt);

  struct rsb_element *elements = genmap_get_elements(h);
  uint i;
  for (i = 0; i < lelt; i++)
    init->data[i] = genmap_get_local_start_index(h) + i + 1;

  genmap_vector_ortho_one(gc, init, lelt);

  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  genmap_laplacian_init(h, gc);

  metric_tic(gc, PRECONDSETUP);
  struct precond *d = precond_setup(2, h, gc);
  metric_toc(gc, PRECONDSETUP);

  metric_tic(gc, PROJECT);
  int ppfi = project(h, gc, d, init, 100, y);
  metric_toc(gc, PROJECT);
  metric_acc(NPROJECT, ppfi);

  precond_free(d);

  genmap_destroy_vector(y);
  genmap_destroy_vector(init);

  return 0;
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

  /* Read the geometry from the .re2 file */
  unsigned int nelt;
  int nv;
  double *coord = NULL;
  int ierr = read_nek_mesh(&nelt, &nv, NULL, &coord, mesh, MPI_COMM_WORLD, 1);

  /* Calculate connectivity */
  long long *vl = NULL;
  int *part = NULL;
  if (ierr == 0) {
    vl = (long long *)calloc(nelt * nv, sizeof(long long));
    int ndim = nv == 8 ? 3 : 2;
    ierr |= parRSB_findConnectivity(vl, coord, nelt, ndim, NULL, 0, tol,
                                    MPI_COMM_WORLD, 0);
  }

  /* Print mesh statistics and find a potentail partition */
  if (ierr == 0) {
    parrsb_part_stat(vl, nelt, nv, MPI_COMM_WORLD);

    parRSB_options options = parrsb_default_options;
    part = (int *)calloc(nelt, sizeof(int));
    ierr |= parRSB_partMesh(part, NULL, vl, coord, nelt, nv, &options,
                            MPI_COMM_WORLD);
  }

  /* Distribute elements to the partition found and print mesh statistics */
  if (ierr == 0) {
    redistribute_elements(nelt, nv, part, vl, MPI_COMM_WORLD);
    parrsb_part_stat(vl, nelt, nv, MPI_COMM_WORLD);
  }

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);
  if (part != NULL)
    free(part);

  MPI_Finalize();

  return 0;
}
