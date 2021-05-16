#include <genmap-iterative.h>
#include <genmap.h>

#include <parRSB.h>

static int project_ilu(genmap_handle h) {
  struct comm *lc = h->local;
  struct comm *gc = h->global;

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;
  assert(np == 1);

  uint lelt = genmap_get_nel(h);

  genmap_vector init;
  genmap_vector_create(&init, lelt);

  uint i;
  for (i = 0; i < lelt; i++)
    init->data[i] = genmap_get_local_start_index(h) + i + 1;
  genmap_vector_ortho_one(gc, init, lelt);

  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  struct precond *d = precond_setup(2, h, gc);

  int ppfi = project(h, gc, d, init, 100, y);

  precond_free(d);

  genmap_destroy_vector(y);

  genmap_destroy_vector(init);

  return ppfi;
}

static int partition(unsigned int *nelt_, int *nv_, long long **vl,
                     double **coord, char *mesh, double tol, MPI_Comm comm) {
  /* Read the geometry from the .re2 file */
  unsigned int nelt = 0;
  int nv = 0;
  int ierr = read_nek_mesh(&nelt, &nv, NULL, coord, mesh, comm, 1);

  /* Calculate connectivity */
  if (ierr == 0) {
    *vl = (long long *)calloc(nelt * nv, sizeof(long long));
    int ndim = (nv == 8) ? 3 : 2;
    ierr |=
        parRSB_findConnectivity(*vl, *coord, nelt, ndim, NULL, 0, tol, comm, 0);
  }

  /* Print mesh statistics and find a potentail partition */
  int *part = NULL;
  if (ierr == 0) {
    parrsb_part_stat(*vl, nelt, nv, comm);

    parRSB_options options = parrsb_default_options;
    part = (int *)calloc(nelt, sizeof(int));
    ierr |= parRSB_partMesh(part, NULL, *vl, *coord, nelt, nv, &options, comm);
  }

  /* Distribute elements to the partition found and print mesh statistics */
  if (ierr == 0) {
    parrsb_distribute_elements(nelt, nv, part, vl, coord, comm);
    parrsb_part_stat(*vl, nelt, nv, comm);
  }

  if (part != NULL)
    free(part);

  *nelt_ = nelt;
  *nv_ = nv;

  return ierr;
}

static genmap_handle setup_laplacian(unsigned int nelt, int nv, long long *vl,
                                     double *coord, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  buffer buf;
  buffer_init(&buf, 1024);

  struct array elems;
  genmap_load_balance(&elems, nelt, nv, coord, vl, &cr, &buf);

  genmap_handle h;
  parRSB_options options = parrsb_default_options;
  genmap_init(&h, comm, &options);

  genmap_set_elements(h, &elems);
  genmap_set_nvertices(h, nv);

  struct comm *gc = genmap_global_comm(h);
  genmap_comm_scan(h, gc);
  genmap_number_faces_and_edges(h, gc);

  genmap_laplacian_init(h, gc);

  buffer_free(&buf);
  crystal_free(&cr);
  comm_free(&c);
  array_free(&elems);

  return h;
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
  partition(&nelt, &nv, &vl, &coord, mesh, tol, MPI_COMM_WORLD);

  genmap_handle h = setup_laplacian(nelt, nv, vl, coord, MPI_COMM_WORLD);

  project_ilu(h);

  genmap_finalize(h);

  if (vl != NULL)
    free(vl);
  if (coord != NULL)
    free(coord);

  MPI_Finalize();

  return 0;
}
