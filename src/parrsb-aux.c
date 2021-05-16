#include <mpi.h>
#include <stdio.h>

#include <gencon.h>
#include <genmap.h>

#define MAXNV 8 /* maximum number of vertices per element */

typedef struct {
  int proc;
  long long vtx[MAXNV];
} elm_data;

int read_nek_mesh(unsigned int *nel, int *nv, long long **vl, double **coord,
                  char *name, MPI_Comm comm, int read) {
  char geom_name[BUFSIZ];
  strncpy(geom_name, name, BUFSIZ);
  strncat(geom_name, ".re2", 5);

  char conn_name[BUFSIZ];
  strncpy(conn_name, name, BUFSIZ);
  strncat(conn_name, ".co2", 5);

  int ndim_;

  struct comm c;
  comm_init(&c, comm);

  Mesh mesh;

  /* Read mesh data */
  if (read & 1) {
    read_geometry(&mesh, geom_name, &c);
    get_vertex_coordinates(coord, mesh);
  }

  /* FIXME: Calculate connectivity instead of the following */
  if (read & 2) {
    assert(read & 1);
    read_connectivity(mesh, conn_name, &c);
    get_vertex_ids(vl, mesh);
  }

  ndim_ = get_mesh_dim(mesh);
  *nel = (unsigned int)get_mesh_nel(mesh);

  mesh_free(mesh);
  comm_free(&c);

  MPI_Bcast(&ndim_, 1, MPI_INT, 0, MPI_COMM_WORLD);
  *nv = (ndim_ == 3) ? 8 : 4;

  return 0;
}

int redistribute_elements(unsigned int nelt, int nv, int *part, long long *vl,
                          MPI_Comm comm) {
  struct array elements;
  array_init(elm_data, &elements, nelt);

  elm_data *data;
  int e, n;
  for (data = elements.ptr, e = 0; e < nelt; ++e) {
    data[e].proc = part[e];
    for (n = 0; n < nv; ++n)
      data[e].vtx[n] = vl[e * nv + n];
  }
  elements.n = nelt;

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  sarray_transfer(elm_data, &elements, proc, 0, &cr);

  nelt = elements.n;
  vl = (long long *)realloc(vl, nv * nelt * sizeof(long long));
  for (data = elements.ptr, e = 0; e < nelt; ++e) {
    for (n = 0; n < nv; ++n)
      vl[e * nv + n] = data[e].vtx[n];
  }

  crystal_free(&cr);
  comm_free(&c);
  array_free(&elements);

  return 0;
}

void parrsb_part_stat(long long *vtx, int nel, int nv, MPI_Comm ce) {
  int i, j;

  struct comm comm;
  int np, id;

  int Nmsg;
  int *Ncomm;

  int nelMin, nelMax;
  int ncMin, ncMax, ncSum;
  int nsMin, nsMax, nsSum;
  int nssMin, nssMax, nssSum;

  struct gs_data *gsh;
  int b;

  int numPoints;
  long long *data;

  comm_init(&comm, ce);
  np = comm.np;
  id = comm.id;

  if (np == 1)
    return;

  numPoints = nel * nv;
  data = (long long *)malloc((numPoints + 1) * sizeof(long long));
  for (i = 0; i < numPoints; i++)
    data[i] = vtx[i];

  gsh = gs_setup(data, numPoints, &comm, 0, gs_pairwise, 0);

  pw_data_nmsg(gsh, &Nmsg);
  Ncomm = (int *)malloc((Nmsg + 1) * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (i = 1; i < Nmsg; ++i) {
    nsMax = Ncomm[i] > Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsMin = Ncomm[i] < Ncomm[i - 1] ? Ncomm[i] : Ncomm[i - 1];
    nsSum += Ncomm[i];
  }
  comm_allreduce(&comm, gs_int, gs_max, &nsMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nsMin, 1, &b);

  nssMin = nsSum;
  nssMax = nsSum;
  nssSum = nsSum;
  comm_allreduce(&comm, gs_int, gs_max, &nssMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nssMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nssSum, 1, &b);

  if (Nmsg)
    nsSum = nsSum / Nmsg;
  else
    nsSum = 0;
  comm_allreduce(&comm, gs_int, gs_add, &nsSum, 1, &b);

  nelMax = nel;
  nelMin = nel;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);

  if (id == 0) {
    printf(" Max neighbors: %d | Min neighbors: %d | Avg neighbors: %lf\n",
           ncMax, ncMin, (double)ncSum / np);
    printf(" Max nvolume: %d | Min nvolume: %d | Avg nvolume: %lf\n", nsMax,
           nsMin, (double)nsSum / np);
    printf(" Max volume: %d | Min volume: %d | Avg volume: %lf\n", nssMax,
           nssMin, (double)nssSum / np);
    printf(" Max elements: %d | Min elements: %d\n", nelMax, nelMin);
    fflush(stdout);
  }

  free(Ncomm);
  comm_free(&comm);
}

#undef MAXNV
