#include "parrsb_impl.h"

#include <math.h>

#define READ_T(coords, buf, T, nv)                                             \
  { memcpy((coords), buf, sizeof(T) * nv); }

#define WRITE_T(dest, val, T, nunits)                                          \
  {                                                                            \
    memcpy(dest, val, sizeof(T) * nunits);                                     \
    dest += sizeof(T) * nunits;                                                \
  }

#define WRITE_INT(dest, val)                                                   \
  { memcpy(dest, &(val), sizeof(int)); }

#define GC_RE2_HEADER_LEN 80
#define GC_CO2_HEADER_LEN 132
#define HEADER_LEN 132

#define check_mpi_call(call, comm)                                             \
  do {                                                                         \
    sint wrk;                                                                  \
    int err = ((call) != MPI_SUCCESS);                                         \
    comm_allreduce(c, gs_int, gs_add, &err, 1, &wrk);                          \
    if (err) return err;                                                       \
  } while (0)

static int re2_header(unsigned *nelt_, unsigned *nv_, ulong *nelgt_,
                      ulong *nelgv_, MPI_File file, struct comm *c) {
  char *buf = (char *)calloc(GC_RE2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  check_mpi_call(MPI_File_read_all(file, buf, GC_RE2_HEADER_LEN, MPI_BYTE, &st),
                 c);

  long long nelgt, nelgv;
  int ndim;
  char version[6];
  sscanf(buf, "%5s %lld %d %lld", version, &nelgt, &ndim, &nelgv);

  int rank = c->id, size = c->np;
  int nelt = nelgt / size, nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  float byte_test = 0;
  check_mpi_call(MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st), c);
  // check_call(fabs(byte_test - 6.543210) > 1e-7, "Byte test failed !", c);

  *nelt_ = nelt, *nv_ = (ndim == 2 ? 4 : 8), *nelgt_ = nelgt, *nelgv_ = nelgv;

  free(buf);
  return 0;
}

static int re2_coord(double **coord_, unsigned int nelt, int nv, MPI_File file,
                     struct comm *c) {
  unsigned ndim = (nv == 4) ? 2 : 3;
  size_t elem_size = nv * ndim * sizeof(double) + sizeof(double);
  size_t header_size = GC_RE2_HEADER_LEN + sizeof(float);

  // Calculate read size for element data on each MPI rank.
  uint rank = c->id;
  size_t read_size = nelt * elem_size + (rank == 0) * header_size;
  char *buf = (char *)calloc(read_size, sizeof(char));
  MPI_Status st;
  check_mpi_call(MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st), c);

  // Allocate coord array.
  size_t coord_size = nelt;
  coord_size = coord_size * nv * ndim;
  double *coord = *coord_ = tcalloc(double, coord_size);

  // Read elements for each rank.
  char *buf0 = buf + (rank == 0) * header_size;
  double x[8], y[8], z[8];
  if (ndim == 3) {
    for (uint i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(y, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(z, buf0, double, nv);
      buf0 += sizeof(double) * nv;

      for (int k = 0; k < nv; k++) {
        coord[i * nv * ndim + k * ndim + 0] = x[k];
        coord[i * nv * ndim + k * ndim + 1] = y[k];
        coord[i * nv * ndim + k * ndim + 2] = z[k];
      }
    }
  } else if (ndim == 2) {
    for (uint i = 0; i < nelt; i++) {
      // skip group id
      buf0 += sizeof(double);
      READ_T(x, buf0, double, nv);
      buf0 += sizeof(double) * nv;
      READ_T(y, buf0, double, nv);
      buf0 += sizeof(double) * nv;

      for (int k = 0; k < nv; k++) {
        coord[i * nv * ndim + k * ndim + 0] = x[k];
        coord[i * nv * ndim + k * ndim + 1] = y[k];
      }
    }
  }

  free(buf);
  return 0;
}

static int re2_boundary(unsigned int *nbcs_, long long **bcs_, int nv,
                        ulong nelgt, MPI_File file, struct comm *c) {
  uint rank = c->id, size = c->np;
  MPI_Comm comm = c->c;

  int ndim = (nv == 4) ? 2 : 3;
  size_t elem_size = nv * ndim * sizeof(double) + sizeof(double);
  size_t header_size = GC_RE2_HEADER_LEN + sizeof(float);

  // Calculate offset for the curve side data.
  MPI_Offset curve_off = header_size + nelgt * elem_size;

  MPI_Status st;
  char bufL[16];
  if (rank == 0)
    MPI_File_read_at(file, curve_off, bufL, sizeof(long), MPI_BYTE, &st);
  check_mpi_call(MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm), c);

  double t;
  READ_T(&t, bufL, long, 1);
  long ncurves = t;

  // Calculate offset for boundary conditions data.
  MPI_Offset bndry_off = curve_off + sizeof(long) + sizeof(long) * 8 * ncurves;
  if (rank == 0)
    MPI_File_read_at(file, bndry_off, bufL, sizeof(long), MPI_BYTE, &st);
  check_mpi_call(MPI_Bcast(bufL, sizeof(long), MPI_BYTE, 0, comm), c);

  READ_T(&t, bufL, long, 1);
  long nbcs = t;

  int nbcsl = nbcs / size;
  int nrem = nbcs - nbcsl * size;
  nbcsl += (rank > (size - 1 - nrem) ? 1 : 0);

  slong out[2][1], bfr[2][1], in = nbcsl;
  comm_scan(out, c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0];

  size_t read_size = nbcsl * sizeof(long) * nv;
  char *buf = calloc(read_size, sizeof(char));
  MPI_Offset offset = bndry_off + sizeof(long) + start * nv * sizeof(long);
  check_mpi_call(
      MPI_File_read_at_all(file, offset, buf, read_size, MPI_BYTE, &st), c);

  struct face_t {
    long eid, fid, e1, e2;
  };

  struct array bfaces;
  array_init(struct face_t, &bfaces, 10);

  double tmp[5];
  char cbc[4], *buf0 = buf;
  struct face_t bface;
  for (sint i = 0; i < nbcsl; i++) {
    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    bface.eid = (long)tmp[0];

    READ_T(tmp, buf0, long, 1);
    buf0 += sizeof(long);
    bface.fid = (long)tmp[0];

    READ_T(tmp, buf0, long, 5);
    buf0 += 5 * sizeof(long);
    READ_T(cbc, buf0, char, 3);
    buf0 += sizeof(long);
    cbc[3] = '\0';

    if (strcmp(cbc, "P  ") == 0) {
      bface.e1 = (long)tmp[0];
      bface.e2 = (long)tmp[1];
      array_cat(struct face_t, &bfaces, &bface, 1);
    }
  }

  nbcs = *nbcs_ = bfaces.n;
  long long *bcs = *bcs_ = tcalloc(long long, 4 * nbcs);

  struct face_t *ptr = bfaces.ptr;
  for (sint i = 0; i < nbcs; i++) {
    bcs[4 * i + 0] = ptr[i].eid;
    bcs[4 * i + 1] = ptr[i].fid;
    bcs[4 * i + 2] = ptr[i].e1;
    bcs[4 * i + 3] = ptr[i].e2;
  }

  array_free(&bfaces);
  free(buf);
  return 0;
}

static int read_geometry(unsigned *nelt, unsigned *nv, double **coord,
                         unsigned *nbcs, long long **bcs, char *fname,
                         struct comm *c) {
  MPI_Info info;
  check_mpi_call(MPI_Info_create(&info), c);

  MPI_File file;
  check_mpi_call(MPI_File_open(c->c, fname, MPI_MODE_RDONLY, info, &file), c);

  ulong nelgt, nelgv;
  int err = re2_header(nelt, nv, &nelgt, &nelgv, file, c);
  err |= re2_coord(coord, *nelt, *nv, file, c);
  err |= re2_boundary(nbcs, bcs, *nv, nelgt, file, c);

  check_mpi_call(MPI_File_close(&file), c);
  check_mpi_call(MPI_Info_free(&info), c);

  return err;
}

static int read_connectivity(unsigned int *nelt_, unsigned *nv_,
                             long long **vl_, char *fname, struct comm *c) {
  uint rank = c->id, size = c->np;
  MPI_Comm comm = c->c;

  MPI_File file;
  check_mpi_call(
      MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &file), c);

  char *buf = (char *)calloc(GC_CO2_HEADER_LEN + 1, sizeof(char));
  MPI_Status st;
  check_mpi_call(MPI_File_read_all(file, buf, GC_CO2_HEADER_LEN, MPI_BYTE, &st),
                 c);

  long long nelgt, nelgv;
  unsigned nv;
  // TODO: Assert version
  char version[6];
  sscanf(buf, "%5s %12lld %12lld %u", version, &nelgt, &nelgv, &nv);

  uint nelt = nelgt / size, nrem = nelgt - nelt * size;
  nelt += (rank > (size - 1 - nrem) ? 1 : 0);

  *nv_ = nv, *nelt_ = nelt;

  float byte_test = 0;
  check_mpi_call(MPI_File_read_all(file, &byte_test, 4, MPI_BYTE, &st), c);
  if (fabs(byte_test - 6.543210) > 1e-7) return 1;

  size_t read_size = nelt * (nv + 1) * sizeof(int);
  size_t header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (rank == 0) read_size += header_size;

  buf = (char *)realloc(buf, read_size * sizeof(char));
  check_mpi_call(MPI_File_read_ordered(file, buf, read_size, MPI_BYTE, &st), c);
  check_mpi_call(MPI_File_close(&file), c);

  char *buf0 = buf + (rank == 0) * header_size;
  long long *vl = *vl_ = tcalloc(long long, nv *nelt);
  int tmp1, tmp2;
  for (uint i = 0; i < nelt; i++) {
    READ_T(&tmp1, buf0, int, 1);
    buf0 += sizeof(int);
    for (unsigned j = 0; j < nv; j++) {
      READ_T(&tmp2, buf0, int, 1);
      buf0 += sizeof(int);
      vl[i * nv + j] = tmp2;
    }
  }

  free(buf);

  return 0;
}

int parrsb_read_mesh(unsigned *nel, unsigned *nv, double **coord,
                     unsigned *nbcs, long long **bcs, char *name,
                     MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  char geom_name[BUFSIZ + 1];
  strncpy(geom_name, name, BUFSIZ);
  strncat(geom_name, ".re2", 5);

  int err = read_geometry(nel, nv, coord, nbcs, bcs, geom_name, &c);

  comm_free(&c);
  return err;
}

int parrsb_read_conn(unsigned *nel, unsigned *nv, long long **vtx, char *name,
                     MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  char con_name[BUFSIZ + 1];
  strncpy(con_name, name, BUFSIZ);
  strncat(con_name, ".co2", 5);
  int err = read_connectivity(nel, nv, vtx, con_name, &c);

  comm_free(&c);
  return err;
}

int parrsb_dump_con(char *name, unsigned nelt, unsigned nv, long long *vl,
                    MPI_Comm comm) {
  const char version[6] = "#v001";
  const float test = 6.54321;

  struct comm c;
  comm_init(&c, comm);
  uint id = c.id;

  char co2_name[BUFSIZ + 1];
  strncpy(co2_name, name, BUFSIZ);
  strncat(co2_name, ".co2", 5);

  MPI_File file;
  int err = MPI_File_open(comm, co2_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);
  if (err) {
    if (id == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s for writing.\n", __FILE__,
              __LINE__, co2_name);
      fflush(stdout);
    }
    return err;
  }

  slong out[2][1], bfr[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, bfr);
  slong start = out[0][0], nelgt = out[1][0];
  slong nelgv = nelgt;

  int write_size = nelt * (nv + 1) * sizeof(int);
  int header_size = GC_CO2_HEADER_LEN + sizeof(float);
  if (id == 0) write_size += header_size;

  char *buf = (char *)calloc(write_size, sizeof(char));
  char *buf0 = buf;

  if (id == 0) {
    sprintf(buf0, "%5s%12d%12d%12d", version, (int)nelgt, (int)nelgv, nv);
    memset(buf0 + strlen(buf0), ' ', GC_CO2_HEADER_LEN - strlen(buf0));
    buf0[GC_CO2_HEADER_LEN] = '\0';
    buf0 += GC_CO2_HEADER_LEN;
    memcpy(buf0, &test, sizeof(float));
    buf0 += sizeof(float);
  }

  for (unsigned i = 0; i < nelt; i++) {
    int temp = start + i + 1;
    WRITE_INT(buf0, temp);
    buf0 += sizeof(int);
    for (unsigned j = 0; j < nv; j++) {
      temp = vl[i * nv + j];
      WRITE_INT(buf0, temp);
      buf0 += sizeof(int);
    }
  }

  MPI_Status st;
  err |= MPI_File_write_ordered(file, buf, write_size, MPI_BYTE, &st);
  err |= MPI_File_close(&file);

  parrsb_barrier(&c);
  comm_free(&c);

  free(buf);

  return err;
}

int parrsb_dump_map(char *name, unsigned nelt, unsigned nv, long long *vtx,
                    MPI_Comm comm) {
  char version[6] = "#v001";
  float test = 6.54321;

  char ma2_name[BUFSIZ + 1];
  strncpy(ma2_name, name, BUFSIZ);
  strncat(ma2_name, ".ma2", 5);

  int nelgt = nelt;
  MPI_Allreduce(&nelt, &nelgt, 1, MPI_INT, MPI_SUM, comm);

  const int npts = nelgt * nv;
  const int depth = (int)log2ll(1.0 * nelgt);
  const int d2 = (int)(pow(2, depth) + 0.5);
  int nactive = nelgt;
  int nrnk = nelgt;
  int noutflow = 0;

  char header[HEADER_LEN];
  sprintf(header, "%5s%12d%12d%12d%12d%12d%12d%12d", version, nelgt, nactive,
          depth, d2, npts, nrnk, noutflow);
  memset(header + strlen(header), ' ', HEADER_LEN - strlen(header));

  MPI_Info infoIn;
  MPI_Info_create(&infoIn);
  MPI_Info_set(infoIn, "access_style", "write_once,random");

  MPI_File file;
  int err = MPI_File_open(comm, ma2_name, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                          infoIn, &file);

  int rank;
  MPI_Comm_rank(comm, &rank);

  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d MPI_File_open failed: %s\n", __FILE__, __LINE__,
              ma2_name);
      fflush(stderr);
    }
    return 1;
  }

  int writeSize = 0;
  if (rank == 0) writeSize = HEADER_LEN * sizeof(char) + sizeof(float);
  writeSize += (nv + 1) * nelt * sizeof(int);

  char *buf = (char *)calloc(writeSize, sizeof(char));
  char *buf0 = buf;
  if (rank == 0) {
    memcpy(buf0, header, HEADER_LEN * sizeof(char));
    buf0 += HEADER_LEN * sizeof(char);
    memcpy(buf0, &test, sizeof(float));
    buf0 += sizeof(float);
  }

  int ivtx[8];
  for (uint i = 0; i < nelt; i++) {
    memcpy(buf0, &rank, sizeof(int));
    buf0 += sizeof(int);

    for (unsigned j = 0; j < nv; j++) ivtx[j] = vtx[i * nv + j];
    memcpy(buf0, ivtx, sizeof(int) * nv);
    buf0 += nv * sizeof(int);
  }

  MPI_Status status;
  err = MPI_File_write_ordered(file, buf, writeSize, MPI_BYTE, &status);
  int errs = (err != 0);

  err = MPI_File_close(&file);
  errs += (err != 0);

  MPI_Info_free(&infoIn);
  free(buf);

  return errs;
}

#undef check_mpi_call

#undef HEADER_LEN
#undef GC_CO2_HEADER_LEN
#undef GC_RE2_HEADER_LEN

#undef WRITE_INT
#undef WRITE_T
#undef READ_T
