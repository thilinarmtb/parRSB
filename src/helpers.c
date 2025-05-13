#include "parrsb_impl.h"

#include <stdarg.h>
#include <sys/resource.h>

int nv_to_ndim(int nv) { return (nv == 8) ? 3 : 2; }

void parrsb_barrier(struct comm *c) {
#if defined(PARRSB_SYNC_BY_REDUCTION)
  sint dummy = c->id, wrk;
  comm_allreduce(c, gs_int, gs_max, &dummy, 1, &wrk);
#else
  comm_barrier(c);
#endif
}

void parrsb_print(const struct comm *c, int verbose, const char *fmt, ...) {
  comm_barrier(c);

  va_list vargs;
  if (c->id == 0 && verbose > 0) {
    va_start(vargs, fmt);
    vprintf(fmt, vargs);
    va_end(vargs);
    printf("\n");
    fflush(stdout);
  }
}

#if defined __GLIBC__

#include <execinfo.h>

// Obtain a backtrace and print it to stdout.
void parrsb_print_stack(void) {
  void *bt[50];
  int bt_size = backtrace(bt, 50);
  if (bt_size == 0) {
    fprintf(stderr, "backtrace(): Obtained 0 stack frames.\n");
    return;
  }

  char **symbols = backtrace_symbols(bt, bt_size);
  fprintf(stderr, "backtrace(): obtained %d stack frames.\n", bt_size);
  for (unsigned i = 0; i < (unsigned)bt_size; i++)
    fprintf(stderr, "%s\n", symbols[i]);
  free(symbols);
}
#else
void parrsb_print_stack(void) {}
#endif // defined __GLIBC__

int log2ll(long long n) {
  int k = 0;
  while (n > 1) n /= 2, k++;

  return k;
}

#if defined(PARRSB_BLAS)

#define fdgetrf GS_FORTRAN_UNPREFIXED(dgetrf, DGETRF)
void fdgetrf(int *M, int *N, double *A, int *lda, int *IPIV, int *info);

#define fdgetri GS_FORTRAN_UNPREFIXED(dgetri, DGETRI)
void fdgetri(int *N, double *A, int *lda, int *ipiv, double *work, int *lwork,
             int *info);

#define fsgetrf GS_FORTRAN_UNPREFIXED(sgetrf, SGETRF)
void fsgetrf(int *M, int *N, float *A, int *lda, int *IPIV, int *info);

#define fsgetri GS_FORTRAN_UNPREFIXED(sgetri, SGETRI)
void fsgetri(int *N, float *A, int *lda, int *ipiv, float *work, int *lwork,
             int *info);

int matrix_inverse(int N, scalar *A) {
  int size = N * N;
  int info;

  int *ipiv = (int *)calloc(N, sizeof(int));
  scalar *work = (scalar *)calloc(N * N, sizeof(scalar));

  if (sizeof(scalar) == sizeof(double))
    fdgetrf(&N, &N, A, &N, ipiv, &info);
  else if (sizeof(scalar) == sizeof(float))
    fsgetrf(&N, &N, (float *)A, &N, ipiv, &info);
  else
    info = N + 1;

  if (info != 0) goto cleanup;

  if (sizeof(scalar) == sizeof(double))
    fdgetri(&N, A, &N, ipiv, work, &size, &info);
  else if (sizeof(scalar) == sizeof(float))
    fsgetri(&N, (float *)A, &N, ipiv, (float *)work, &size, &info);
  else
    info = N + 1;

cleanup:
  free(ipiv);
  free(work);
  return info;
}

#else
int matrix_inverse(int N, scalar *A) {
  (void)A;
  fprintf(stderr, "Error: Compile parRSB with BLAS enabled!\n");
  return N + 1;
}
#endif // PARRSB_BLAS

int parrsb_dist_mesh(unsigned int *nelt_, long long **vl_, double **coord_,
                     int *part, unsigned nv, MPI_Comm comm) {
  typedef struct {
    int proc;
    long long vtx[MAXNV];
    double coord[MAXNV * MAXDIM];
  } elem_data;

  uint nelt = *nelt_;
  struct array elements;
  array_init(elem_data, &elements, nelt);

  elem_data data;
  long long *vl = *vl_;
  uint e, n;
  for (e = 0; e < nelt; ++e) {
    data.proc = part[e];
    for (n = 0; n < nv; ++n) data.vtx[n] = vl[e * nv + n];
    array_cat(elem_data, &elements, &data, 1);
  }
  assert(elements.n == nelt);

  unsigned ndim = nv_to_ndim(nv);
  elem_data *ed = elements.ptr;
  double *coord = (coord_ == NULL ? NULL : *coord_);
  if (coord != NULL) {
    for (e = 0; e < nelt; e++)
      for (n = 0; n < ndim * nv; n++) ed[e].coord[n] = coord[e * ndim * nv + n];
  }

  struct comm c;
  comm_init(&c, comm);

  struct crystal cr;
  crystal_init(&cr, &c);

  sarray_transfer(elem_data, &elements, proc, 0, &cr);
  *nelt_ = nelt = elements.n;
  ed = elements.ptr;

  vl = *vl_ = (long long *)realloc(*vl_, nv * nelt * sizeof(long long));
  for (e = 0; e < nelt; ++e)
    for (n = 0; n < nv; ++n) vl[e * nv + n] = ed[e].vtx[n];

  if (coord != NULL) {
    coord = *coord_ =
        (double *)realloc(*coord_, ndim * nv * nelt * sizeof(double));
    for (e = 0; e < nelt; ++e) {
      for (n = 0; n < ndim * nv; ++n) coord[e * ndim * nv + n] = ed[e].coord[n];
    }
  }

  crystal_free(&cr);
  comm_free(&c);
  array_free(&elements);

  return 0;
}

void parrsb_get_part_stat(int *nc, int *ns, int *nss, int *nel, long long *vtx,
                          int nelt, int nv, MPI_Comm ce) {
  struct comm comm;
  comm_init(&comm, ce);

  uint np = comm.np;
  if (np == 1) return;

  size_t Npts = nelt * nv;
  slong *data = (slong *)malloc((Npts + 1) * sizeof(slong));
  for (size_t i = 0; i < Npts; i++) data[i] = vtx[i];
  struct gs_data *gsh = gs_setup(data, Npts, &comm, 0, gs_pairwise, 0);

  int Nmsg;
  pw_data_nmsg(gsh, &Nmsg);

  int *Ncomm = (int *)malloc((Nmsg + 1) * sizeof(int));
  pw_data_size(gsh, Ncomm);

  gs_free(gsh);
  free(data);

  sint nelMin, nelMax, nelSum;
  sint ncMin, ncMax, ncSum;
  sint nsMin, nsMax, nsSum;
  sint nssMin, nssMax, nssSum;
  sint b;

  ncMax = Nmsg;
  ncMin = Nmsg;
  ncSum = Nmsg;
  comm_allreduce(&comm, gs_int, gs_max, &ncMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &ncMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &ncSum, 1, &b);

  nsMax = Ncomm[0];
  nsMin = Ncomm[0];
  nsSum = Ncomm[0];
  for (int i = 1; i < Nmsg; ++i) {
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

  nelMax = nelt;
  nelMin = nelt;
  nelSum = nelt;
  comm_allreduce(&comm, gs_int, gs_max, &nelMax, 1, &b);
  comm_allreduce(&comm, gs_int, gs_min, &nelMin, 1, &b);
  comm_allreduce(&comm, gs_int, gs_add, &nelSum, 1, &b);

  free(Ncomm);
  comm_free(&comm);

  if (nc != NULL) {
    nc[0] = ncMin;
    nc[1] = ncMax;
    nc[2] = ncSum;
  }

  if (ns != NULL) {
    ns[0] = nsMin;
    ns[1] = nsMax;
    ns[2] = nsSum;
  }

  if (nss != NULL) {
    nss[0] = nssMin;
    nss[1] = nssMax;
    nss[2] = nssSum;
  }

  if (nel != NULL) {
    nel[0] = nelMin;
    nel[1] = nelMax;
    nel[2] = nelSum;
  }
}

#define WRITE_T(dest, val, T, nunits)                                          \
  do {                                                                         \
    memcpy(dest, (val), sizeof(T) * nunits);                                   \
    dest += sizeof(T) * nunits;                                                \
  } while (0)

int parrsb_vector_dump(const char *fname, scalar *y, struct rsb_element *elm,
                       uint nelt, unsigned nv, struct comm *c) {
  MPI_File file;
  sint err = MPI_File_open(c->c, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                           MPI_INFO_NULL, &file);
  slong wrk[2][1];
  comm_allreduce(c, gs_int, gs_add, &err, 1, wrk);

  uint rank = c->id;
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error opening file: %s\n", __FILE__, __LINE__,
              fname);
      fflush(stderr);
    }
    return err;
  }

  slong out[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong nelgt = out[1][0];

  int ndim = nv_to_ndim(nv);
  uint write_size = ((ndim + 1) * sizeof(double) + sizeof(slong)) * nelt;
  if (rank == 0) write_size += sizeof(long) + sizeof(int); // for nelgt and ndim

  char *bfr, *bfr0;
  bfr = bfr0 = (char *)calloc(write_size, sizeof(char));
  if (rank == 0) {
    WRITE_T(bfr0, &nelgt, slong, 1);
    WRITE_T(bfr0, &ndim, int, 1);
  }

  uint i;
  for (i = 0; i < nelt; i++) {
    WRITE_T(bfr0, &elm[i].globalId, ulong, 1);
    WRITE_T(bfr0, &elm[i].coord[0], double, ndim);
    WRITE_T(bfr0, &y[i], double, 1);
  }

  MPI_Status st;
  err = MPI_File_write_ordered(file, bfr, write_size, MPI_BYTE, &st);
  comm_allreduce(c, gs_int, gs_add, &err, 1, wrk);
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error writing file: %s.\n", __FILE__, __LINE__,
              fname);
      fflush(stdout);
    }
    return err;
  }

  err = MPI_File_close(&file);
  comm_scan(out, c, gs_int, gs_add, &err, 1, wrk);
  if (err) {
    if (rank == 0) {
      fprintf(stderr, "%s:%d Error closing file: %s.\n", __FILE__, __LINE__,
              fname);
      fflush(stdout);
    }
    return err;
  }

  parrsb_barrier(c);
  free(bfr);
  return err;
}

#undef WRITE_T
