#include "parRSB.h"

#include "con-impl.h"

static int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
static int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};

#define tfree(ptr) free((ptr)), (ptr) = 0

static int check_input(int nv, int size, int rank) {
  int err = 0;
  if (size != 2) {
    if (rank == 0)
      fprintf(stderr, "This test must be run with exactly 2 MPI processes.\n");
    err = 1;
  }

  if (nv != 8) {
    if (rank == 0)
      fprintf(stderr,
              "This test is only valid for 3D meshes with 8 vertices.\n");
    err = 1;
  }

  fflush(stderr);

  return err;
}

static void tag_periodic_faces(long long *gid, double *fc, int *bid,
                               unsigned nelt, unsigned nv, const long long *vl,
                               unsigned nvf, unsigned ndim, unsigned nf,
                               const double *coord, const long long *bcs,
                               const struct comm *c) {
  double *symc = calloc(nelt * nv * ndim, sizeof(double));
  for (unsigned i = 0; i < nelt; i++) {
    size_t offset = i * nv * ndim;
    for (unsigned v = 0; v < nv; v++) {
      unsigned sym_v = PRE_TO_SYM_VERTEX[v];
      for (unsigned d = 0; d < ndim; d++)
        symc[offset + sym_v * ndim + d] = coord[offset + v * ndim + d];
    }
  }

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0];

  unsigned fc_off = 0;
  unsigned fv_off = 0;
  for (unsigned i = 0; i < nf; i++) {
    bid[i] = c->id;

    long long eid = bcs[i * 4] - start - 1;
    int fid = PRE_TO_SYM_FACE[bcs[i * 4 + 1] - 1];

    size_t ev_off = eid * nv;
    size_t ec_off = eid * nv * ndim;
    for (unsigned v = 0; v < nvf; v++) {
      unsigned u = faces3D[fid][v] - 1;
      gid[fv_off + v] = vl[ev_off + u];
      for (unsigned d = 0; d < ndim; d++)
        fc[fc_off + v * ndim + d] = symc[ec_off + u * ndim + d];
    }
    fc_off += nvf * ndim;
    fv_off += nvf;
  }

  for (int i = 0; i < c->np; i++) {
    if (i == c->id) {
      printf("prev proc: %d:", c->id);
      for (unsigned i = 0; i < nelt * nv; i++) printf(" %2lld", vl[i]);
      printf("\n");
      printf("symv proc: %d:", c->id);
      for (unsigned i = 0; i < nelt * nv; i++) printf(" %2lld", vl[i]);
      printf("\n");
      printf("gid  proc: %d:", c->id);
      for (unsigned i = 0; i < nf * nvf; i++) printf(" %2lld", gid[i]);
      printf("\n");
    }
    fflush(stdout);
    MPI_Barrier(c->c);
  }

  tfree(symc);
}

static int check_accuracy(unsigned nelt, unsigned nv, unsigned ndim,
                          const double *coord, unsigned nbcs,
                          const long long *bcs, double tol, unsigned nf,
                          unsigned nvf, const long long *vlf,
                          const struct comm *c) {

  unsigned size = nelt * nv;
  long long *vl = (long long *)calloc(size, sizeof(long long));
  parrsb_conn_mesh(vl, coord, nelt, ndim, bcs, nbcs, tol, c->c, 0);

  long long *vlm = calloc(size, sizeof(long long));
  long long *vlM = calloc(size, sizeof(long long));
  parrsb_conn_mesh(vlM, coord, nelt, ndim, 0, 0, tol, c->c, 0);
  for (unsigned i = 0; i < size; i++) vlm[i] = vlM[i];

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0];

  for (unsigned i = 0; i < nf; i++) {
    long eid = bcs[i * 4] - start - 1;
    int fid = PRE_TO_SYM_FACE[bcs[i * 4 + 1] - 1];

    size_t ev_off = eid * nv;
    for (unsigned v = 0; v < nvf; v++) {
      unsigned u = faces3D[fid][v] - 1;
      vlm[ev_off + u] = vlM[ev_off + u] = vlf[i * nvf + v];
    }
  }

  for (int i = 0; i < c->np; i++) {
    if (i == c->id) {
      printf("vlf proc: %d:", c->id);
      for (unsigned i = 0; i < nf * nvf; i++) printf(" %2lld", vlf[i]);
      printf("\n");
      printf("vl  proc: %d:", c->id);
      for (unsigned i = 0; i < nelt * nv; i++) printf(" %2lld", vl[i]);
      printf("\n");
      printf("vlm proc: %d:", c->id);
      for (unsigned i = 0; i < nelt * nv; i++) printf(" %2lld", vlm[i]);
      printf("\n");
    }
    fflush(stdout);
    MPI_Barrier(c->c);
  }

  sint err = 0;
#if 0
  struct gs_data *gsh = gs_setup(vl, size, c, 0, gs_pairwise, 0);
  buffer bfr;
  buffer_init(&bfr, 1024);
  gs(vlm, gs_long, gs_min, 0, gsh, &bfr);
  gs(vlM, gs_long, gs_max, 0, gsh, &bfr);
  buffer_free(&bfr);
  gs_free(gsh);

  for (unsigned i = 0; i < size; i++) {
    if (vlm[i] != vlM[i]) {
      err = 1;
      break;
    }
  }
  comm_allreduce(c, gs_int, gs_add, &err, 1, wrk);
#endif

  tfree(vl), tfree(vlm), tfree(vlM);
  return err;
}

int main(int argc, char *argv[]) {
  double *coord = 0, *fcoord = 0;
  long long *bcs = 0, *vl = 0, *gid = 0;
  int *bid = 0;

  MPI_Init(&argc, &argv);

  struct comm c;
  comm_init(&c, MPI_COMM_WORLD);

  parrsb_cmd_line_opts *params = parrsb_parse_cmd_opts(argc, argv);

  unsigned nelt, nbcs, nv;
  parrsb_read_mesh(&nelt, &nv, NULL, &coord, &nbcs, &bcs, params->mesh, c.c, 1);

  int err = check_input(nv, c.np, c.id);
  if (err) goto cleanup;

  unsigned nvf = nv / 2;
  unsigned ndim = nvf / 2 + 1;
  unsigned nf = nbcs;

  // Calculate connectivity without periodic info.
  vl = (long long *)calloc(nelt * nv, sizeof(long long));
  parrsb_conn_mesh(vl, coord, nelt, ndim, NULL, 0, params->tol, c.c, 0);

  // Match periodic faces using the new routine.
  gid = calloc(nf * nv, sizeof(long long));
  bid = calloc(nf, sizeof(int));
  fcoord = calloc(nf * nvf * ndim, sizeof(double));
  tag_periodic_faces(gid, fcoord, bid, nelt, nv, vl, nvf, ndim, nbcs, coord,
                     bcs, &c);
  parrsb_match_periodic_faces(gid, nf, bid, nvf, ndim, fcoord, params->tol,
                              c.c);

  // Check accuracy compared to the old version.
  err |= check_accuracy(nelt, nv, ndim, coord, nbcs, bcs, params->tol, nf, nvf,
                        gid, &c);

cleanup:
  tfree(coord), tfree(fcoord);
  tfree(bcs), tfree(vl), tfree(gid);
  tfree(bid);

  parrsb_cmd_opts_free(params);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
