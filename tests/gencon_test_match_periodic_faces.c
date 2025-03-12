#include "parRSB.h"

#include "con-impl.h"

static int PRE_TO_SYM_VERTEX[GC_MAX_VERTICES] = {0, 1, 3, 2, 4, 5, 7, 6};
static int PRE_TO_SYM_FACE[GC_MAX_FACES] = {2, 1, 3, 0, 4, 5};

int main(int argc, char *argv[]) {
  scalar *pre_coord = NULL, *sym_coord = NULL, *face_coord = NULL;
  slong *bcs = NULL, *gid = NULL;
  sint *bid = NULL;

  MPI_Init(&argc, &argv);

  parrsb_cmd_line_opts *params = parrsb_parse_cmd_opts(argc, argv);

  struct comm c;
  comm_init(&c, MPI_COMM_WORLD);

  uint nelt, nbcs, nv;
  int err = parrsb_read_mesh(&nelt, &nv, NULL, &pre_coord, &nbcs, &bcs,
                             params->mesh, c.c, 1);

  int rank, size;
  MPI_Comm_rank(c.c, &rank);
  MPI_Comm_size(c.c, &size);

  if (size != 2) {
    fprintf(stderr, "This test must be run with exactly 2 MPI processes.\n");
    err = 1;
    goto cleanup;
  }

  if (nv != 8) {
    fprintf(stderr, "This test is only valid for 3D meshes with 8 vertices.\n");
    err = 1;
    goto cleanup;
  }

  uint nvf = nv / 2;
  uint ndim = nvf / 2 + 1;
  uint nf = nbcs;

  if (nelt == 0 || nf == 0) goto match;

  sym_coord = calloc(nelt * nv * ndim, sizeof(scalar));
  for (uint i = 0; i < nelt; i++) {
    size_t offset = i * nv * ndim;
    for (uint v = 0; v < nv; v++) {
      uint sym_v = PRE_TO_SYM_VERTEX[v];
      for (uint d = 0; d < ndim; d++)
        sym_coord[offset + sym_v * ndim + d] = pre_coord[offset + v * ndim + d];
    }
  }
  free(pre_coord);

  slong out[2][1], wrk[2][1], in = nelt;
  comm_scan(out, &c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0];

  bid = calloc(nf, sizeof(sint));
  face_coord = calloc(nf * nvf * ndim, sizeof(scalar));

  uint offset_ = 0;
  for (uint i = 0; i < nf; i++) {
    bid[i] = rank;

    int element_id = bcs[i * 4] - start - 1;
    int face_id = PRE_TO_SYM_FACE[bcs[i * 4 + 1] - 1];
    uint offset = element_id * nv * ndim;
    for (uint v = 0; v < nvf; v++) {
      uint u = faces3D[face_id][v] - 1;
      for (uint d = 0; d < ndim; d++)
        face_coord[offset_ + v * ndim + d] = sym_coord[offset + u * ndim + d];
    }
    offset_ += nvf * ndim;
  }
  free(bcs), free(sym_coord);

  gid = calloc(nf * nv, sizeof(slong));
match:
  parrsb_match_periodic_faces(gid, nf, bid, nvf, ndim, face_coord, params->tol,
                              c.c);

cleanup:
  free(gid), free(bid), free(bcs);
  free(pre_coord), free(sym_coord), free(face_coord);
  parrsb_cmd_opts_free(params);
  comm_free(&c);
  MPI_Finalize();

  return err;
}
