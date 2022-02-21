#include "genmap-impl.h"
#include "mat.h"
#include <ctype.h>
#include <math.h>
#include <occa.h>

#define BLK_SIZE 512

static occaDevice device;
static occaJson props;
static occaKernel zero, copy, sum, scale, dot, norm, add2s1, add2s2, addc;
static occaKernel lplcn00, lplcn01;

// Laplacian and CG work arrays
struct par_mat *M = NULL;
static occaMemory o_adj_off, o_adj_idx, o_adj_val; // Adjacency matrix
static occaMemory o_diag_idx, o_diag_val;          // Diagonal
static occaMemory o_p, o_w, o_r, o_rr, o_x, o_y;   // CG

// Work arrays
size_t wrk_size, nblocks;
static scalar *wrk = NULL;
static occaMemory o_wrk;

struct gs_data *gsh = NULL; // gs for host communication

// Following definitions should match the ones in laplacian.c
struct laplacian {
  int type, nv;
  uint nel;
  void *data;
};

struct csr_laplacian {
  struct par_mat *M;
  struct gs_data *gsh;
  scalar *buf;
};

int occa_init(const char *backend_, int device_id, int platform_id,
              struct comm *c) {
  // Initialize backend
  char backend[BUFSIZ] = {'\0'};
  size_t len = strnlen(backend_, BUFSIZ);
  int i;
  for (i = 0; i < len; i++)
    backend[i] = tolower(backend_[i]);

  char *fmt_ocl = "{mode: 'OpenCL', device_id: %d, platform_id: %d}";
  char *fmt_cuda = "{mode: 'CUDA', device_id: %d}";
  char *fmt_serial = "{mode: 'Serial'}";

  char fmt[BUFSIZ] = {'\0'};
  if (strncmp(backend, "opencl", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_ocl, device_id, platform_id);
  else if (strncmp(backend, "cuda", BUFSIZ) == 0)
    snprintf(fmt, BUFSIZ, fmt_cuda, device_id);
  else
    snprintf(fmt, BUFSIZ, fmt_serial);

  device = occaCreateDeviceFromString(fmt);

  // OCCA Kernels
  char okl[PATH_MAX];
  char *okl_dir = getenv("PARRSB_OKL_DIR");
  if (okl_dir != NULL)
    strncpy(okl, okl_dir, PATH_MAX);
  else if (c->id == 0) {
    printf("PARRSB_OKL_DIR is not defined !\n");
    exit(1);
  }
  strncat(okl, "/occa.okl", PATH_MAX);

  props = occaCreateJson();
  occaJsonObjectSet(props, "defines/BLK_SIZE", occaInt(512));
  occaJsonObjectSet(props, "defines/uint", occaString("unsigned int"));
  occaJsonObjectSet(props, "defines/ulong", occaString("unsigned long"));
  occaJsonObjectSet(props, "defines/scalar", occaString("double"));

  if (c->id == 0) {
    zero = occaDeviceBuildKernel(device, okl, "zero", props);
    copy = occaDeviceBuildKernel(device, okl, "copy", props);
    sum = occaDeviceBuildKernel(device, okl, "sum", props);
    scale = occaDeviceBuildKernel(device, okl, "scale", props);
    dot = occaDeviceBuildKernel(device, okl, "dot", props);
    norm = occaDeviceBuildKernel(device, okl, "norm", props);
    add2s1 = occaDeviceBuildKernel(device, okl, "add2s1", props);
    add2s2 = occaDeviceBuildKernel(device, okl, "add2s2", props);
    addc = occaDeviceBuildKernel(device, okl, "addc", props);
    lplcn00 = occaDeviceBuildKernel(device, okl, "laplacian_v00", props);
    lplcn01 = occaDeviceBuildKernel(device, okl, "laplacian_v01", props);
  }

  comm_barrier(c);

  zero = occaDeviceBuildKernel(device, okl, "zero", props);
  copy = occaDeviceBuildKernel(device, okl, "copy", props);
  sum = occaDeviceBuildKernel(device, okl, "sum", props);
  scale = occaDeviceBuildKernel(device, okl, "scale", props);
  dot = occaDeviceBuildKernel(device, okl, "dot", props);
  norm = occaDeviceBuildKernel(device, okl, "norm", props);
  add2s1 = occaDeviceBuildKernel(device, okl, "add2s1", props);
  add2s2 = occaDeviceBuildKernel(device, okl, "add2s2", props);
  addc = occaDeviceBuildKernel(device, okl, "addc", props);
  lplcn00 = occaDeviceBuildKernel(device, okl, "laplacian_v00", props);
  lplcn01 = occaDeviceBuildKernel(device, okl, "laplacian_v01", props);

  return 0;
}

int occa_lanczos_init(struct comm *c, struct laplacian *l, int niter) {
  // CG arrays
  uint lelt = l->nel;
  o_p = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_w = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_r = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_y = occaDeviceMalloc(device, sizeof(scalar) * lelt, NULL, occaDefault);
  o_rr = occaDeviceMalloc(device, sizeof(scalar) * lelt * (niter + 1), NULL,
                          occaDefault);

  // Work array size
  nblocks = (lelt + BLK_SIZE - 1) / BLK_SIZE;
  wrk_size = lelt * (niter + 1) > nblocks ? lelt * (niter + 1) : nblocks;

  // Laplacian
  struct csr_laplacian *L = (struct csr_laplacian *)l->data;
  M = L->M;
  assert(IS_CSR(M));
  assert(IS_DIAG(M));
  if (M->rn > 0) {
    o_adj_off =
        occaDeviceMalloc(device, sizeof(uint) * (M->rn + 1), NULL, occaDefault);
    occaCopyPtrToMem(o_adj_off, M->adj_off, occaAllBytes, 0, occaDefault);

    uint nadj = M->adj_off[M->rn];
    o_adj_idx =
        occaDeviceMalloc(device, sizeof(uint) * nadj, NULL, occaDefault);
    occaCopyPtrToMem(o_adj_idx, M->adj_idx, occaAllBytes, 0, occaDefault);

    o_adj_val =
        occaDeviceMalloc(device, sizeof(scalar) * nadj, NULL, occaDefault);
    occaCopyPtrToMem(o_adj_val, M->adj_val, occaAllBytes, 0, occaDefault);

    o_diag_val =
        occaDeviceMalloc(device, sizeof(scalar) * M->rn, NULL, occaDefault);
    occaCopyPtrToMem(o_diag_val, M->diag_val, occaAllBytes, 0, occaDefault);

    o_diag_idx =
        occaDeviceMalloc(device, sizeof(uint) * M->rn, NULL, occaDefault);
    occaCopyPtrToMem(o_diag_idx, M->diag_idx, occaAllBytes, 0, occaDefault);

    o_x = occaDeviceMalloc(device, sizeof(scalar) * M->cn, NULL, occaDefault);
    if (wrk_size < M->cn + M->rn)
      wrk_size = M->cn + M->rn;
  }

  wrk = tcalloc(scalar, wrk_size);
  o_wrk =
      occaDeviceMalloc(device, sizeof(scalar) * wrk_size, NULL, occaDefault);

  slong *ids = (slong *)tcalloc(slong, M->cn);
  for (uint i = 0; i < M->cn; i++)
    ids[i] = -M->cols[i];
  for (uint i = 0; i < M->rn; i++)
    ids[M->diag_idx[i]] *= -1;

  gsh = gs_setup(ids, M->cn, c, 0, gs_crystal_router, 0);
  free(ids);

  return 0;
}

static int vec_ortho(occaMemory o_in, ulong nelg, uint lelt, struct comm *c) {
  occaKernelRun(sum, occaUInt(lelt), o_in, o_wrk);
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  scalar tmp = 0.0, rtr;
  uint i;
  for (i = 0; i < nblocks; i++)
    tmp += wrk[i];
  comm_allreduce(c, gs_double, gs_add, &tmp, 1, &rtr);
  tmp /= nelg;

  occaKernelRun(addc, o_in, occaDouble(-1.0 * tmp), occaInt(lelt));

  return 0;
}

static scalar vec_norm(occaMemory o_a, uint lelt, struct comm *c) {
  occaKernelRun(norm, o_wrk, occaInt(lelt), o_a);
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  scalar pp = 0.0, tmp;
  uint i;
  for (i = 0; i < nblocks; i++)
    pp += wrk[i];
  comm_allreduce(c, gs_double, gs_add, &pp, 1, &tmp);

  return pp;
}

static scalar vec_dot(occaMemory o_a, occaMemory o_b, uint lelt,
                      struct comm *c) {
  occaKernelRun(dot, o_wrk, occaInt(lelt), o_a, o_b);
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  scalar pap = 0.0, tmp;
  uint i;
  for (i = 0; i < nblocks; i++)
    pap += wrk[i];
  comm_allreduce(c, gs_double, gs_add, &pap, 1, &tmp);

  return pap;
}

static scalar vec_normalize(occaMemory o_scaled, uint indx, occaMemory o_a,
                            uint lelt, struct comm *c) {
  occaKernelRun(norm, o_wrk, occaInt(lelt), o_a);
  occaCopyMemToPtr(wrk, o_wrk, occaAllBytes, 0, occaDefault);

  scalar rtr = 0, tmp;
  uint i;
  for (i = 0; i < nblocks; i++)
    rtr += wrk[i];
  comm_allreduce(c, gs_double, gs_add, &rtr, 1, &tmp);
  scalar rnorm = sqrt(rtr);
  scalar rni = 1.0 / rnorm;

  occaKernelRun(scale, o_scaled, occaUInt(indx), o_a, occaDouble(rni),
                occaUInt(lelt));

  return rnorm;
}

static void laplacian_op(occaMemory o_w, occaMemory o_p, uint lelt,
                         occaMemory o_wrk, buffer *bfr) {
  occaCopyMemToPtr((void *)(wrk + M->cn), o_p, occaAllBytes, 0, occaDefault);

  for (uint i = 0; i < M->rn; i++)
    wrk[M->diag_idx[i]] = wrk[M->cn + i];
  gs(wrk, gs_double, gs_add, 0, gsh, bfr);

  // Fix this: occaAllbytes --> M->cn * sizeof(scalar)
  occaCopyPtrToMem(o_wrk, wrk, sizeof(scalar) * M->cn, 0, occaDefault);

  occaKernelRun(lplcn01, o_w, occaUInt(lelt), o_adj_off, o_adj_idx, o_adj_val,
                o_diag_idx, o_diag_val, o_wrk);
}

int occa_lanczos_aux(genmap_vector diag, genmap_vector upper, genmap_vector *rr,
                     uint lelt, ulong nelg, int niter, genmap_vector f,
                     struct laplacian *gl, struct comm *gsc, buffer *bfr) {
  assert(f->size == lelt);

  // vec_create_zeros(p, ...)
  occaKernelRun(zero, o_p, occaUInt(lelt));

  // vec_copy(r, f)
  occaCopyPtrToMem(o_r, f->data, occaAllBytes, 0, occaDefault);

  // vec_ortho(gsc, r, nelg);
  vec_ortho(o_r, nelg, lelt, gsc);

  // vec_scale(rr[0], r, rni);
  scalar rtr = vec_norm(o_r, lelt, gsc);
  scalar rnorm = sqrt(rtr);
  scalar rni = 1.0 / rni;
  occaKernelRun(scale, o_rr, occaUInt(0), o_r, occaDouble(rni), occaUInt(lelt));

  scalar eps = 1.e-5;
  scalar rtol = rnorm * eps;

  scalar rtz1 = 1.0, rtz2;
  scalar pap = 0.0, pap_old;
  scalar alpha, beta;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    // add2s1(p, r, beta, n)
    occaKernelRun(add2s1, o_p, o_r, occaDouble(beta), occaUInt(lelt));

    scalar pp = vec_norm(o_p, lelt, gsc);

    // orthogonalize
    vec_ortho(o_p, nelg, lelt, gsc);

    // laplacian
    laplacian_op(o_w, o_p, lelt, o_wrk, bfr);

    // pap = vec_dot(w, p);
    pap_old = pap;
    pap = vec_dot(o_w, o_p, lelt, gsc);

#if 1
    if (gsc->id == 0)
      printf("occa iter = %d beta = %lf pp = %lf pap = %lf\n", iter, beta, pp,
             pap);
#endif

    alpha = rtz1 / pap;

    // add2s2(r, w, -1.0 * alpha, n);
    occaKernelRun(add2s2, o_r, o_w, occaDouble(-1.0 * alpha), occaInt(lelt));

    // vec_scale(rr[iter + 1], r, rni);
    rtr = vec_norm(o_r, lelt, gsc);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;
    occaKernelRun(scale, o_rr, occaUInt((iter + 1) * lelt), o_r,
                  occaDouble(rni), occaUInt(lelt));

    if (iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  occaCopyMemToPtr(wrk, o_rr, occaAllBytes, 0, occaDefault);
  uint i;
  for (i = 0; i < iter + 1; i++)
    memcpy(rr[i]->data, &wrk[i * lelt], sizeof(scalar) * lelt);

  metric_acc(TOL_FINAL, rnorm);
  metric_acc(TOL_TARGET, rtol);

  return iter;
}

int occa_lanczos_free() {
  occaFree(&o_p);
  occaFree(&o_w);
  occaFree(&o_r);
  occaFree(&o_y);
  occaFree(&o_x);
  occaFree(&o_rr);
  occaFree(&o_wrk);
  occaFree(&o_adj_off);
  occaFree(&o_adj_idx);
  occaFree(&o_adj_val);
  occaFree(&o_diag_idx);
  occaFree(&o_diag_val);

  if (wrk != NULL)
    free(wrk);
  if (gsh != NULL)
    gs_free(gsh);

  return 0;
}

int occa_free() {
  occaFree(&copy);
  occaFree(&sum);
  occaFree(&scale);
  occaFree(&dot);
  occaFree(&norm);
  occaFree(&add2s1);
  occaFree(&add2s2);
  occaFree(&addc);

  occaFree(&device);

  return 0;
}
