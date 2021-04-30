#include <genmap-impl.h>
#include <genmap-multigrid.h>

int log2i(sint i) {
  sint k = 1, l = 0;
  while (k <= i)
    k *= 2, l++;
  return l - 1;
}

void mg_check(genmap_handle h, struct comm *c) {
  struct mg_data csr, gs;
  mg_setup(h, c, 0, &csr);
  mg_setup(h, c, 1, &gs);

  int nlevels = mg_get_nlevels(&csr);
  assert(nlevels == mg_get_nlevels(&gs));

  uint *csr_off = mg_get_level_off(&csr);
  uint *gs_off = mg_get_level_off(&gs);

  uint i;
  for (i = 0; i < nlevels + 1; i++)
    assert(csr_off[i] == gs_off[i]);

  uint n1 = csr_off[1] - csr_off[0];
  uint n2 = csr_off[2] - csr_off[1];

  GenmapScalar *csr_r, *gs_r;
  GenmapCalloc(2 * n1, &csr_r);
  GenmapCalloc(2 * n1, &gs_r);

  for (i = 0; i < n1; i++)
    gs_r[i] = csr_r[i] = i + 1.0;

  mg_restrict(csr_r, 0, &csr, &h->buf);
  mg_restrict(gs_r, 0, &gs, &h->buf);

  for (i = 0; i < n1 + n2; i++) {
    if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
      printf("Error 1.0: %d %lf %lf %u %u %u\n", c->id, gs_r[i], csr_r[i], i,
             n1, n2);
      assert(0);
    }
  }

  for (i = 0; i < n1 + n2; i++)
    gs_r[i] = csr_r[i] = i + 1.0;

  mg_interpolate(csr_r, 0, &csr, &h->buf);
  mg_interpolate(gs_r, 0, &gs, &h->buf);

  for (i = 0; i < n1 + n2; i++) {
    if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
      printf("Error 1.1: %lf %lf %u %u %u\n", c->id, gs_r[i], csr_r[i], i, n1,
             n2);
      assert(0);
    }
  }

  for (i = 0; i < n1; i++)
    gs_r[i] = csr_r[i] = i + 1.0;

  mg_diagonal_scaling(csr_r + n1, csr_r, 0.2, 0, &csr);
  mg_diagonal_scaling(gs_r + n1, gs_r, 0.2, 0, &gs);

  for (i = 0; i < 2 * n1; i++) {
    if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
      printf("Error 2.0: %d %lf %lf %u %u %u\n", c->id, gs_r[i], csr_r[i], i,
             n1, n2);
      assert(0);
    }
  }

  // for (i = 0; i < n2; i++)
  //   gs_r[i] = csr_r[i] = i + 1.0;

  // mg_diagonal_scaling(csr_r + n2, csr_r, 0.2, 1, &csr);
  // mg_diagonal_scaling(gs_r + n2, gs_r, 0.2, 1, &gs);

  // for (i = 0; i < 2 * n2; i++) {
  //   if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
  //     printf("Error 2.1: %d %lf %lf %u %u %u\n", c->id, gs_r[i], csr_r[i], i,
  //     n1,
  //            n2);
  //     assert(0);
  //   }
  // }

  for (i = 0; i < n1; i++)
    csr_r[i] = gs_r[i] = 1.0;

  mg_operator(csr_r + n1, csr_r, 0, &csr, &h->buf);
  mg_operator(gs_r + n1, gs_r, 0, &gs, &h->buf);

  for (i = n1; i < 2 * n1; i++) {
    if (fabs(gs_r[i]) > 1e-10 || fabs(csr_r[i]) > 1e-10) {
      printf("Error 3.0: %d %lf %lf %u %u\n", c->id, gs_r[i], csr_r[i], i, n1);
      assert(0);
    }
  }

  for (i = 0; i < n1; i++)
    csr_r[i] = gs_r[i] = i + 1.0;

  mg_operator(csr_r + n1, csr_r, 0, &csr, &h->buf);
  mg_operator(gs_r + n1, gs_r, 0, &gs, &h->buf);

  for (i = n1; i < 2 * n1; i++) {
    if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
      printf("Error 3.1: %d %lf %lf %u %u\n", c->id, gs_r[i], csr_r[i], i, n1);
      assert(0);
    }
  }

  for (i = 0; i < n2; i++)
    csr_r[i] = gs_r[i] = i + 1.0;

  mg_operator(csr_r + n2, csr_r, 1, &csr, &h->buf);
  mg_operator(gs_r + n2, gs_r, 1, &gs, &h->buf);

  for (i = n2; i < 2 * n2; i++) {
    if (fabs(gs_r[i] - csr_r[i]) > 1e-10) {
      printf("Error 3.2: %d %lf %lf %u %u\n", c->id, gs_r[i], csr_r[i], i, n2);
      assert(0);
    }
  }

  GenmapFree(csr_r);
  GenmapFree(gs_r);
}

void mg_setup(genmap_handle h, struct comm *c, int type, struct mg_data *d) {
  if (type == 0)
    mg_setup_csr(h, c, d);
  else if (type == 1)
    mg_setup_gs(h, c, d);
  d->h = h;
}

int mg_get_nlevels(struct mg_data *d) { return d->get_nlevels(d); }

uint *mg_get_level_off(struct mg_data *d) { return d->get_level_off(d); }

int mg_get_nsmooth(struct mg_data *d, int level) {
  return d->get_nsmooth(d, level);
}

GenmapScalar mg_get_sigma(struct mg_data *d, int level) {
  return d->get_sigma(d, level);
}

void mg_diagonal_scaling(GenmapScalar *v, GenmapScalar *u, GenmapScalar sigma,
                         int level, struct mg_data *d) {
  d->diagonal_scaling(d, level, v, u, sigma);
}

void mg_operator(GenmapScalar *v, GenmapScalar *u, int level, struct mg_data *d,
                 buffer *buf) {
  d->G(d, level, v, u, buf);
}

void mg_restrict(GenmapScalar *v, int level, struct mg_data *d, buffer *buf) {
  d->rstrct(d, level, v, buf);
}

void mg_interpolate(GenmapScalar *v, int level, struct mg_data *d,
                    buffer *buf) {
  d->intrp(d, level, v, buf);
}

void mg_coarse_solve(GenmapScalar *u, GenmapScalar *r, struct mg_data *d) {
  d->coarse(d, u, r);
}

void mg_free(struct mg_data *d) { d->free(d); }
