#include <genmap-impl.h>
#include <genmap-multigrid.h>

void vcycle(GenmapScalar *u1, GenmapScalar *rhs, struct mg_data *d,
            buffer *buf) {
  int nlevels = mg_get_nlevels(d);
  uint *lvl_off = mg_get_level_off(d);

  // TODO: Allocate from buffer space?
  uint off = lvl_off[nlevels];
  GenmapScalar *s;
  GenmapMalloc(off, &s);
  GenmapScalar *Gs;
  GenmapMalloc(off, &Gs);
  GenmapScalar *r;
  GenmapMalloc(off, &r);
  GenmapScalar *u;
  GenmapMalloc(off, &u);

  uint i;
  for (i = 0; i < off; i++)
    s[i] = Gs[i] = r[i] = u[i] = 0.0;
  for (i = 0; i < lvl_off[1]; i++)
    r[i] = rhs[i];

  GenmapScalar *diag;
  GenmapScalar sigma;
  uint n, j;
  int nsmooth, lvl;
  for (lvl = 0; lvl < nlevels - 1; lvl++) {
    off = lvl_off[lvl];
    n = lvl_off[lvl + 1] - off;
    nsmooth = mg_get_nsmooth(d, lvl);
    sigma = mg_get_sigma(d, lvl);
    diag = mg_get_diagonal(d, lvl);

    // u = sigma*D*rhs
    for (j = 0; j < n; j++)
      u[off + j] = sigma * r[off + j] / diag[j];

    // Gs = G*u
    mg_operator(Gs + off, u + off, lvl, d, buf);

    // r = rhs - Gs
    for (j = 0; j < n; j++)
      r[off + j] = r[off + j] - Gs[off + j];

    for (i = 0; i < nsmooth; i++) {
      sigma = sigma + 0.066666 / nsmooth;
      // s = sigma*D*r, u = u + s
      for (j = 0; j < n; j++) {
        s[off + j] = sigma * r[off + j] / diag[j];
        u[off + j] += s[off + j];
      }

      // Gs = G*s
      mg_operator(Gs + off, s + off, lvl, d, buf);

      // r = r - Gs
      for (j = 0; j < n; j++)
        r[off + j] = r[off + j] - Gs[off + j];
    }

    // Restrict to coarser level
    mg_restrict(r + off, lvl, d, buf);
  }

  // Solve at the coarsest level
  off = lvl_off[nlevels - 1];
  n = lvl_off[nlevels] - off;

  if (n == 1) {
    diag = mg_get_diagonal(d, nlevels - 1);
    if (fabs(diag[0]) > sqrt(GENMAP_TOL))
      u[off] = r[off] / diag[0];
    else
      u[off] = 0.0;
    r[off] = u[off];
  }

  GenmapScalar over = 1.33333;
  for (lvl = nlevels - 2; lvl >= 0; lvl--) {
    off = lvl_off[lvl];
    // Je = J*e
    mg_interpolate(r + off, lvl, d, buf);

    // u = u + over * Je
    n = lvl_off[lvl + 1] - off;
    for (j = 0; j < n; j++)
      r[off + j] = over * r[off + j] + u[off + j];
  }

  // avoid this
  for (i = 0; i < lvl_off[1]; i++)
    u1[i] = r[i];

  GenmapFree(s);
  GenmapFree(Gs);
  GenmapFree(r);
  GenmapFree(u);
}
