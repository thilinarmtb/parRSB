// FIXME: genmap-impl.h is only for GenmapFree
#include <genmap-impl.h>
#include <genmap-multigrid-gs.h>

void mg_setup_aux_gs(struct mg_data_gs *d, slong *wrk, buffer *buf) {
  struct comm *c = d->c;

  uint nelt = d->level_off[1] - d->level_off[0];

  slong in;
  slong out[2][1], bf[2][1];

  uint size0, size1;
  size0 = nelt;

  sint np0, np1;
  np0 = c->np;

  sint j, lvl;
  for (lvl = 1; lvl < d->nlevels; lvl++) {
    in = size0;
    comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
    slong start0 = out[0][0];
    slong ng0 = out[1][0];

    slong ng1 = (ng0 + 1) / 2;
    int is_odd = 0;
    if (ng0 > 2 && (ng0 % 2) == 1) {
      ng1 -= 1;
      is_odd = 1;
    }

    np1 = (ng1 >= np0) ? np0 : ng1;

    size1 = (c->id < np1) ? ng1 / np1 : 0;
    sint nrem = (c->id < np1) ? ng1 - size1 * np1 : 0;
    /* Following should match set_owner in genmap-multigrid-csr-aux.c */
    size1 += (c->id < nrem) ? 1 : 0;

    in = size1;
    comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
    slong start1 = out[0][0];
    assert(ng1 == out[1][0]);

    for (j = 0; j < size0; j++)
      wrk[j] = -((start0 + j + 1) + 1) / 2;
    if (size0 > 0 && is_odd == 1 && (start0 + size0) == ng0)
      wrk[j - 1] += 1;
    for (j = 0; j < size1; j++)
      wrk[size0 + j] = start1 + j + 1;

    d->J[lvl - 1] = gs_setup(wrk, size0 + size1, c, 0, gs_pairwise, 0);

    d->nsmooth[lvl] = 2;
    d->sigma[lvl] = 0.6;
    d->J[lvl] = NULL;

    d->level_off[lvl + 1] = d->level_off[lvl] + size1;

    size0 = size1;
    np0 = np1;
  }

  /* Set J0 -> Current level to Level 0 interpolation, 0th level and
   * last level does not need J0 */
  for (lvl = d->nlevels - 1; lvl > 0; lvl--) {
    size1 = d->level_off[lvl + 1] - d->level_off[lvl];

    in = size1;
    comm_scan(out, c, gs_long, gs_add, &in, 1, bf);
    slong start1 = out[0][0];

    uint cur_off = d->level_off[lvl];

    for (j = 0; j < size1; j++)
      wrk[cur_off + j] = start1 + j + 1;

    // Apply R^T till we get to Level0
    for (j = lvl - 1; j >= 0; j--)
      gs(wrk + d->level_off[j], gs_long, gs_add, 0, d->J[j], buf);

    // Setup J0 now that we know the corresponding ids at level0
    for (j = 0; j < nelt; j++)
      wrk[j] *= -1;
    for (j = 0; j < size1; j++)
      wrk[nelt + j] = wrk[cur_off + j];

    d->J0[lvl] = gs_setup(wrk, nelt + size1, c, 0, gs_pairwise, 0);
  }
}
