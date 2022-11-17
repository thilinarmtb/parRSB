#include "sort.h"

extern unsigned get_proc_bin(uint id, uint np);

static void sfree(void *p, const char *file, unsigned line) {
  if (p)
    free(p);
}
#define tfree(p) sfree(p, __FILE__, __LINE__)

static void reorder_dofs_local(uint *levels, uint s, uint e, uint nv, uint lvl,
                               const long long *ids, struct comm *c,
                               buffer *bfr, sint *wrk) {
  if (e <= s + 1) {
    for (uint i = s * nv; i < (s + 1) * nv; i++)
      if (ids[i] > 0 && levels[i] == 0)
        levels[i] = lvl;
    return;
  }

  uint size = (e - s) * nv;
  struct gs_data *gsh = gs_setup(&ids[s * nv], size, c, 0, gs_pairwise, 0);

  // Identify the dofs on the interface.
  uint mid = (s + e) / 2;
  for (uint i = s * nv; i < mid * nv; i++)
    wrk[i] = 0;
  for (uint i = mid * nv; i < e * nv; i++)
    wrk[i] = 1;

  gs(&wrk[s * nv], gs_int, gs_add, 0, gsh, bfr);

  for (uint i = mid * nv; i < e * nv; i++)
    wrk[i] = 0;

  gs(&wrk[s * nv], gs_int, gs_add, 0, gsh, bfr);

  for (uint i = s * nv; i < e * nv; i++) {
    if (wrk[i] > 0 && ids[i] > 0 && levels[i] == 0)
      levels[i] = lvl;
  }
  gs_free(gsh);

  // Recursively, go down numbering the other levels.
  reorder_dofs_local(levels, s, mid, nv, lvl - 1, ids, c, bfr, wrk);
  reorder_dofs_local(levels, mid, e, nv, lvl - 1, ids, c, bfr, wrk);
  return;
}

void parrsb_reorder_dofs(long long *nid, unsigned n, unsigned nv,
                         const long long *ids, const MPI_Comm comm) {
  struct comm ci;
  comm_init(&ci, comm);

  buffer bfr;
  buffer_init(&bfr, n);

  uint *levels = tcalloc(uint, n);
  for (uint i = 0; i < n; i++)
    levels[i] = nid[i] = 0;

  // Let's identify the levels of the dofs as we go down the RSB partitioning
  // tree.
  struct comm c;
  comm_split(&ci, n > 0, ci.id, &c);
  if (n > 0) {
    sint *wrk = tcalloc(sint, n);
    unsigned lvl = 1e6;
    while (c.np > 1) {
      struct gs_data *gsh = gs_setup(ids, n, &c, 0, gs_pairwise, 0);

      unsigned bin = get_proc_bin(c.id, c.np);
      for (uint i = 0; i < n; i++)
        wrk[i] = bin;
      gs(wrk, gs_int, gs_add, 0, gsh, &bfr);

      if (bin) {
        for (uint i = 0; i < n; i++)
          wrk[i] = 0;
      }
      gs(wrk, gs_int, gs_add, 0, gsh, &bfr);

      for (uint i = 0; i < n; i++) {
        if (wrk[i] > 0 && ids[i] > 0 && levels[i] == 0)
          levels[i] = lvl;
      }

      gs_free(gsh);

      struct comm t;
      comm_split(&c, bin, c.id, &t);
      comm_free(&c);
      comm_dup(&c, &t);
      comm_free(&t);

      lvl--;
    }
    // Now identify the levels of the local dofs.
    uint ne = n / nv;
    reorder_dofs_local(levels, 0, ne, nv, lvl, ids, &c, &bfr, wrk);

    tfree(wrk);
  }
  comm_free(&c);

  // Sanity check
  for (uint i = 0; i < n; i++) {
    if (ids[i] > 0 && levels[i] == 0) {
      printf("Some dofs are not marked !\n");
      fflush(stdout);
      exit(1);
    }
  }

  // Number ids based on the level. We send the ids to % p to make sure the
  // numbering is consistent and continuous.
  struct dof_t {
    ulong nid, id;
    uint level, seq, p;
  };

  struct array dofs;
  array_init(struct dof_t, &dofs, n);

  struct dof_t dof;
  for (uint i = 0; i < n; i++) {
    if (ids[i]) {
      dof.id = ids[i], dof.p = dof.id % ci.np;
      dof.seq = i, dof.level = levels[i];
      array_cat(struct dof_t, &dofs, &dof, 1);
    }
  }
  tfree(levels);

  sarray_sort(struct dof_t, dofs.ptr, dofs.n, level, 0, &bfr);

  sint min = INT_MAX, max = -INT_MAX;
  if (dofs.n > 0) {
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    min = pd[0].level, max = pd[dofs.n - 1].level;
  }

  slong wrk[2][1];
  comm_allreduce(&ci, gs_int, gs_max, &max, 1, wrk);
  comm_allreduce(&ci, gs_int, gs_min, &min, 1, wrk);
  // FIXME: +1 is wrong if problem size is 0.
  sint nlvls = max - min + 1;

  struct crystal cr;
  crystal_init(&cr, &ci);

  sarray_transfer(struct dof_t, &dofs, p, 1, &cr);

  // Sanity check.
  sarray_sort_2(struct dof_t, dofs.ptr, dofs.n, id, 1, level, 0, &bfr);
  if (dofs.n > 0) {
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    uint i = 1, j = 0;
    while (i < dofs.n) {
      while (i < dofs.n && pd[j].id == pd[i].id)
        assert(pd[j].level == pd[i].level), i++;
      j = i++;
    }
  }

  sarray_sort_2(struct dof_t, dofs.ptr, dofs.n, level, 0, id, 1, &bfr);

  struct dof_t *pd = (struct dof_t *)dofs.ptr;
  uint idx = 0;
  slong ng = 0;
  for (uint lvl = 0; lvl < nlvls; lvl++) {
    sint l = INT_MAX;
    if (idx < dofs.n)
      l = pd[idx].level;
    comm_allreduce(&ci, gs_int, gs_min, &l, 1, wrk);

    ulong id = 0;
    uint idx1 = idx, k = 0, k1 = 0;
    for (; idx1 < dofs.n && pd[idx1].level == l; idx1++) {
      if (pd[idx1].id != id)
        id = pd[idx1].id, k++;
    }

    slong out[2][1], in = k;
    comm_scan(out, &ci, gs_long, gs_add, &in, 1, wrk);
    slong s = out[0][0];

    id = 0;
    for (; idx < idx1; idx++) {
      if (pd[idx].id != id)
        id = pd[idx].id, k1++;
      pd[idx].nid = ng + s + k1;
    }
    assert(k == k1);

    ng += out[1][0];
  }

  sarray_transfer(struct dof_t, &dofs, p, 0, &cr);

  // Sanity check.
  sarray_sort_2(struct dof_t, dofs.ptr, dofs.n, id, 1, nid, 1, &bfr);
  if (dofs.n > 0) {
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    uint i = 1, j = 0;
    while (i < dofs.n) {
      while (i < dofs.n && pd[j].id == pd[i].id)
        assert(pd[j].nid == pd[i].nid), i++;
      j = i, i++;
    }
  }

  sarray_sort(struct dof_t, dofs.ptr, dofs.n, seq, 0, &bfr);

  if (dofs.n > 0) {
    struct dof_t *pd = (struct dof_t *)dofs.ptr;
    for (uint i = 0; i < dofs.n; i++)
      nid[pd[i].seq] = pd[i].nid;
  }

  array_free(&dofs), buffer_free(&bfr);
  crystal_free(&cr), comm_free(&ci);
}

#define fparrsb_reorder_dofs                                                   \
  FORTRAN_UNPREFIXED(fparrsb_order_dofs, FPARRSB_ORDER_DOFS)
void fparrsb_reorder_dofs(long long *nid, int *n, int *nv, long long *ids,
                          MPI_Fint *comm, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*comm);
  parrsb_reorder_dofs(nid, *n, *nv, ids, c);
  *err = 0;
}

#undef tfree
