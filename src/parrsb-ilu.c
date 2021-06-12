#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-sort.h>

#include <parRSB.h>

struct row_info {
  ulong gid;
  uint dest;
  uint owner;
  uint proc;
};

struct csr_entry {
  ulong r, c;
  uint proc;
  GenmapScalar v;
};

/* TODO: Do a binary search */
static int find_gid(ulong id, ulong *gid, uint n) {
  uint i;
  for (i = 0; i < n; i++)
    if (gid[i] == id)
      return 1;

  return 0;
}

static struct array *find_outbound_rows(struct array *rows, int p, int level,
                                        unsigned int *level_off,
                                        struct csr_mat_ *M, struct comm *c,
                                        buffer *buf) {
  uint rn = M->rn;
  ulong *col = M->col;

  struct array rows_aux;
  array_init(struct row_info, &rows_aux, 10);

  struct row_info ri;
  ri.dest = c->id;

  uint s = level_off[level + 1];
  uint e = level_off[level];

  uint i, j;
  for (i = s; i < e; i++) {
    ulong id = M->row_id[i];
    for (j = M->row_off[i]; j < M->row_off[i + 1] && col[j] <= id; j++) {
      ri.gid = col[j];
      ri.proc = col[j] % c->np;
      ri.owner = find_gid(col[j], M->row_id, rn);
      array_cat(struct row_info, &rows_aux, &ri, 1);
    }
  }

  /* Remove local duplicates in rows_aux */
  uint n = rows_aux.n;
  array_init(struct row_info, rows, n);
  if (n > 0) {
    sarray_sort(struct row_info, rows_aux.ptr, rows_aux.n, gid, 1, buf);

    struct row_info *ptr = rows_aux.ptr;
    array_cat(struct row_info, rows, &ptr[0], 1);
    slong cur_gid = ptr[0].gid;

    for (i = 1; i < rows_aux.n; i++) {
      if (ptr[i].gid != cur_gid) {
        array_cat(struct row_info, rows, &ptr[i], 1);
        cur_gid = ptr[i].gid;
      }
    }
  }
  array_free(&rows_aux);

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct row_info, rows, proc, 0, &cr);
  n = rows->n;

  if (n > 0) {
    sarray_sort(struct row_info, rows->ptr, n, gid, 1, buf);
    struct row_info *ptr = rows->ptr;

    uint cur_idx = 0;
    slong cur_gid = ptr[0].gid;
    sint cur_owner = ptr[0].owner > 0 ? ptr[0].dest : -1;

    for (i = 1; i < n; i++) {
      if (ptr[i].gid != cur_gid) {
        assert(cur_owner >= 0);
        for (j = cur_idx; j < i; j++)
          ptr[j].proc = cur_owner;
        cur_gid = ptr[i].gid;
        cur_idx = i;
        cur_owner = -1;
      }
      if (ptr[i].owner > 0)
        cur_owner = ptr[i].dest;
    }
    assert(cur_owner >= 0);
    for (j = cur_idx; j < n; j++)
      ptr[j].proc = cur_owner;
  }

  sarray_transfer(struct row_info, rows, proc, 0, &cr);
  n = rows->n;

  array_init(struct row_info, &rows_aux, n);
  array_cat(struct row_info, &rows_aux, rows->ptr, n);
  array_free(rows);
  array_init(struct row_info, rows, 10);

  /* Get rid of sends to own processor */
  if (n > 0) {
    struct row_info *ptr = rows_aux.ptr;
    for (i = 0; i < n; i++) {
      if (ptr[i].dest == p)
        array_cat(struct row_info, rows, &ptr[i], 1);
    }
  }

  array_free(&rows_aux);
  crystal_free(&cr);

  return rows;
}

static int send_outbound_rows(struct array *entries, struct array *rows,
                              struct csr_mat_ *M, struct comm *c, buffer *buf) {
  sarray_sort_2(struct row_info, rows->ptr, rows->n, gid, 1, dest, 0, buf);
  struct row_info *ptr = rows->ptr;
  uint n = rows->n;

  array_init(struct csr_entry, entries, 10);

  struct csr_entry t;
  uint i, j, s;
  for (i = s = 0; i < n; i++) {
    ulong gid = ptr[i].gid;
    while (s < M->rn && M->row_id[s] < gid)
      s++;
    assert(M->row_id[s] == gid);

    t.r = gid;
    t.proc = ptr[i].dest;
    for (j = M->row_off[s]; j < M->row_off[s + 1]; j++) {
      t.c = M->col[j];
      t.v = M->v[j];
      array_cat(struct csr_entry, entries, &t, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct csr_entry, entries, proc, 1, &cr);
  crystal_free(&cr);

  return 0;
}

static int ilu0_level(struct csr_mat_ *M, struct csr_mat_ *N, int lvl,
                      unsigned int *lvl_off) {
  const uint rn = M->rn;
  assert(rn > 1);

  const uint *off = M->row_off;
  double *v = M->v;
  const ulong *col = M->col;

  double a_kk, a_kj, a_ik;
  uint i, ii, j, k, kk, ik, ij;

  uint start = lvl_off[lvl + 1];
  uint end = lvl_off[lvl];

  for (ii = start; ii < end; ii++) { /* Go over number of rows */
    i = M->row_id[ii];
    /* TODO: Only loop over non-zeros */
    for (kk = M->row_off[ii]; kk < M->row_off[ii + 1] && M->col[kk] < i; kk++) {
      k = M->col[kk];
      if (csr_mat_get_global(&a_kk, NULL, M, k, k) != 0)
        assert(csr_mat_get_global(&a_kk, NULL, N, k, k) == 0);
      if (fabs(a_kk) < 1e-12)
        continue;

      assert(csr_mat_get_local(&a_ik, &ik, M, ii + 1, k) == 0);
      if (fabs(a_ik) < 1e-12)
        continue;

      /* a_ik = a_ik / a_kk */
      a_ik = v[ik] = v[ik] / a_kk;

      /* a_ij = a_ij - a_ik * a_kj */
      for (ij = ik + 1; ij < off[ii + 1]; ij++) { /* Go over the columns */
        j = col[ij];
        if (csr_mat_get_global(&a_kj, NULL, M, k, j) != 0)
          assert(csr_mat_get_global(&a_kj, NULL, N, k, j) == 0);
        v[ij] = v[ij] - a_ik * a_kj;
      }
    }
  }

  return 0;
}

static int ilu0_aux(struct csr_mat_ *M, unsigned int nlevels,
                    unsigned int *level_off, MPI_Comm *comms, buffer *buf) {
  struct csr_mat_ *N = NULL;
  sint i = nlevels - 1;
  ilu0_level(M, N, i, level_off);

  for (i = nlevels - 2; i >= 0; i--) {
    struct comm c;
    comm_init(&c, comms[i]);

    /* Do ILU serially among the processors in the current level */
    int p;
    for (p = 0; p < c.np; p++) {
      /* Receive rows from other processors */
      struct array entries, rows;
      find_outbound_rows(&rows, p, i, level_off, M, &c, buf);
      send_outbound_rows(&entries, &rows, M, &c, buf);

      /* Do the local ilu */
      if (c.id == p) {
        csr_mat_setup(&N, &entries, NULL, buf);
        ilu0_level(M, N, i, level_off);
        csr_mat_free(N);
      }

      array_free(&entries);
      array_free(&rows);
    }

    comm_free(&c);
  }

  return 0;
}

void parrsb_ilu0(unsigned int nlevels, unsigned int *level_off, MPI_Comm *comms,
                 struct csr_mat_ *M) {
  buffer buf;
  buffer_init(&buf, 1024);

  ilu0_aux(M, nlevels, level_off, comms, &buf);

  buffer_free(&buf);
}
