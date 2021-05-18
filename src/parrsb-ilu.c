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

/* Do a binary search */
static int find_gid(ulong id, ulong *gid, uint n) {
  uint i;
  for (i = 0; i < n; i++)
    if (gid[i] == id)
      return 1;

  return 0;
}

static struct array *find_outbound_rows(struct csr_mat_ *M, struct comm *c,
                                        buffer *buf) {
  uint rn = M->rn;
  ulong *col = M->col;

  struct array rows_aux;
  array_init(struct row_info, &rows_aux, 10);

  struct row_info ri;
  ri.dest = c->id;

  uint i, j;
  for (i = 0; i < rn; i++) {
    ulong id = M->row_id[i];
    for (j = M->row_off[i]; j < M->row_off[i + 1] && col[j] <= id; j++) {
      ri.gid = col[j];
      ri.proc = col[j] % c->np;
      ri.owner = find_gid(col[j], M->row_id, rn);
      array_cat(struct row_info, &rows_aux, &ri, 1);
    }
  }

  /* Remove local duplicates in rows_aux */
  struct array *rows = tcalloc(struct array, 1);
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
      if (ptr[i].dest != c->id)
        array_cat(struct row_info, rows, &ptr[i], 1);
    }
  }

  array_free(&rows_aux);
  crystal_free(&cr);

  return rows;
}

static int schedule_outbound_rows(uint **send_off_, ulong **send_ids_,
                                  uint **send_procs_, struct array *rows,
                                  uint nlevels, unsigned int *level_off,
                                  ulong *row_id, buffer *buf, int rank) {
  sarray_sort_2(struct row_info, rows->ptr, rows->n, gid, 1, dest, 0, buf);

  uint *send_off = *send_off_ = tcalloc(uint, nlevels + 1);
  uint off = 0;
  send_off[nlevels] = off;

  struct row_info *ptr = rows->ptr;
  uint n = rows->n;

  uint j;
  sint i, e;
  for (i = nlevels - 1; i >= 0; i--) {
    e = level_off[i] - 1;
    if (e >= 0)
      for (; off < n && ptr[off].gid <= row_id[e]; off++)
        ;
    send_off[i] = off;
  }
  assert(send_off[0] == n);

  ulong *send_ids = *send_ids_ = tcalloc(ulong, n);
  uint *send_procs = *send_procs_ = tcalloc(uint, n);

  for (i = 0; i < n; i++) {
    send_ids[i] = ptr[i].gid;
    send_procs[i] = ptr[i].dest;
  }

  return 0;
}

static int schedule_incoming_rows(uint **send_off_, ulong **send_ids_,
                                  uint **send_procs_, struct array *rows,
                                  uint nlevels, unsigned int *level_off,
                                  ulong *row_id, buffer *buf, int rank) {
  sarray_sort_2(struct row_info, rows->ptr, rows->n, gid, 1, dest, 0, buf);

  uint *send_off = *send_off_ = tcalloc(uint, nlevels + 1);
  uint off = 0;
  send_off[nlevels] = off;

  struct row_info *ptr = rows->ptr;
  uint n = rows->n;

  uint j;
  sint i, e;
  for (i = nlevels - 1; i >= 0; i--) {
    e = level_off[i] - 1;
    if (e >= 0)
      for (; off < n && ptr[off].gid <= row_id[e]; off++)
        ;
    send_off[i] = off;
  }
  assert(send_off[0] == n);

  ulong *send_ids = *send_ids_ = tcalloc(ulong, n);
  uint *send_procs = *send_procs_ = tcalloc(uint, n);

  for (i = 0; i < n; i++) {
    send_ids[i] = ptr[i].gid;
    send_procs[i] = ptr[i].dest;
  }

  return 0;
}

static int send_outbound_rows(struct array *entries, struct csr_mat_ *M,
                              int level, uint *send_off, ulong *send_ids,
                              uint *send_procs, struct comm *c) {
  uint s = send_off[level + 1];
  uint e = send_off[level];

  struct csr_entry t;
  uint i, j, idx = 0;
  for (i = s; i < e; i++) {
    ulong id = send_ids[i];
    for (; idx < M->rn && M->row_id[idx] < id; idx++)
      ;
    assert(M->row_id[idx] == id);

    t.r = id;
    t.proc = send_procs[i];
    for (j = M->row_off[idx]; j < M->row_off[idx + 1]; j++) {
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
  uint i, ii, j, k, kk;

  uint start = lvl_off[lvl + 1];
  uint end = lvl_off[lvl];

  for (ii = start; ii < end; ii++) { /* Go over number of rows */
    i = M->row_id[ii];
    /* TODO: Only loop over non-zeros */
    for (k = 1; k < i; k++) {
      if (csr_mat_get_global(&a_kk, NULL, M, k, k) != 0)
        assert(csr_mat_get_global(&a_kk, NULL, N, k, k) == 0);
      if (fabs(a_kk) < 1e-12)
        continue;

      assert(csr_mat_get_local(&a_ik, &kk, M, ii + 1, k) == 0);
      if (fabs(a_ik) < 1e-12)
        continue;

      /* a_ik = a_ik / a_kk */
      a_ik = v[kk] = v[kk] / a_kk;

      /* a_ij = a_ij - a_ik * a_kj */
      for (kk++; kk < off[ii + 1]; kk++) { /* Go over the columns */
        j = col[kk];
        if (csr_mat_get_global(&a_kj, NULL, M, k, j) != 0)
          assert(csr_mat_get_global(&a_kj, NULL, N, k, j) == 0);
        v[kk] = v[kk] - a_ik * a_kj;
      }
    }
  }

  return 0;
}

void ilu0_setup(uint **send_off, ulong **send_ids, uint **send_procs,
                int nlevels, unsigned int *level_off, struct csr_mat_ *M,
                struct comm *c, buffer *buf) {
  struct array *rows = find_outbound_rows(M, c, buf);
  schedule_outbound_rows(send_off, send_ids, send_procs, rows, nlevels,
                         level_off, M->row_id, buf, c->id);

  array_free(rows);
  free(rows);
}

static int ilu0_aux(struct csr_mat_ *M, unsigned int nlevels,
                    unsigned int *level_off, uint *send_off, ulong *send_ids,
                    uint *send_procs, struct comm *c, buffer *buf) {
  struct csr_mat_ *N = NULL;
  sint i;
  for (i = nlevels - 1; i >= 1; i--) {
    /* Do the local ilu */
    if (c->id == 0)
      printf("Doing Level = %d\n", i);
    ilu0_level(M, N, i, level_off);

    /* Send/recv rows from other processors */
    struct array entries;
    array_init(struct csr_entry, &entries, 10);
    send_outbound_rows(&entries, M, i, send_off, send_ids, send_procs, c);
    if (N != NULL)
      csr_mat_free(N);
    csr_mat_setup(&N, &entries, c, buf);
    array_free(&entries);
  }

  return 0;
}

void parrsb_ilu0(unsigned int nlevels, unsigned int *level_off,
                 struct csr_mat_ *M, MPI_Comm comm) {
  struct comm c;
  comm_init(&c, comm);

  buffer buf;
  buffer_init(&buf, 1024);

  uint *send_off = NULL;
  ulong *send_ids = NULL;
  uint *send_procs = NULL;
  ilu0_setup(&send_off, &send_ids, &send_procs, nlevels, level_off, M, &c,
             &buf);

  ilu0_aux(M, nlevels, level_off, send_off, send_ids, send_procs, &c, &buf);

  if (send_off != NULL)
    free(send_off);
  if (send_ids != NULL)
    free(send_ids);
  if (send_procs != NULL)
    free(send_procs);

  buffer_free(&buf);
  comm_free(&c);
}
