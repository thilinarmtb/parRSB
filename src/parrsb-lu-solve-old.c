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

struct rhs_entry {
  ulong r;
  uint proc;
  GenmapScalar v;
  uint old;
};

/* TODO: Do a binary search */
static int find_gid(ulong id, ulong *gid, uint n) {
  uint i;
  for (i = 0; i < n; i++)
    if (gid[i] == id)
      return 1;

  return 0;
}

static struct array *find_outbound_rows(struct array *rows, int p, int fw,
                                        int level, unsigned int *level_off,
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

  sint i, j;
  if (fw == 1)
    for (i = s; i < e; i++) {
      ulong id = M->row_id[i];
      for (j = M->row_off[i]; j < M->row_off[i + 1] && col[j] <= id; j++) {
        ri.gid = col[j];
        ri.proc = col[j] % c->np;
        ri.owner = find_gid(col[j], M->row_id, rn);
        array_cat(struct row_info, &rows_aux, &ri, 1);
        // printf("p = %d level = %d row = %ld id = %ld\n", p, level, id,
        // col[j]);
      }
    }
  else
    for (i = s; i < e; i++) {
      ulong id = M->row_id[i];
      for (j = M->row_off[i + 1] - 1; j >= (int)M->row_off[i] && col[j] >= id;
           j--) {
        ri.gid = col[j];
        ri.proc = col[j] % c->np;
        ri.owner = find_gid(col[j], M->row_id, rn);
        array_cat(struct row_info, &rows_aux, &ri, 1);
        // printf("p = %d level = %d row = %ld id = %ld\n", p, level, id,
        // col[j]);
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
#if 0
    if (level == 1 && p == 0)
      printf("cur_gid1=%ld\n", cur_gid);
#endif

    for (i = 1; i < n; i++) {
      if (ptr[i].gid != cur_gid) {
        assert(cur_owner >= 0);
        for (j = cur_idx; j < i; j++)
          ptr[j].proc = cur_owner;
        cur_gid = ptr[i].gid;
        cur_idx = i;
        cur_owner = -1;
#if 0
        if (level == 1 && p == 0)
          printf("cur_gid1=%ld\n", cur_gid);
#endif
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
      if (ptr[i].dest == p) {
        array_cat(struct row_info, rows, &ptr[i], 1);
#if 0
        if (level == 1 && p == 0)
          printf("cur_gid2=%ld\n", ptr[i].gid);
#endif
      }
    }
  }

  array_free(&rows_aux);
  crystal_free(&cr);

  return rows;
}

static int send_outbound_rows(struct array *entries, struct array *rows,
                              struct csr_mat_ *M, double *b, struct comm *c,
                              int level, int p, buffer *buf) {
  sarray_sort_2(struct row_info, rows->ptr, rows->n, gid, 1, dest, 0, buf);
  struct row_info *ptr = rows->ptr;
  uint n = rows->n;

  array_init(struct rhs_entry, entries, 10);

  struct rhs_entry t;
  uint i, j, s;
  for (i = s = 0; i < n; i++) {
    ulong gid = ptr[i].gid;
    while (s < M->rn && M->row_id[s] < gid)
      s++;
    assert(M->row_id[s] == gid);

    t.r = gid;
    t.proc = ptr[i].dest;
    t.v = b[s];
    array_cat(struct rhs_entry, entries, &t, 1);
#if 0
    if (level == 1 && p == 0)
      printf("gid1=%ld %d\n", ptr[i].gid, ptr[i].dest);
#endif
  }

  struct crystal cr;
  crystal_init(&cr, c);
  sarray_transfer(struct rhs_entry, entries, proc, 1, &cr);
  crystal_free(&cr);

#if 0
  struct rhs_entry *pr = entries->ptr;
  for (i = 0; i < entries->n; i++)
    if (level == 1 && p == 0)
      printf("%d %d %u gid11=%ld\n", c->id, i, entries->n, pr[i].r);
#endif

  return 0;
}

static void concat_outbound_rows(struct array *all, struct rhs_entry *ptr,
                                 uint n, int level, int p, struct comm *c,
                                 buffer *buf) {
  struct rhs_entry *ptr0 = all->ptr;
  uint i;
  for (i = 0; i < all->n; i++)
    ptr0[i].old = 1;
  for (i = 0; i < n; i++) {
    ptr[i].old = 0;
    // printf("p = %d level = %d r = %ld, v = %lf\n", p, level, ptr[i].r,
    // ptr[i].v);
  }

  struct array temp;
  array_init(struct rhs_entry, &temp, 10);
#if 0
  if (level == 1 && p == 0)
    printf("%d temp.n=%u all->n=%u n=%u\n", c->id, temp.n, all->n, n);
#endif

  array_cat(struct rhs_entry, &temp, ptr, n);
  array_cat(struct rhs_entry, &temp, all->ptr, all->n);
  sarray_sort_2(struct rhs_entry, temp.ptr, temp.n, r, 1, old, 0, buf);

  array_free(all);
  array_init(struct rhs_entry, all, 10);

  if (temp.n > 0) {
    struct rhs_entry *ptr = temp.ptr;
    ulong gid = ptr[0].r;
    array_cat(struct rhs_entry, all, &ptr[0], 1);
#if 0
    if (level == 1 && p == 0)
      printf("gid2=%ld\n", gid);
#endif

    for (i = 1; i < temp.n; i++)
      if (ptr[i].r != gid) {
        array_cat(struct rhs_entry, all, &ptr[i], 1);
        gid = ptr[i].r;
#if 0
        if (level == 1 && p == 0)
          printf("gid2=%ld\n", gid);
#endif
      }
  }

  array_free(&temp);
}

static double get_rhs(ulong gid, uint n, ulong *ids, double *y,
                      struct array *entries, int level, int rank) {
  uint i;
  for (i = 0; i < n; i++)
    if (ids[i] == gid)
      return y[i];

  struct rhs_entry *ptr = entries->ptr;
  for (i = 0; i < entries->n; i++)
    if (ptr[i].r == gid)
      return ptr[i].v;

  printf("%d %d not found: %lu\n", level, rank, gid);

  return 0.0;
}

static int fw_solve_level(double *y, struct csr_mat_ *A, double *b,
                          struct array *entries, int level,
                          unsigned int *level_off, int p) {
  /* Forward substitution with L - Lower matrix (diagonal = all 1s)
   * L(Ux) = L(y) = b */
  sint i, j;
  for (i = level_off[level + 1]; i < level_off[level]; i++) {
    y[i] = b[i];
    /* Remove the dot product: L[i, j] * y[j], j < i */
    for (j = A->row_off[i]; j < A->row_off[i + 1] && A->col[j] < A->row_id[i];
         j++) {
      double x = get_rhs(A->col[j], A->rn, A->row_id, y, entries, level, p);
      y[i] -= A->v[j] * x;
      // printf("level = %d row = %ld p = %d b[i] = %lf, v=%lf x = %lf y[i] =
      // %lf\n",
      //        level, A->row_id[i], p, b[i], A->v[j], x, y[i]);
    }
    // printf("forward level = %d row = %ld p = %d b[i] = %lf, y[i] = %lf\n",
    //        level, A->row_id[i], p, b[i], y[i]);
  }

  return 0;
}

static int bw_solve_level(double *x, struct csr_mat_ *A, double *y,
                          struct array *entries, int level,
                          unsigned int *level_off, int p) {
  /* Back substitution with U, Ux = y */
  assert(level_off[level] <= A->rn && "Ooops !");
  assert(level_off[level + 1] <= A->rn && "Ooops !");

  sint i, j;
  for (i = level_off[level] - 1; i >= (int)level_off[level + 1]; i--) {
    x[i] = y[i];
    /* Remove the dot product: L[i, j] * x[j], j > i */
    for (j = A->row_off[i + 1] - 1;
         j >= A->row_off[i] && A->col[j] > A->row_id[i]; j--)
      x[i] -=
          A->v[j] * get_rhs(A->col[j], A->rn, A->row_id, x, entries, level, p);
    if (fabs(A->v[j] - 0) > 1e-12)
      x[i] /= A->v[j];
    // printf("backward level = %d row = %ld p = %d y[i] = %lf, x[i] = %lf\n",
    //        level, A->row_id[i], p, y[i], x[i]);
  }

  return 0;
}

static int lu_solve_aux(double *x, struct csr_mat_ *A, double *b, int nlevels,
                        unsigned int *level_off, MPI_Comm *comms, buffer *buf) {
  uint n = A->rn;
  double *y = tcalloc(double, n);

  sint i = nlevels - 1;
  fw_solve_level(y, A, b, NULL, i, level_off, 0);

  struct array all;
  array_init(struct rhs_entry, &all, 10);

  for (i = nlevels - 2; i >= 0; i--) {
    struct comm c;
    comm_init(&c, comms[i]);

    /* Do LU solve in serial among the processors in the current level */
    int p;
    for (p = 0; p < c.np; p++) {
      /* Receive rows from other processors */
      struct array entries, rows;
      find_outbound_rows(&rows, p, 1, i, level_off, A, &c, buf);
      send_outbound_rows(&entries, &rows, A, y, &c, i, p, buf);
      concat_outbound_rows(&all, entries.ptr, entries.n, i, p, &c, buf);

      if (c.id == p)
        fw_solve_level(y, A, b, &all, i, level_off, p);

      array_free(&entries);
      array_free(&rows);
    }

    comm_free(&c);
  }

  array_free(&all);
  array_init(struct rhs_entry, &all, 10);

  for (i = 0; i <= nlevels - 2; i++) {
    struct comm c;
    comm_init(&c, comms[i]);

    /* Do LU solve in serial among the processors in the current level */
    int p;
    for (p = c.np - 1; p >= 0; p--) {
      /* Receive rows from other processors */
      struct array entries, rows;
      find_outbound_rows(&rows, p, 0, i, level_off, A, &c, buf);
      send_outbound_rows(&entries, &rows, A, x, &c, i, p, buf);
      concat_outbound_rows(&all, entries.ptr, entries.n, i, p, &c, buf);

      if (c.id == p)
        bw_solve_level(x, A, y, &all, i, level_off, p);

      array_free(&entries);
      array_free(&rows);
    }

    comm_free(&c);
  }

  bw_solve_level(x, A, y, &all, i, level_off, 0);

  array_free(&all);
  free(y);

  return 0;
}

void parrsb_lu_solve(double *x, struct csr_mat_ *M, double *b,
                     unsigned int nlevels, unsigned int *level_off,
                     MPI_Comm *comms) {
  uint n = M->rn;
  buffer buf;
  buffer_init(&buf, sizeof(double) * n);
  lu_solve_aux(x, M, b, nlevels, level_off, comms, &buf);
  buffer_free(&buf);
}
