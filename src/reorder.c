#include "sort.h"

extern unsigned get_proc_bin(uint id, uint np);

#define MIN(a, b) ((a) < (b) ? (a) : (b))

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

  sarray_sort_2(struct dof_t, dofs.ptr, dofs.n, level, 0, id, 1, &bfr);

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

int FACE_VERTICES[6][4] = {{0, 2, 4, 6}, {1, 3, 5, 7}, {0, 1, 4, 5},
                           {2, 3, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7}};

void id_elements_with_unconnected_faces(int *frontier, uint ne, unsigned nv,
                                        const long long *vids, buffer *bfr) {
  if (nv != 8) {
    fprintf(stderr, "Only 3D hexes (nv = 8) are supported !\n");
    fflush(stderr);
    exit(1);
  }
  unsigned nf = 6;
  unsigned nvf = nv / 2;

  // If an element has at least one face on the boundary (i.e., not connected
  // to any other face, the element is on the frontier. Refer src/con.c for
  // symmetric numbering of vertices and faces. A face can only have only one
  // other face neighbor.
  struct face_t {
    uint seq, connected;
    ulong vtx[4];
  };

  struct array faces;
  array_init(struct face_t, &faces, 6 * ne + 1);

  struct face_t face = {.connected = 0};
  ulong tmp;
  for (uint e = 0; e < ne; e++) {
    face.seq = e;
    for (unsigned f = 0; f < nf; f++) {
      // set face vertices
      for (unsigned v = 0; v < nvf; v++)
        face.vtx[v] = vids[e * nv + FACE_VERTICES[f][v]];
      // sort face vertices
      for (unsigned v = 0; v < nvf - 1; v++) {
        for (unsigned v1 = v + 1; v1 < nvf; v1++) {
          if (face.vtx[v1] < face.vtx[v])
            tmp = face.vtx[v], face.vtx[v] = face.vtx[v1], face.vtx[v1] = tmp;
        }
      }
      array_cat(struct face_t, &faces, &face, 1);
    }
  }

  sarray_sort_3(struct face_t, faces.ptr, faces.n, vtx[0], 1, vtx[1], 1, vtx[2],
                1, bfr);
  struct face_t *pf = (struct face_t *)faces.ptr;
  uint i = 0;
  while (i + 1 < faces.n) {
    if (pf[i].vtx[0] == pf[i + 1].vtx[0] && pf[i].vtx[0] == pf[i + 1].vtx[1] &&
        pf[i].vtx[2] == pf[i + 1].vtx[2]) {
      pf[i].connected = pf[i + 1].connected = 1, i++;
    }
    i++;
  }

  sarray_sort(struct face_t, faces.ptr, faces.n, seq, 0, bfr);
  pf = (struct face_t *)faces.ptr;
  for (uint e = 0; e < ne; e++) {
    frontier[e] = 0;
    for (unsigned f = 0; f < nf; f++) {
      if (!pf[e * nf + f].connected) {
        frontier[e] = 1;
        break;
      }
    }
  }

  array_free(&faces);
}

struct elem_t {
  ulong eid, vid[8];
  double xyz[8 * 3], mat[8 * 8], mask[8];
  uint p;
};

void find_vtx_frontiers(int *frontier, long long *vtx, struct array *uelems,
                        uint ne, const long long *eids, unsigned nv,
                        struct comm *ci, buffer *bfr) {
  struct id_t {
    ulong id;
  };

  struct array ids;
  array_init(struct id_t, &ids, ne + 1);

  struct id_t id;
  for (uint i = 0; i < ne; i++) {
    id.id = eids[i];
    array_cat(struct id_t, &ids, &id, 1);
  }

  sarray_sort(struct id_t, ids.ptr, ids.n, id, 1, bfr);

  // Mark everything as in frontier first.
  uint nu = uelems->n;
  for (uint i = 0; i < nu * nv; i++)
    frontier[i] = 0;

  // Unmark them if the element was part of input.
  struct elem_t *pu = (struct elem_t *)uelems->ptr;
  struct id_t *pi = (struct id_t *)ids.ptr;
  uint i = 0, j = 0;
  while (i < ne) {
    for (; j < nu && pu[j].eid < pi[i].id; j++)
      ;
    for (; j < nu && pu[j].eid == pi[i].id; j++) {
      for (unsigned v = 0; v < nv; v++)
        frontier[j * nv + v] = 1;
    }
    i++;
  }

  if (uelems->n > 0) {
    struct elem_t *pe = (struct elem_t *)uelems->ptr;
    for (uint i = 0; i < uelems->n; i++) {
      for (unsigned v = 0; v < nv; v++)
        vtx[i * nv + v] = pe[i].vid[v];
    }
    struct comm c;
    comm_split(ci, ci->id, ci->id, &c);

    struct gs_data *gsh = gs_setup(vtx, uelems->n * nv, &c, 0, gs_pairwise, 0);
    gs(frontier, gs_int, gs_max, 0, gsh, bfr);

    gs_free(gsh), comm_free(&c);
  }

  array_free(&ids);
}

void parrsb_fetch_nbrs(unsigned *nei, long long *eids, unsigned nv,
                       long long *vids, double *xyz, int *frontier,
                       double *mask, double *mat, const MPI_Comm comm,
                       unsigned maxne) {
  size_t ne = *nei, size = nv;
  size *= ne;

  struct vtx_t {
    ulong id;
    uint p, o, seq;
  };

  struct array vtxs;
  array_init(struct vtx_t, &vtxs, size);

  struct comm ci;
  comm_init(&ci, comm);

  struct vtx_t vt;
  for (uint e = 0; e < ne; e++) {
    for (unsigned v = 0; v < nv; v++) {
      vt.id = vids[e * nv + v], vt.o = ci.id, vt.p = vt.id % ci.np, vt.seq = e;
      array_cat(struct vtx_t, &vtxs, &vt, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, &ci);

  sarray_transfer(struct vtx_t, &vtxs, p, 1, &cr);

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort(struct vtx_t, vtxs.ptr, vtxs.n, id, 1, &bfr);

  struct array vtx2p;
  array_init(struct vtx_t, &vtx2p, vtxs.n);

  struct vtx_t *pv = (struct vtx_t *)vtxs.ptr;
  uint s = 0;
  while (s < vtxs.n) {
    uint e = s + 1;
    while (e < vtxs.n && pv[s].id == pv[e].id)
      e++;
    for (uint i = s; i < e; i++) {
      vt = pv[i];
      for (uint j = s; j < e; j++) {
        vt.o = pv[j].o;
        array_cat(struct vtx_t, &vtx2p, &vt, 1);
      }
    }
    s = e;
  }
  array_free(&vtxs);

  sarray_transfer(struct vtx_t, &vtx2p, p, 0, &cr);
  sarray_sort_2(struct vtx_t, vtx2p.ptr, vtx2p.n, seq, 0, o, 0, &bfr);

  uint *offs = tcalloc(uint, ne + 1);
  uint *proc = tcalloc(uint, ne * 4), nproc = ne * 4 + 1;
  pv = (struct vtx_t *)vtx2p.ptr, s = 0;
  while (s < vtx2p.n) {
    uint seq = pv[s].seq, e = s + 1;
    while (e < vtx2p.n && seq == pv[e].seq)
      e++;

    // FIXME: Make work a dynamic array.
    uint work[1024], n = 1;
    work[0] = pv[s].o;
    for (uint i = s + 1; i < e; i++) {
      if (work[n - 1] != pv[i].o)
        work[n] = pv[i].o, n++;
    }

    offs[seq + 1] = offs[seq] + n;
    if (offs[seq + 1] >= nproc) {
      nproc = offs[seq + 1] + nproc / 2 + 1;
      proc = trealloc(uint, proc, nproc);
    }
    for (uint i = 0, j = offs[seq]; i < n; i++)
      proc[j + i] = work[i];
    s = e;
  }
  array_free(&vtx2p);

  struct array elems;
  array_init(struct elem_t, &elems, size);

  unsigned nd = (nv == 8) ? 3 : 2;
  struct elem_t et;
  for (uint e = 0; e < ne; e++) {
    et.eid = eids[e];
    for (unsigned v = 0; v < nv; v++) {
      et.vid[v] = vids[e * nv + v], et.mask[v] = mask[e * nv + v];
      for (unsigned d = 0; d < nd; d++)
        et.xyz[v * nd + d] = xyz[e * nv * nd + v * nd + d];
      for (unsigned u = 0; u < nv; u++)
        et.mat[v * nv + u] = mat[e * nv * nv + v * nv + u];
    }
    for (uint s = offs[e]; s < offs[e + 1]; s++) {
      et.p = proc[s];
      array_cat(struct elem_t, &elems, &et, 1);
    }
  }
  tfree(offs), tfree(proc);

  sarray_transfer(struct elem_t, &elems, p, 1, &cr);
  crystal_free(&cr);

  // Get rid of the duplicates.
  struct array uelems;
  array_init(struct elem_t, &uelems, elems.n / 2 + 1);

  sarray_sort(struct elem_t, elems.ptr, elems.n, eid, 1, &bfr);
  if (elems.n > 0) {
    struct elem_t *pe = (struct elem_t *)elems.ptr;

    array_cat(struct elem_t, &uelems, &pe[0], 1);
    uint i = 0, j = 1;
    while (j < elems.n) {
      if (pe[j].eid != pe[i].eid) {
        array_cat(struct elem_t, &uelems, &pe[j], 1);
        i = j;
      }
      j++;
    }
  }
  array_free(&elems);

  sint err = (uelems.n > maxne), wrk;
  comm_allreduce(&ci, gs_int, gs_add, &err, 1, &wrk);
  if (err > 0) {
    if (ci.id == 0) {
      fprintf(stderr, "maxne = %u is not large enough !\n", maxne);
      fflush(stderr);
    }
    buffer_free(&bfr), array_free(&uelems), comm_free(&ci);
    exit(1);
  }

  find_vtx_frontiers(frontier, vids, &uelems, ne, eids, nv, &ci, &bfr);

  *nei = uelems.n;
  if (uelems.n > 0) {
    struct elem_t *pe = (struct elem_t *)uelems.ptr;
    for (uint e = 0; e < uelems.n; e++) {
      eids[e] = pe[e].eid;
      for (unsigned v = 0; v < nv; v++) {
        mask[e * nv + v] = pe[e].mask[v];
        for (unsigned d = 0; d < nd; d++)
          xyz[e * nv * nd + v * nd + d] = pe[e].xyz[v * nd + d];
        for (unsigned u = 0; u < nv; u++)
          mat[e * nv * nv + v * nv + u] = pe[e].mat[v * nv + u];
      }
    }
  }
  array_free(&uelems), comm_free(&ci);

  buffer_free(&bfr);
}

#define fparrsb_fetch_nbrs                                                     \
  FORTRAN_UNPREFIXED(fparrsb_fetch_nbrs, FPARRSB_FETCH_NBRS)
void fparrsb_fetch_nbrs(int *nei, long long *eids, int *nv, long long *vids,
                        double *xyz, int *frontier, double *mask, double *mat,
                        MPI_Fint *comm, int *maxne, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*comm);
  unsigned ne = *nei;
  parrsb_fetch_nbrs(&ne, eids, *nv, vids, xyz, frontier, mask, mat, c, *maxne);
  *nei = ne;
  *err = 0;
}

#undef tfree
#undef MIN
