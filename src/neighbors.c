#include "sort.h"
#include <err.h>

static void sfree(void *p, const char *file, unsigned line) {
  if (p)
    free(p);
}
#define tfree(p) sfree(p, __FILE__, __LINE__)

struct elem_t {
  ulong eid, vid[8];
  double xyz[8 * 3], mat[8 * 8], mask[8];
  int frontier[8];
  uint p, seq;
};

int find_frontier(struct array *uelems, uint ne, const long long *eids,
                  unsigned nv, struct comm *ci, buffer *bfr) {
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

  // Initialize the frontier to -1. We will copy the frontier value for all the
  // input elements as is and for the new elements, it will be max(frontier
  // value of input elements) + 1.
  uint nu = uelems->n;
  int *frontier = tcalloc(int, nu *nv);
  for (uint i = 0; i < nu * nv; i++)
    frontier[i] = -1;

  struct elem_t *pu = (struct elem_t *)uelems->ptr;
  struct id_t *pi = (struct id_t *)ids.ptr;
  int maxf = -1;
  uint i = 0, j = 0;
  while (i < ne) {
    for (; j < nu && pu[j].eid < pi[i].id; j++)
      ;
    for (; j < nu && pu[j].eid == pi[i].id; j++) {
      for (unsigned v = 0; v < nv; v++) {
        frontier[j * nv + v] = pu[j].frontier[v];
        if (frontier[j * nv + v] > maxf)
          maxf = frontier[j * nv + v];
      }
    }
    i++;
  }

  maxf++;
  for (uint e = 0; e < nu; e++) {
    for (unsigned v = 0; v < nv; v++) {
      if (frontier[e * nv + v] == -1)
        frontier[e * nv + v] = maxf;
    }
  }

  for (uint e = 0; e < nu; e++) {
    for (unsigned v = 0; v < nv; v++)
      pu[e].frontier[v] = frontier[e * nv + v];
  }

  tfree(frontier), array_free(&ids);

  return maxf;
}

void parrsb_fetch_nbrs(unsigned *nei, long long *eids, unsigned nv,
                       long long *vids, double *xyz, double *mask, double *mat,
                       unsigned *nwi, int *frontier, const MPI_Comm comm,
                       unsigned maxne) {
  struct vtx_t {
    ulong id;
    uint p, o, seq;
  };

  size_t ne = *nei, size = nv;
  size *= ne;

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
  uint *proc = tcalloc(uint, ne * 4 + 1), nproc = ne * 4 + 1;
  uint *work = tcalloc(uint, 1024), nwork = 1024;
  pv = (struct vtx_t *)vtx2p.ptr, s = 0;
  uint seq, e;
  while (s < vtx2p.n) {
    seq = pv[s].seq, e = s + 1;
    while (e < vtx2p.n && seq == pv[e].seq)
      e++;

    if (nwork < e - s + 1)
      nwork = e - s + 1, work = trealloc(uint, work, nwork);
    uint n = 1;
    work[0] = pv[s].o;
    for (uint i = s + 1; i < e; i++) {
      if (work[n - 1] != pv[i].o)
        work[n] = pv[i].o, n++;
    }

    offs[seq + 1] = offs[seq] + n;
    if (offs[seq + 1] >= nproc)
      nproc = 3 * offs[seq + 1] / 2 + 1, proc = trealloc(uint, proc, nproc);
    for (uint i = 0, j = offs[seq]; i < n; i++)
      proc[j + i] = work[i];
    s = e;
  }
  array_free(&vtx2p);

  struct array elems;
  array_init(struct elem_t, &elems, size);

  unsigned nw = *nwi;
  if (nw == 0) {
    for (size_t i = 0; i < size; i++)
      frontier[i] = 0;
  }

  unsigned nd = (nv == 8) ? 3 : 2;
  struct elem_t et;
  for (uint e = 0; e < ne; e++) {
    et.eid = eids[e];
    for (unsigned v = 0; v < nv; v++) {
      et.vid[v] = vids[e * nv + v];
      et.mask[v] = mask[e * nv + v];
      et.frontier[v] = frontier[e * nv + v];
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
  tfree(offs), tfree(proc), tfree(work);

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
      if (pe[j].eid != pe[i].eid)
        array_cat(struct elem_t, &uelems, &pe[j], 1), i = j;
      j++;
    }
  }
  array_free(&elems);

  sint err = (uelems.n > maxne), maxr = uelems.n, wrk;
  comm_allreduce(&ci, gs_int, gs_add, &err, 1, &wrk);
  comm_allreduce(&ci, gs_int, gs_max, &maxr, 1, &wrk);
  if (err > 0) {
    if (ci.id == 0) {
      fprintf(stderr, "maxne = %u is too small ! Try maxne larger than %d\n",
              maxne, maxr);
      fflush(stderr);
    }
    buffer_free(&bfr), array_free(&uelems), comm_free(&ci);
    exit(1);
  }

  *nwi = find_frontier(&uelems, ne, eids, nv, &ci, &bfr);

  // Sort elements first by maximum wave number of it vertices, then by global
  // element id.
  struct elem_t *pu = (struct elem_t *)uelems.ptr;
  for (uint e = 0; e < uelems.n; e++) {
    unsigned maxf = 0;
    for (unsigned v = 0; v < nv; v++) {
      if (pu[e].frontier[v] > maxf)
        maxf = pu[e].frontier[v];
    }
    pu[e].seq = maxf;
  }
  sarray_sort_2(struct elem_t, uelems.ptr, uelems.n, seq, 0, eid, 1, &bfr);
  buffer_free(&bfr);

  // Set output arrays.
  *nei = uelems.n;
  pu = (struct elem_t *)uelems.ptr;
  for (uint e = 0; e < uelems.n; e++) {
    eids[e] = pu[e].eid;
    for (unsigned v = 0; v < nv; v++) {
      vids[e * nv + v] = pu[e].vid[v];
      mask[e * nv + v] = pu[e].mask[v];
      frontier[e * nv + v] = pu[e].frontier[v];
      for (unsigned d = 0; d < nd; d++)
        xyz[e * nv * nd + v * nd + d] = pu[e].xyz[v * nd + d];
      for (unsigned u = 0; u < nv; u++)
        mat[e * nv * nv + v * nv + u] = pu[e].mat[v * nv + u];
    }
  }

  struct comm c;
  comm_split(&ci, ci.id, ci.id, &c);
  struct gs_data *gsh = gs_setup(vids, uelems.n * nv, &c, 0, gs_pairwise, 0);
  gs(frontier, gs_int, gs_min, 0, gsh, &bfr);
  gs_free(gsh), comm_free(&c);

  comm_free(&ci), array_free(&uelems);
}

struct eid_t {
  ulong eid;
  uint e;
};

static int binary_search(ulong eid, struct eid_t *pe, uint n) {
  if (n == 0)
    return -1;

  uint l = 0, u = n - 1;
  while (u - l > 1) {
    uint mid = (u + l) / 2;
    if (pe[mid].eid == eid)
      return mid;
    else if (pe[mid].eid < eid)
      l = mid;
    else
      u = mid;
  }

  if (pe[l].eid == eid)
    return l;
  if (pe[u].eid == eid)
    return u;
  return -1;
}

void parrsb_fetch_nbrs_v2(unsigned *nei, long long *eids, unsigned nv,
                          long long *vids, double *xyz, double *mask,
                          double *mat, unsigned nw, int *wids, MPI_Comm comm,
                          unsigned max_ne) {
  size_t ne = *nei;

  // 1. Find neighbor elements of input elements based on vertex connectivity.
  struct vtx_t {
    ulong vid, eid, nid;
    uint p, np;
  };

  struct array vtxs;
  array_init(struct vtx_t, &vtxs, ne * nv);

  struct comm c;
  comm_init(&c, comm);

  struct vtx_t vt = {.np = c.id};
  for (uint e = 0; e < ne; e++) {
    vt.eid = eids[e];
    for (unsigned v = 0; v < nv; v++) {
      vt.vid = vids[e * nv + v], vt.p = vt.vid % c.np;
      array_cat(struct vtx_t, &vtxs, &vt, 1);
    }
  }

  struct crystal cr;
  crystal_init(&cr, &c);
  sarray_transfer(struct vtx_t, &vtxs, p, 1, &cr);

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort(struct vtx_t, vtxs.ptr, vtxs.n, vid, 1, &bfr);

  struct array vtx2e;
  array_init(struct vtx_t, &vtx2e, vtxs.n);

  struct vtx_t *pv = (struct vtx_t *)vtxs.ptr;
  uint s = 0, e;
  while (s < vtxs.n) {
    e = s + 1;
    while (e < vtxs.n && pv[s].vid == pv[e].vid)
      e++;
    for (uint i = s; i < e; i++) {
      vt = pv[i];
      for (uint j = s; j < e; j++) {
        vt.np = pv[j].p, vt.nid = pv[j].eid;
        array_cat(struct vtx_t, &vtx2e, &vt, 1);
      }
    }
    s = e;
  }
  array_free(&vtxs);

  sarray_transfer(struct vtx_t, &vtx2e, p, 0, &cr);
  sarray_sort_2(struct vtx_t, vtx2e.ptr, vtx2e.n, eid, 1, nid, 1, &bfr);

  // 2. Build element to neighbor map and element to processor map for input
  // elements.
  uint *offs = tcalloc(uint, max_ne + 1), max_nbrs = 27 * max_ne;
  ulong *elist = tcalloc(ulong, max_ne), *nbrs = tcalloc(ulong, max_nbrs);
  uint *wlist = tcalloc(uint, max_ne), *proc = tcalloc(uint, max_nbrs);
  uint *plist = tcalloc(uint, max_ne);

  uint cnt = 0;
  pv = (struct vtx_t *)vtx2e.ptr, s = 0, offs[0] = 0;
  while (s < vtx2e.n) {
    elist[cnt] = pv[s].eid, wlist[cnt] = 0, plist[cnt] = c.id;

    e = s + 1;
    while (e < vtx2e.n && pv[s].eid == pv[e].eid)
      e++;

    uint s0 = offs[cnt];
    if (max_nbrs < s0 + e - s) {
      max_nbrs = 3 * s0 / 2 + e - s;
      proc = trealloc(uint, proc, max_nbrs);
      nbrs = trealloc(ulong, nbrs, max_nbrs);
    }

    nbrs[s0] = pv[s].nid, proc[s0] = pv[s].np, s0++;
    for (uint i = s + 1; i < e; i++) {
      if (nbrs[s0 - 1] != pv[i].nid)
        nbrs[s0] = pv[i].nid, proc[s0] = pv[i].np, s0++;
    }
    cnt++, offs[cnt] = s0, s = e;
  }
  // Sanity check.
  assert(cnt == ne);
  array_free(&vtx2e);

  // 3. Put all local elements in frontier array and sort by element id.
  // We will keep updating this and the map as we update the frontier.
  struct array input, frontier;
  array_init(struct elem_t, &input, ne);
  array_init(struct elem_t, &frontier, 3 * ne / 2);

  struct eid_t et;
  for (uint e = 0; e < ne; e++) {
    et.eid = eids[e], et.e = e;
    array_cat(struct eid_t, &frontier, &et, 1);
  }
  sarray_sort(struct eid_t, frontier.ptr, frontier.n, eid, 1, &bfr);
  array_cat(struct eid_t, &input, &frontier, frontier.n);

  // 4. Update the frontier by finding new neighbor elements from the previous
  // frontier.
  struct req_t {
    ulong eid;
    uint p, seq;
  };

  struct array rqsts;
  array_init(struct req_t, &rqsts, ne);

  struct res_t {
    ulong eid, nid;
    uint p, np;
  };

  struct array respns;
  array_init(struct res_t, &respns, rqsts.n * 10);

  uint fs = 0, fe = ne;
  for (uint i = 1; i <= nw; i++) {
    // Find all the new elements appearing in the map in last wave.
    for (uint i = fs; i < fe; i++) {
      for (uint s = offs[i], e = offs[i + 1]; s < e; s++) {
        if (binary_search(nbrs[s], frontier.ptr, frontier.n) == -1) {
          struct eid_t et = {.eid = nbrs[s]};
          array_cat(struct eid_t, &frontier, &et, 1);
          // FIXME: This is bad. Fix it.
          sarray_sort(struct eid_t, frontier.ptr, frontier.n, eid, 1, &bfr);

          struct req_t rt = {.eid = nbrs[s], .p = proc[s]};
          array_cat(struct req_t, &rqsts, &rt, 1);
        }
      }
    }

    // Get the neighbors of the new elements.
    sarray_transfer(struct req_t, &rqsts, p, 1, &cr);
    sarray_sort(struct req_t, rqsts.ptr, rqsts.n, eid, 1, &bfr);
    struct req_t *pr = (struct req_t *)rqsts.ptr;

    for (uint i = 0; i < rqsts.n; i++) {
      int idx = binary_search(pr[i].eid, input.ptr, input.n);
      if (idx < 0 || idx >= ne)
        errx(EXIT_FAILURE, "Couldn't find element: %lld on processor: %d.",
             pr[i].eid, c.id);

      struct res_t rt = {.eid = pr[i].eid, .p = pr[i].p};
      for (uint s = offs[idx], e = offs[idx + 1]; s < e; s++) {
        rt.nid = nbrs[s], rt.np = proc[s];
        array_cat(struct res_t, &respns, &rt, 1);
      }
    }

    sarray_transfer(struct res_t, &respns, p, 1, &cr);
    sarray_sort_2(struct res_t, respns.ptr, respns.n, eid, 1, nid, 1, &bfr);

    // Update the map with the new elements and their neighbors.
    struct res_t *prs = (struct res_t *)respns.ptr;
    fs = fe, s = 0;
    while (s < respns.n) {
      if (fe >= max_ne)
        errx(EXIT_FAILURE, "max_ne: %u is too small.", max_ne);

      elist[fe] = prs[s].eid, plist[fe] = prs[s].p, wlist[fe] = nw, fe++;
      e = s + 1;
      while (e < respns.n && prs[s].eid == prs[e].eid)
        e++;

      offs[fe] = offs[fe - 1] + e - s;
      if (max_nbrs < offs[fe]) {
        max_nbrs = 3 * offs[fe] / 2 + 1;
        proc = trealloc(uint, proc, max_nbrs);
        nbrs = trealloc(ulong, nbrs, max_nbrs);
      }

      for (uint i = 0; i < e - s; i++) {
        proc[offs[fe - 1] + i] = prs[s + i].np;
        nbrs[offs[fe - 1] + i] = prs[s + i].nid;
      }
    }
    rqsts.n = respns.n = 0;
  }
  array_free(&respns), array_free(&frontier), array_free(&input);
  tfree(offs), tfree(proc), tfree(nbrs);

  for (uint i = 0; i < fe; i++) {
    struct req_t rt = {.eid = elist[i], .p = plist[i], .seq = i};
    array_cat(struct req_t, &rqsts, &rt, 1);
  }

  sarray_transfer(struct req_t, &rqsts, p, 1, &cr);
  sarray_sort(struct req_t, rqsts.ptr, rqsts.n, eid, 1, &bfr);
  struct req_t *pr = (struct req_t *)rqsts.ptr;

  struct array elements;
  array_init(struct elem_t, &elements, rqsts.n);

  struct elem_t elmt;
  for (uint i = 0, j = 0; i < rqsts.n; i++) {
    while (j < ne && eids[j] < pr[i].eid)
      j++;
    // Sanity check.
    assert(j < ne && eids[j] == pr[i].eid);

    elmt.eid = eids[j], elmt.p = pr[i].p;
    for (unsigned v = 0; v < nv; v++) {
      elmt.vid[v] = vids[j * nv + v];
      elmt.mask[v] = mask[j * nv + v];
      for (unsigned u = 0; u < nv; u++)
        elmt.mat[v * nv + u] = mat[j * nv * nv + v * nv + u];
    }

    array_cat(struct elem_t, &elements, &elmt, 1);
  }
  array_free(&rqsts);

  sarray_transfer(struct elem_t, &elements, p, 1, &cr);
  sarray_sort(struct elem_t, elements.ptr, elements.n, seq, 0, &bfr);
  struct elem_t *pe = (struct elem_t *)elements.ptr;

  *nei = ne = fe;
  for (uint i = 0; i < ne; i++) {
    eids[i] = elist[i], wids[i] = wlist[i];
    for (unsigned v = 0; v < nv; v++) {
      eids[i * nv + v] = pe[i].vid[v];
      mask[i * nv + v] = pe[i].mask[v];
      for (unsigned u = 0; u < nv; u++)
        mat[i * nv * nv + v * nv + u] = pe[i].mat[v * nv + u];
    }
  }
  array_free(&elements);

  tfree(elist), tfree(wlist), tfree(plist);
  buffer_free(&bfr), crystal_free(&cr), comm_free(&c);

  return;
}

#define fparrsb_fetch_nbrs                                                     \
  FORTRAN_UNPREFIXED(fparrsb_fetch_nbrs, FPARRSB_FETCH_NBRS)
void fparrsb_fetch_nbrs(int *nei, long long *eids, int *nv, long long *vids,
                        double *xyz, double *mask, double *mat, int *nwi,
                        int *frontier, MPI_Fint *comm, int *maxne, int *err) {
  *err = 1;
  MPI_Comm c = MPI_Comm_f2c(*comm);
  unsigned ne = *nei, nw = *nwi;
  parrsb_fetch_nbrs(&ne, eids, *nv, vids, xyz, mask, mat, &nw, frontier, c,
                    *maxne);
  *nei = ne, *nwi = nw, *err = 0;
}

void parrsb_remove_frontier(unsigned *nei, long long *eids, unsigned nv,
                            long long *vids, double *xyz, double *mask,
                            double *mat, unsigned *nwi, int *frontier) {
  struct entry_t {
    ulong r, c;
    double v;
    uint seq, keep;
  };

  unsigned nw = *nwi;
  if (nw == 0)
    return;

  size_t ne = *nei, size = ne;
  size *= nv * nv;

  struct array entries;
  array_init(struct entry_t, &entries, size);

  unsigned *keep = tcalloc(unsigned, ne);

  struct entry_t et;
  for (unsigned e = 0; e < ne; e++) {
    keep[e] = 1;
    for (unsigned v = 0; v < nv; v++) {
      if (frontier[e * nv + v] == nw)
        keep[e] = 0;
    }

    for (unsigned v = 0; v < nv; v++) {
      for (unsigned u = 0; u < nv; u++) {
        if (frontier[e * nv + v] < nw && frontier[e * nv + u] < nw) {
          et.r = vids[e * nv + v], et.c = vids[e * nv + u], et.keep = keep[e];
          et.seq = e * nv * nv + v * nv + u, et.v = mat[et.seq];
          array_cat(struct entry_t, &entries, &et, 1);
        }
      }
    }
  }

  buffer bfr;
  buffer_init(&bfr, 1024);
  sarray_sort_3(struct entry_t, entries.ptr, entries.n, r, 1, c, 1, keep, 0,
                &bfr);

  struct array unique;
  array_init(struct entry_t, &unique, size);

  struct entry_t *pe = (struct entry_t *)entries.ptr;
  uint s = 0;
  while (s < entries.n) {
    uint e = s;
    double v = 0;
    while (e < entries.n && pe[s].r == pe[e].r && pe[s].c == pe[e].c) {
      if (!pe[e].keep) {
        v += pe[e].v;
      } else {
        pe[e].v += v, v = 0;
        array_cat(struct entry_t, &unique, &pe[e], 1);
      }
      e++;
    }
    // Sanity check.
    assert(pe[e - 1].keep == 1);
    s = e;
  }
  array_free(&entries);

  sarray_sort(struct entry_t, unique.ptr, unique.n, seq, 0, &bfr);
  buffer_free(&bfr);

  unsigned nd = (nv == 8) ? 3 : 2;

  uint nk = 0;
  struct entry_t *pu = (struct entry_t *)unique.ptr;
  for (unsigned e = 0; e < ne; e++) {
    if (keep[e]) {
      eids[nk] = eids[e];
      for (unsigned v = 0; v < nv; v++) {
        vids[nk * nv + v] = vids[e * nv + v];
        mask[nk * nv + v] = mask[e * nv + v];
        frontier[nk * nv + v] = frontier[e * nv + v];
        for (unsigned d = 0; d < nd; d++)
          xyz[nk * nv * nd + v * nd + d] = xyz[e * nv * nd + v * nd + d];
        for (unsigned u = 0; u < nv; u++)
          mat[nk * nv * nv + v * nv + u] = pu[nk * nv * nv + v * nv + u].v;
      }
      nk++;
    }
  }
  // Sanity check.
  assert(nk == unique.n / (nv * nv));

  *nwi = nw - 1, *nei = nk;

  array_free(&unique), tfree(keep);
}

#define fparrsb_remove_frontier                                                \
  FORTRAN_UNPREFIXED(fparrsb_remove_frontier, FPARRSB_REMOVE_FRONTIER)
void fparrsb_remove_frontier(int *nei, long long *eids, int *nv,
                             long long *vids, double *xyz, double *mask,
                             double *mat, int *nwi, int *frontier, int *err) {
  *err = 1;
  unsigned ne = *nei, nw = *nwi;
  parrsb_remove_frontier(&ne, eids, *nv, vids, xyz, mask, mat, &nw, frontier);
  *nei = ne, *nwi = nw, *err = 0;
}

#undef tfree
