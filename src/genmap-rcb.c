#include <sort.h>
#include <float.h>
#include <parRSB.h>
#include <genmap-impl.h>

void get_axis_len_local(double *min, double *max, void *elems, uint nel, int ndim)
{
  sint i;
  for(i = 0; i < ndim; i++)
    min[i] = DBL_MAX, max[i] = -DBL_MAX;

  size_t unit_size;
  unsigned char *type=elems;
  if (*type == GENMAP_RCB_ELEMENT) {
    unit_size = sizeof(struct rcb_element);
  } else if(*type == GENMAP_RSB_ELEMENT) {
    unit_size = sizeof(struct rsb_element);
  }

  struct rcb_element *elem;
  for(i = 0; i < nel; i++) {
    elem = (struct rcb_element *) ((char *)elems + i*unit_size);
    if (elem->coord[0] < min[0])
      min[0] = elem->coord[0];
    if (elem->coord[0] > max[0])
      max[0] = elem->coord[0];

    if (elem->coord[1] < min[1])
      min[1] = elem->coord[1];
    if (elem->coord[1] > max[1])
      max[1] = elem->coord[1];
  }

  if (ndim == 3) {
    for(i = 0; i < nel; i++) {
      elem = (struct rcb_element *) ((char *)elems + i*unit_size);
      if (elem->coord[2] < min[2])
        min[2] = elem->coord[2];
      if (elem->coord[2] > max[2])
        max[2] = elem->coord[2];
    }
  }
}

void get_axis_len(double *length, struct array *a, struct comm *c, int ndim)
{
  double min[MAXDIM], max[MAXDIM], buf[MAXDIM];

  get_axis_len_local(min, max, a->ptr, a->n, ndim);
  comm_allreduce(c, gs_double, gs_min, min, MAXDIM, buf);
  comm_allreduce(c, gs_double, gs_max, max, MAXDIM, buf);

  sint i;
  for(i = 0; i < ndim; i++)
    length[i] = max[i] - min[i];
}

int rcb_level(struct comm *c, struct array *a, int ndim){
  double length[MAXDIM];

  if (c->np == 1)
    return 0;

  metric_tic(c, AXISLEN);
  get_axis_len(length, a, c, ndim);
  metric_toc(c, AXISLEN);

  int axis1 = 0, d;
  for (d = 1; d < ndim; d++)
    if (length[d] > length[axis1])
      axis1 = d;

  metric_tic(c, PARSORT);
  unsigned char *type = a->ptr;
  if (*type == GENMAP_RCB_ELEMENT) {
    switch(axis1){
      case 0:
        parallel_sort(struct rcb_element, a, coord[0], gs_double, 0, 1, c);
        break;
      case 1:
        parallel_sort(struct rcb_element, a, coord[1], gs_double, 0, 1, c);
        break;
      case 2:
        parallel_sort(struct rcb_element, a, coord[2], gs_double, 0, 1, c);
        break;
      default:
        break;
    }
  } else if (*type == GENMAP_RSB_ELEMENT) {
    switch(axis1){
      case 0:
        parallel_sort(struct rsb_element, a, coord[0], gs_double, 0, 1, c);
        break;
      case 1:
        parallel_sort(struct rsb_element, a, coord[1], gs_double, 0, 1, c);
        break;
      case 2:
        parallel_sort(struct rsb_element, a, coord[2], gs_double, 0, 1, c);
        break;
      default:
        break;
    }
  }
  metric_toc(c, PARSORT);

  return 0;
}

int genmap_rcb(genmap_handle h) {
  GenmapComm local_c = GenmapGetLocalComm(h);
  struct comm *lc = &local_c->gsc;

  int ndim = (h->nv == 8) ? 3 : 2;

  rcb(lc, h->elements, ndim);

  struct rcb_element *eptr = h->elements->ptr;
  int e;
  for(e = 0; e < h->elements->n; e++)
    eptr[e].seq = e;

  return 0;
}

int rcb(struct comm *ci, struct array *elements, int ndim) {
  struct comm c;
  comm_dup(&c, ci);

  int size = c.np;
  int rank = c.id;

  while (size > 1) {
    rcb_level(&c, elements, ndim);

    int p = (size + 1) / 2;
    int bin = (rank >= p);

    comm_ext comm_rcb;
#ifdef MPI
    MPI_Comm_split(c.c, bin, rank, &comm_rcb);
#else
    comm_rcb = 1;
#endif

    comm_free(&c);
    comm_init(&c, comm_rcb);

#ifdef MPI
    MPI_Comm_free(&comm_rcb);
#endif

    size = c.np;
    rank = c.id;
  }

  buffer bfr;
  buffer_init(&bfr, 1024);

  rcb_local(elements, 0, elements->n, ndim, &bfr);

  buffer_free(&bfr);

  comm_free(&c);

  return 0;
}

void rcb_local(struct array *a, uint start, uint end, int ndim, buffer *buf) {
  sint size = end-start;
  assert(size >= 0);

  if (size <= 2)
    return;

  size_t unit_size;
  unsigned char *type = a->ptr;
  if (*type == GENMAP_RCB_ELEMENT) {
    unit_size = sizeof(struct rcb_element);
  } else if(*type==GENMAP_RSB_ELEMENT) {
    unit_size = sizeof(struct rsb_element);
  }

  void *st = (void *)a->ptr + unit_size*start;
  double length[3], min[3], max[3];
  get_axis_len_local(min, max, st, size, ndim);

  sint i;
  for(i = 0; i < ndim; i++)
    length[i] = max[i] - min[i];

  int axis = 0;
  if (fabs(length[axis]) < fabs(length[1]))
    axis = 1;
  if (ndim == 3)
    if (fabs(length[axis]) < fabs(length[2]))
      axis = 2;

  if (*type == GENMAP_RCB_ELEMENT) {
    switch (axis) {
      case 0:
        sarray_sort(struct rcb_element, st, size, coord[0], 3, buf);
        break;
      case 1:
        sarray_sort(struct rcb_element, st, size, coord[1], 3, buf);
        break;
      case 2:
        sarray_sort(struct rcb_element, st, size, coord[2], 3, buf);
        break;
      default:
        break;
    }
  } else if (*type==GENMAP_RSB_ELEMENT) {
    switch (axis) {
      case 0:
        sarray_sort(struct rsb_element, st, size, coord[0], 3, buf);
        break;
      case 1:
        sarray_sort(struct rsb_element, st, size, coord[1], 3, buf);
        break;
      case 2:
        sarray_sort(struct rsb_element, st, size, coord[2], 3, buf);
        break;
      default:
        break;
    }
  }

  uint mid = (start+end)/2;
  rcb_local(a, start, mid, ndim, buf);
  rcb_local(a, mid  , end, ndim, buf);
}
