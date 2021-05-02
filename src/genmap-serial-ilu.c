#include <genmap-impl.h>
#include <genmap-iterative.h>
#include <genmap-precond.h>

int genmap_serial_ilu(genmap_handle h) {
  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_number_faces_and_edges(h, gc);
  genmap_comm_scan(h, gc);

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int np = gc->np;
  assert(np == 1);

  metric_tic(lc, RCB);
  rcb(lc, h->elements, ndim, &h->buf);
  metric_toc(lc, RCB);

  uint lelt = genmap_get_nel(h);
  genmap_vector init;
  genmap_vector_create(&init, lelt);

  struct rsb_element *elements = genmap_get_elements(h);
  uint i;
  for (i = 0; i < lelt; i++)
    init->data[i] = genmap_get_local_start_index(h) + i + 1;

  genmap_vector_ortho_one(gc, init, lelt);

  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  genmap_laplacian_init(h, gc);

  metric_tic(gc, PRECONDSETUP);
  struct precond *d = precond_setup(2, h, gc);
  metric_toc(gc, PRECONDSETUP);

  metric_tic(gc, PROJECT);
  int ppfi = project(h, gc, d, init, 100, y);
  metric_toc(gc, PROJECT);
  metric_acc(NPROJECT, ppfi);

  precond_free(d);

  genmap_destroy_vector(y);
  genmap_destroy_vector(init);

  return 0;
}
