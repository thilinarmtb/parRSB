/*
Parition mesh using Nek5000's vertex connectivity (con) file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "conReader.h"

#include "gslib.h"
#include "parRSB.h"
#include "quality.h"

#define MAXNV 8 /* maximum number of vertices per element */
typedef struct {
  int proc;
  long long id;
  long long vtx[MAXNV];
} elm_data;

int main(int argc, char *argv[]) {
  struct comm comm;
  struct crystal cr;
  struct array eList;
  elm_data *data;
  struct con con;

  int ierr;
  int e, n, nel, nv;
  int options[3];

  MPI_Init(&argc, &argv);
  comm_init(&comm, MPI_COMM_WORLD);

  ierr = conRead(argv[1], &con, comm.c);
  if(ierr) goto quit;

  nv  = con.nv;
  nel = con.nel;
  int *part = (int*) malloc(nel * sizeof(int));

  options[0] = 1; /* use custom options */
  options[1] = 3; /* debug level        */
  options[2] = 0; /* not used           */

  ierr = parRSB_partMesh(part, con.vl, nel, nv, options, comm.c);
  if(ierr) goto quit;

  /* redistribute data */
  array_init(elm_data, &eList, nel), eList.n = nel;
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    data[e].proc = part[e];
    data[e].id   = con.el[e];
    for(n = 0; n < nv; ++n) {
      data[e].vtx[n] = con.vl[e * nv + n];
    }
  }
  free(part);
  conFree(&con);

#if 0
  crystal_init(&cr, &comm);
  sarray_transfer(elm_data, &eList, proc, 0, &cr);
  crystal_free(&cr);
  nel = eList.n;

  long long *el = (long long*) malloc(nel * sizeof(long long));
  long long *vl = (long long*) malloc(nv * nel * sizeof(long long));
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    el[e] = data[e].id;
    for(n = 0; n < nv; ++n) {
      vl[e * nv + n] = data[e].vtx[n];
    }
  }

  printPartStat(vl, nel, nv, comm.c);

  free(el);
  free(vl);
#endif

  buffer buf; buffer_init(&buf,1024);
  sarray_sort(elm_data,eList.ptr,(unsigned int)nel,id,1,&buf);
  int *eli = (int*) malloc(nel * sizeof(int));
  int *vli = (int*) malloc(nv * nel * sizeof(int));
  for(data = eList.ptr, e = 0; e < nel; ++e) {
    eli[e] = data[e].proc;
    for(n = 0; n < nv; ++n) {
      vli[e * nv + n] = data[e].vtx[n];
    }
  }
  parRSBWriteMapFile(nel,nv,eli,vli,"example.ma2",MPI_COMM_WORLD);

  free(eli);
  free(vli);
  buffer_free(&buf);

  array_free(&eList);
  comm_free(&comm);

quit:
  MPI_Finalize();
  return ierr;
}
