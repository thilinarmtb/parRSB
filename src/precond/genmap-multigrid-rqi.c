#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-multigrid-precon.h>

// Input z should be orthogonal to 1-vector, have unit norm.
// RQI should not change z.
int rqi(genmap_handle h, genmap_comm c, mgData d, GenmapVector z, int max_iter, GenmapVector y)
{
  assert(z->size==y->size);

  uint lelt = z->size;
  GenmapVector err;
  GenmapCreateVector(&err,lelt);

  struct comm *gsc = &c->gsc;
  int rank = genmap_comm_rank(GenmapGetGlobalComm(h));

  GenmapLong nelg = genmap_get_global_nel(h);

  // Grammian
  GenmapScalar *Z, *GZ, *M, *rhs, *v, *buf;
  GenmapMalloc(max_iter*lelt, &Z);
  GenmapMalloc(lelt, &GZ);
  GenmapMalloc(max_iter*max_iter, &M);
  GenmapMalloc(max_iter, &rhs);
  GenmapMalloc(max_iter, &v);
  GenmapMalloc(max_iter*max_iter, &buf);

  metric_tic(gsc, PROJECT);
  int ppfi = project(h, c, d, z, 100, y);
  metric_toc(gsc, PROJECT);
  metric_acc(NPROJECT, ppfi);

  uint i, j, k, l;
  for (i = 0; i < 50; i++) {
    GenmapScalar norm = GenmapDotVector(y,y);
    comm_allreduce(gsc, gs_double, gs_add, &norm, 1, buf);
    GenmapScalar normi = 1.0/sqrt(norm);

    GenmapAxpbyVector(z, z, 0.0, y, normi);
    GenmapOrthogonalizebyOneVector(c, z, nelg);

    int N = i + 1;

    if (h->options->rsb_grammian == 1) {
      metric_tic(gsc,GRAMMIAN);
      //if k>1;
      //  Z(:,k)=z-Z(:,1:k-1)*(Z(:,1:k-1)'*z);
      //  Z(:,k)=Z(:,k)/norm(Z(:,k));
      //end;
      if (i>0) {
        // rhs = Z[1:k-1,:]*z
        for (j=0; j<i; j++) {
          rhs[j]=0.0;
          for(l=0; l<lelt; l++)
            rhs[j]+=Z[j*lelt+l]*z->data[l];
        }
        // Global reduction rhs[j]
        comm_allreduce(gsc,gs_double,gs_add,rhs,i,buf);

        //Z[k,:] = z[:] - Z[:,1:lelt]*rhs[:]
        for (l=0; l<lelt; l++)
          Z[i*lelt+l] = z->data[l];

        for (j=0; j<i; j++) {
          for (l=0; l<lelt; l++)
            Z[i*lelt+l] = Z[i*lelt+l] - rhs[j]*Z[j*lelt+l];
        }

        //Z[k,:]= Z[k,:]/||Z[k,:]||
        norm=0.0;
        for (l=0; l<lelt; l++)
          norm+=Z[i*lelt+l]*Z[i*lelt+l];

        comm_allreduce(gsc,gs_double,gs_add,&norm,1,buf);
        norm = 1.0/sqrt(norm);

        for (l=0; l<lelt; l++)
          Z[i*lelt+l]*=norm;

        //M=Z(1:k,:)*G*Z(1:k,:);
        for (j=0; j<N; j++) {
          GenmapLaplacian(h,c,&Z[j*lelt],GZ);
          for (k=0; k<N; k++) {
            M[k*N+j]=0.0;
            for (l=0; l<lelt; l++)
              M[k*N+j] += Z[k*lelt+l]*GZ[l];
          }
        }

        // Global reduction of M
        comm_allreduce(gsc,gs_double,gs_add,M,N*N,buf);

        // Inverse power iterarion on M
        genmap_inverse_power(v, N, M, 0);

        for (j=0; j<lelt; j++)
          z->data[j] = 0.0;

        for (j=0; j<N; j++) {
          for(k=0; k<lelt; k++)
            z->data[k]+=Z[j*lelt+k]*v[j];
        }
        GenmapOrthogonalizebyOneVector(c,z,nelg);
      } else {
        //Z(k,:) = z;
        for(l=0; l<lelt; l++)
          Z[i*lelt+l]=z->data[l];
      }
      metric_toc(gsc,GRAMMIAN);
    }

    metric_tic(gsc, PROJECT);
    ppfi = project(h, c, d, z, 100, y);
    metric_toc(gsc, PROJECT);
    metric_acc(NPROJECT, ppfi);

    GenmapOrthogonalizebyOneVector(c, y, nelg);

    GenmapScalar lambda = GenmapDotVector(y,z);
    GenmapGop(c,&lambda,1,GENMAP_SCALAR,GENMAP_SUM);

    GenmapAxpbyVector(err, y, 1.0, z, -lambda);
    GenmapScalar norme = GenmapDotVector(err, err);
    GenmapGop(c, &norme, 1, GENMAP_SCALAR, GENMAP_SUM);
    norme = sqrt(norme);

    metric_acc(END + i, norme);

    if (ppfi == 1)
      break;
  }

  GenmapFree(Z);
  GenmapFree(GZ);
  GenmapFree(M);
  GenmapFree(rhs);
  GenmapFree(v);
  GenmapFree(buf);

  GenmapDestroyVector(err);

  return i + 1;
}
