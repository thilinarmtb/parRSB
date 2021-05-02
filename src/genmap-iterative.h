#ifndef _GENMAP_ITERATIVE_H_
#define _GENMAP_ITERATIVE_H_

#include <math.h>
#include <stdio.h>

#include <genmap-impl.h>
#include <genmap-partition.h>
#include <genmap-precond.h>

int project(genmap_handle h, struct comm *c, struct precond *d, genmap_vector r,
            int max_iter, genmap_vector x);

int flex_cg(genmap_handle h, struct comm *c, struct precond *d, genmap_vector r,
            int max_iter, genmap_vector x);

int rqi(genmap_handle h, struct comm *c, genmap_vector z, int max_iter,
        genmap_vector fiedler);

#endif
