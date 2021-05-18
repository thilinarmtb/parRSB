#ifndef _GENMAP_ELEMENTS_H_
#define _GENMAP_ELEMENTS_H_

#include <genmap-types.h>

/*
 Preprocessor Corner notation:      Symmetric Corner notation:

         4+-----+3    ^ s                    3+-----+4    ^ s
         /     /|     |                      /     /|     |
        /     / |     |                     /     / |     |
      8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
       |     | /     /                     |     | /     /
       |     |/     /                      |     |/     /
      5+-----+6    t                      5+-----+6    t



                   i) Preprocessor notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    3   / |     |
                             4--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     1

                  ii) Symmetric notation:

                                     +--------+     ^ S
                                    /        /|     |
                                   /    4   / |     |
                             1--> /        /  |     |
                                 +--------+ 2 +     +----> R
                                 |        |  /     /
                                 |    6   | /     /
                                 |        |/     /
                                 +--------+     T
                                     3

*/

/* Upper bounds for faces, vertices, and edges */
#define GC_MAX_VERTICES 8
#define GC_MAX_EDGES 12
#define GC_MAX_FACES 6

/* Upper bounds for faces */
#define GC_MAX_FACE_VERTICES 4
#define GC_MAX_EDGE_VERTICES 2

/* Upper bound for number of dimensions */
#define GC_MAX_DIM 3
#define GC_MAX_NEIGHBORS 3

/* 1 - indexed */
extern int faces3D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];
extern int faces2D[GC_MAX_FACES][GC_MAX_FACE_VERTICES];

/* 0 - indexed */
extern int edges3D[GC_MAX_EDGES][GC_MAX_EDGE_VERTICES];

/* rcb_element is used for rcb and rib */
struct rcb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[GC_MAX_DIM];
};

/* rsb_element should be a superset of rcb_element */
struct rsb_element {
  int type;
  GenmapInt proc;
  GenmapInt origin;
  GenmapInt seq;
  GenmapLong globalId;
  GenmapScalar coord[GC_MAX_DIM];
  GenmapScalar fiedler;
  GenmapLong vertices[GC_MAX_VERTICES + GC_MAX_EDGES + GC_MAX_FACES];
  GenmapInt part;
  int level;
};

#endif
