#ifndef DIM
#error "DIM is undefined"
#endif

#define transform_points GS_TOKEN_PASTE(transform_points_, DIM)
static void transform_points(scalar *coordo, sint nf, const sint *const bid,
                             sint nv, const scalar *const coordi,
                             const scalar R[4][3]) {
  const size_t nc = (size_t)nf * nv;
  const scalar t[3] = {R[3][0], R[3][1], R[3][2]};

  for (uint i = 0; i < nc; i++) {
    const scalar *coordi_i = &coordi[i * DIM];
    scalar *coordo_i = &coordo[i * DIM];

    if (bid[i] == 0) {
      for (int j = 0; j < DIM; j++) coordo_i[j] = coordi_i[j];
    } else if (bid[i] == 1) {
      for (int j = 0; j < DIM; j++)
        for (int k = 0; k < DIM; k++) coordo_i[k] += R[j][k] * coordi_i[k];
    }
  }
}
