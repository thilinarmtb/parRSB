#include "con-impl.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

static scalar normi(scalar *a, sint n) {
  scalar norm = 0;
  for (sint i = 0; i < n; i++)
    if (norm < fabs(a[i])) norm = fabs(a[i]);
  return norm;
}

static int test_transform_face_00(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++)
    for (sint j = 0; j < 3; j++)
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {0.0, 0.0, 0.0};

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 3; i++)
    for (sint j = 0; j < 3; j++) R_err[3 * i + j] = R[i][j] - R_expected[i][j];

  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > tol);
  err |= (normi(R_err, 9) > tol);
  return err;
}

static int test_transform_face_01(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;
      face1[i][j] += 10.0;
    }
  }

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {10.0, 10.0, 10.0};

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 3; i++)
    for (sint j = 0; j < 3; j++) R_err[3 * i + j] = R[i][j] - R_expected[i][j];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > tol);
  err |= (normi(R_err, 9) > tol);
  return err;
}

static int test_transform_face_02(const scalar tol) {
  scalar face1[4][3], face0[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = face0[i][j] = rand() / (scalar)RAND_MAX;
      face1[i][j] += j;
    }
  }

  scalar R_expected[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  scalar t_expected[3] = {0.0, 1.0, 2.0};

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 3; i++)
    for (sint j = 0; j < 3; j++) R_err[3 * i + j] = R[i][j] - R_expected[i][j];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > tol);
  err |= (normi(R_err, 9) > tol);
  return err;
}

static int test_transform_face_03(const scalar tol) {
  scalar face0[4][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}};

  // Rotate by 60 degrees around the x-axis and translate by (1, 2, 3).
  scalar theta = M_PI / 3.0;
  scalar R_expected[3][3] = {{1.0, 0.0, 0.0},
                             {0.0, cos(theta), -sin(theta)},
                             {0.0, sin(theta), cos(theta)}};
  scalar t_expected[3] = {1.0, 2.0, 3.0};

  scalar face1[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face1[i][j] = 0;
      for (sint k = 0; k < 3; k++)
        face1[i][j] += R_expected[j][k] * face0[i][k];
      face1[i][j] += t_expected[j];
    }
  }

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar R_err[9], t_err[3];
  for (sint i = 0; i < 3; i++)
    for (sint j = 0; j < 3; j++) R_err[3 * i + j] = R[i][j] - R_expected[i][j];
  for (sint i = 0; i < 3; i++) t_err[i] = t[i] - t_expected[i];

  sint err = 0;
  err |= (normi(t_err, 3) > tol);
  err |= (normi(R_err, 9) > tol);
  return err;
}

static int test_transform_face_04(const scalar tol) {
  scalar face0[4][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 0.0}};
  scalar face1[4][3] = {
      {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, {1.0, 0.0, 1.0}};

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar face2[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face2[i][j] = 0;
      for (sint k = 0; k < 3; k++) face2[i][j] += R[j][k] * face0[i][k];
      face2[i][j] += t[j];
    }
  }

  scalar f_err[12];
  for (sint i = 0; i < 4; i++)
    for (sint j = 0; j < 3; j++) f_err[3 * i + j] = face2[i][j] - face1[i][j];

  return (normi(f_err, 12) > tol);
}

static int test_transform_face_05(const scalar tol) {
  scalar face0[4][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 1.0, 0.0}};
  scalar face1[4][3] = {
      {1.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};

  scalar R[3][3], t[3];
  transform_face(R, t, face1, face0, 4, 3, tol);

  scalar face2[4][3];
  for (sint i = 0; i < 4; i++) {
    for (sint j = 0; j < 3; j++) {
      face2[i][j] = 0;
      for (sint k = 0; k < 3; k++) face2[i][j] += R[j][k] * face0[i][k];
      face2[i][j] += t[j];
    }
  }

  scalar f_err[12];
  for (sint i = 0; i < 4; i++)
    for (sint j = 0; j < 3; j++) f_err[3 * i + j] = face2[i][j] - face1[i][j];

  return (normi(f_err, 12) > tol);
}

int test_transform_face(scalar tol) {
  sint errs = 0;

  errs |= test_transform_face_00(tol);
  errs |= test_transform_face_01(tol);
  errs |= test_transform_face_02(tol);
  errs |= test_transform_face_03(tol);
  errs |= test_transform_face_04(tol);
  errs |= test_transform_face_05(tol);

  return errs;
}

int main(int argc, char *argv[]) { return test_transform_face(1e-14); }
