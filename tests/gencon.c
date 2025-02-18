#include "con-impl.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  test_automatic_periodic_face_match(NULL, 0, NULL, 0, 0, NULL, 1e-15, comm);

  MPI_Finalize();
}
