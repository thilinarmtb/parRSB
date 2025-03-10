#include "con-impl.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  test_transform_face(1e-14, comm);

  MPI_Finalize();
}
