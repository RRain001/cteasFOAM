// M0 MPI toolchain test only. No OpenFOAM solver or numerical validation.
#include <mpi.h>
#include <iostream>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = -1;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int value = rank + 1;
    int sum = 0;
    MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    const bool ok = size == 2 && sum == 3;
    if (rank == 0)
    {
        std::cout << (ok ? "PASS" : "FAIL")
                  << ": MPI ranks=" << size << " allreduce=" << sum << '\n';
    }
    MPI_Finalize();
    return ok ? 0 : 1;
}
