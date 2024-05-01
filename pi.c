// Параллельное вычисление числа пи

#include <stdio.h>
#include <mpi.h>
#include <stdint.h> // int64_t
#include <stdlib.h> // atoll

void calcBounds (int rank, int size, int64_t iters, int64_t *from, int64_t *to);

int main (int argc, char *argv[]) {
    MPI_Init (&argc, &argv);

    int rank, size;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    int64_t iters = (int64_t)(8e7), from, to;
    calcBounds (rank, size, iters, &from, &to);

    printf ("rank %d: from = %10ld, to = %10ld\n", rank, from, to);

    double ans = 1.0;
    for (unsigned long i = from; i < to; ++i) {
        ans = ans * 2*i * 2*i / (2*i - 1) / (2*i + 1);
    }

    printf ("rank %d: calculated %.19lf\n", rank, ans);

    if (rank > 0)
        MPI_Send (&ans, 1, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
    else {
        for (int peer_rank = 1; peer_rank < size; ++peer_rank) {
            double peer_ans;
            MPI_Recv (&peer_ans, 1, MPI_DOUBLE, peer_rank, 42, MPI_COMM_WORLD, NULL);
            ans *= peer_ans;
        }
        printf ("Final answer: %.19lf\n", 2 * ans);
    }

    MPI_Finalize ();


}

void calcBounds (int rank, int size, int64_t iters, int64_t *from, int64_t *to) {
    double leftd  = (double)  rank      / size * (iters - 1);
    double rightd = (double) (rank + 1) / size * (iters - 1);
    *from  = -(int64_t)(-leftd - 0.5) + 1; // round up
    *to = (int64_t)(rightd + 0.5) + 1; // round down
}