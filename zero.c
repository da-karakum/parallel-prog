// Написать программу, которая измеряет задержку передачи сообщения между
// двумя узлами сети с помощью технологии MPI
// Вывод:
// rank 0: measured latency 40749 nsec
// rank 1: received Cat

#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#define BUFFER_LEN 4
#define MESSAGE_1 "Cat"
#define MESSAGE_2 "Dog"

#define arrsize(a) (sizeof(a) / sizeof(*a))
static_assert (arrsize(MESSAGE_1) <= BUFFER_LEN);
static_assert (arrsize(MESSAGE_2) <= BUFFER_LEN);

long diff_nsec (const struct timespec *, const struct timespec *);


int main (int argc, char *argv[]) {
    MPI_Init (&argc, &argv);
    int size, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    struct timespec sent, received;

    if (rank == 0) {
        clock_gettime (CLOCK_MONOTONIC, &sent);
        MPI_Send (MESSAGE_1, BUFFER_LEN, MPI_CHAR, 1, 42, MPI_COMM_WORLD);
        MPI_Recv (&received, sizeof (received), MPI_BYTE, 1, 42, 
            MPI_COMM_WORLD, NULL);

        long elapsed = diff_nsec (&sent, &received);

        printf ("rank 0: measured latency %ld nsec\n", elapsed);
    } else {
        char buf [BUFFER_LEN];
        MPI_Recv (buf, BUFFER_LEN, MPI_CHAR, 0, 42, MPI_COMM_WORLD, NULL);
        clock_gettime (CLOCK_MONOTONIC, &received);
        MPI_Send (&received, sizeof (received), MPI_BYTE, 0, 42, 
            MPI_COMM_WORLD);

        printf ("rank 1: received %s\n", buf);
    }
    
    MPI_Finalize ();
}

long diff_nsec (const struct timespec *from, const struct timespec *to) {
    long diff_sec = to->tv_sec - from->tv_sec;
    long diff_nsec = to->tv_nsec - from->tv_nsec;
    return diff_sec * (int)(1e9 + 0.5) + diff_nsec;
}
