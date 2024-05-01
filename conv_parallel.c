// Параллельная программа решения уравнения переноса по схеме прямоугольник

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#ifndef NDEBUG
#define log(...) printf (__VA_ARGS__)
#else
#define log(...)
#endif

typedef double calc_t;


calc_t const T = M_PI;
calc_t const X = 2 * M_PI;

calc_t f (calc_t t, calc_t x) {
#   ifndef NO_DELAY
    for (int i = 0; i < 10000; ++i)
        ;
#   endif
    return 0; //t * x;
}

calc_t phi (calc_t x) {
    return x < M_PI / 2 ? sin (2 * x) : 0;
}

calc_t psi (calc_t t) {
    return t < M_PI / 4 ? sin (4 * t) : 0;
}

calc_t rt (calc_t lt, calc_t rb, calc_t lb, 
           calc_t aux_sum, calc_t aux_diff, calc_t rhs) {
    calc_t ltrb = lt - rb;
    return lb + (rhs - ltrb * aux_diff) / aux_sum;
}

#define acc(k, m) (net[(k) * sizeX + (m) - left])


void solve (calc_t *net, int sizeT, int left, int right, calc_t tau, calc_t h, 
            int rank, int comm_size) {
    
    int const sizeX = right - left;
    calc_t const aux_sum  = 1 / (2 * tau) + 1 / (2 * h);
    calc_t const aux_diff = 1 / (2 * tau) - 1 / (2 * h);
    int k, m;
    calc_t t, x;
    
    // Calculating t=0 boundary
    for (m = left, x = left * h; m < right; ++m, x += h)
        acc (0, m) = phi (x);

    // Array for rightmost vertical of left peer
    calc_t *from_peer = NULL;
    if (rank != 0) {
        from_peer = (calc_t *) malloc (sizeT * sizeof (*from_peer));
        from_peer[0] = phi ((left - 1) * h);
    }

    
    for (k = 1, t = tau; k < sizeT; ++k, t += tau) {
        m = left, x = 0;
        if (rank == 0)
            acc (k, m) = psi (t);
        else {
            MPI_Recv (from_peer + k, sizeof (*from_peer), MPI_BYTE, rank - 1,
                        42, MPI_COMM_WORLD, NULL);
            acc (k, m) = rt (from_peer[k], acc(k-1, m), from_peer[k-1], 
                             aux_sum, aux_diff, f (t - tau / 2, x + h / 2));
        }


        for (m = left + 1, x = h; m < right; ++m, x += h)
            acc (k, m) = rt (acc(k, m-1), acc(k-1, m), acc(k-1, m-1), 
                             aux_sum, aux_diff, f (t - tau / 2, x + h / 2));

        if (rank < comm_size - 1)
            MPI_Bsend (&acc (k, right - 1), sizeof(calc_t), MPI_BYTE, rank + 1,
                        42, MPI_COMM_WORLD);
    }

    free (from_peer);

    return;
}

void calcBounds (int rank, int comm_size, int sizeX, int *left, int *right);
int maxBound (int size_global, int comm_size);

int main (int argc, char *argv[]) {
    MPI_Init (&argc, &argv);
    int comm_size, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    int sizeT = atoi (argv[1]);
    int sizeX_global = atoi (argv[2]);

    if (sizeX_global < comm_size) {
        if (rank == 0) 
            printf ("Error: sizeX must not be less than np\n");
        return 1;
    }

    calc_t tau = T / (sizeT - 1);
    calc_t   h = X / (sizeX_global - 1);

    int left, right, sizeX;
    calcBounds (rank, comm_size, sizeX_global, &left, &right);
    sizeX = right - left;

    calc_t *net;
    if (rank != 0)
        net = (calc_t *) malloc (sizeT * (right - left) * sizeof (calc_t));
    else
        net = (calc_t *) malloc (sizeT * maxBound (sizeX_global, comm_size) * sizeof (calc_t));

    if (net == NULL) {
        perror ("malloc");
        return 1;
    }

    solve (net, sizeT, left, right, tau, h, rank, comm_size);

#   ifndef NO_RESULTS

    if (rank != 0) {
        MPI_Send (net, sizeT * sizeX * sizeof(calc_t), MPI_BYTE, 0, 42, 
            MPI_COMM_WORLD);
        goto cleanup;
    }

    calc_t *net_global = 
        (calc_t *) malloc (sizeT * sizeX_global * sizeof(calc_t));

    int k, m;

    for (k = 0; k < sizeT; ++k)
        for (m = 0; m < right; ++m)
            net_global[k * sizeX_global + m] = acc (k, m);

    for (int peer_rank = 1; peer_rank < comm_size; ++peer_rank) {
        int left, right;
        calcBounds (peer_rank, comm_size, sizeX_global, &left, &right);
        size_t sizeX_peer = right - left; 
        MPI_Recv (net, sizeT * sizeX_peer * sizeof(calc_t), MPI_BYTE, 
            peer_rank, 42, MPI_COMM_WORLD, NULL);
        for (k = 0; k < sizeT; ++k)
            for (m = left; m < right; ++m)
                net_global[k * sizeX_global + m] = net[k * sizeX_peer + m];
    }

    // t axis points
    calc_t t = 0;
    for (int k = 0; k < sizeT; ++k, t += tau)
        printf ("%.3f ", t);
    printf ("\n");

    // x axis points
    calc_t x = left * h;
    for (int m = 0; m < sizeX_global; ++m, x += h)
        printf ("%.3f ", x);
    printf ("\n");

    // net solution
    for (int k = 0; k < sizeT; ++k) {
        for (int m = 0; m < sizeX_global; ++m)
            printf ("%.3f ", net_global[k * sizeX_global + m]);
        printf ("\n");
    }

    free (net_global);

#   endif // !NO_RESULTS

    cleanup:
    free (net);
    MPI_Finalize ();
}

void calcBounds (int rank, int comm_size, int sizeX, int *left, int *right) {
    double leftd  = (double)  rank      / comm_size * sizeX;
    double rightd = (double) (rank + 1) / comm_size * sizeX;
    *left  = -(int)(-leftd - 0.5); // round up
    *right = (int)(rightd + 0.5); // round down
}

int maxBound (int size_global, int comm_size) {
    return size_global / comm_size + 1;
}
