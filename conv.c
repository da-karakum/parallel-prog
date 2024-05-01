// Последовательная программа решения уравнения переноса по схеме прямоугольник

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double calc_t;


calc_t const T = M_PI;
calc_t const X = 2 * M_PI;

calc_t f (calc_t t, calc_t x) {
    return 0; //t * x;
}

calc_t phi (calc_t x) {
    return x < X / 4 ? sin (2 * x) : 0;
}

calc_t psi (calc_t t) {
    return t < M_PI / 4 ? sin (4 * t) : 0;
}

calc_t rt (calc_t lt, calc_t rb, calc_t lb, 
           calc_t aux_sum, calc_t aux_diff, calc_t rhs) {
    calc_t ltrb = lt - rb;
    return lb + (rhs - ltrb * aux_diff) / aux_sum;
}

#define acc(k, m) (net[(k) * sizeX + (m)])


void solve (calc_t *net, int sizeT, int sizeX, calc_t tau, calc_t h, calc_t t0,
            calc_t x0) {

    calc_t aux_sum  = 1 / (2 * tau) + 1 / (2 * h);
    calc_t aux_diff = 1 / (2 * tau) - 1 / (2 * h);

    calc_t t = t0;
    for (int k = 0; k < sizeT; ++k, t += tau)
        acc (k, 0) = psi (t);
    
    calc_t x = x0;
    for (int m = 0; m < sizeX; ++m, x += h)
        acc (0, m) = phi (x);

    t = t0 + tau;
    x = x0 + h;
    for (int k = 1; k < sizeT; ++k, t += tau) 
        for (int m = 1; m < sizeX; ++m, x += h)
            acc (k, m) = rt (acc(k, m-1), acc(k-1, m), acc(k-1, m-1), 
                             aux_sum, aux_diff, f (t - tau / 2, x + h / 2));

    return;
}

int main (int argc, char *argv[]) {
    int sizeT = atoi (argv[1]);
    int sizeX = atoi (argv[2]);

    calc_t tau = T / (sizeT - 1);
    calc_t   h = X / (sizeX - 1);

    calc_t *net = (calc_t *) malloc (sizeT * sizeX * sizeof (calc_t));

    if (net == NULL) {
        perror ("malloc");
        return 1;
    }

    solve (net, sizeT, sizeX, tau, h, 0, 0);

    // t axis points
    calc_t t = 0;
    for (int k = 0; k < sizeT; ++k, t += tau)
        printf ("%.3f ", t);
    printf ("\n");

    // x axis points
    calc_t x = 0;
    for (int m = 0; m < sizeX; ++m, x += h)
        printf ("%.3f ", x);
    printf ("\n");

    // net solution
    for (int k = 0; k < sizeT; ++k) {
        for (int m = 0; m < sizeX; ++m)
            printf ("%.3f ", acc (k, m));
        printf ("\n");
    }


    free (net);
}