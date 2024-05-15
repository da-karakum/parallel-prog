/*  Параллельная программа решения уравнения переноса 
    u_t + u_x = f(t,x), 0 < t < T, 0 < x < X
    U(0,x) = phi(x)
    u(t,0) = psi(t)
    разными схемами. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <pthread.h>
#include <fcntl.h>
// -lm -lpthread

typedef double calc_t;


calc_t const T = M_PI;
calc_t const X = 2 * M_PI;

calc_t f (calc_t t, calc_t x) {
    return 0; //t * x;
}

calc_t phi (calc_t x) {
    return x > 0 && x < M_PI / 2 ? sin (2 * x) : 0;
}

calc_t psi (calc_t t) {
    return 0;
}

calc_t solution (calc_t t, calc_t x) {
    return phi (x - t);
}

// Схема прямоугольник
// lt------rt
// |        |
// |        |
// |        |
// lb------rb
calc_t rt (calc_t lt, calc_t rb, calc_t lb, 
           calc_t aux_sum, calc_t aux_diff, calc_t rhs) {
    calc_t ltrb = lt - rb;
    return lb + (rhs - ltrb * aux_diff) / aux_sum;
}

// Схема явный уголок
//         corner
//         |
//         |
//         |
// l-------r
calc_t corner (calc_t l, calc_t r, calc_t f, calc_t tau, calc_t tauoverh) {
    return tau * f + (tauoverh) * l + (1 - tauoverh) * r;
}

// Схема Лакса-Вендроффа с аппроксимационной вязкостью
// lt-----lv
//        |
//        |
//        |
// l------m------r
// устойчива при tau <= h
calc_t lv (calc_t l, calc_t m, calc_t r, calc_t f, 
        calc_t kl, calc_t km, calc_t kr, calc_t kf) {
    return kf*f + kl*l + km*m + kr*r;
}

#define acc(k, m) (net[(k) * sizeX + (m)])

int semid;
calc_t *net;
int sizeT;
int sizeX;
calc_t tau;
calc_t h;
int nproc;

typedef struct {
    int id;
    int left;
    int right;
} ThreadData;

static int sem_lock (int num) {
    if (__builtin_expect(nproc == 0, 0))
        return 0;
    struct sembuf oper = { 
        .sem_flg = 0,
        .sem_num = num,
        .sem_op = -1,
    };
    return semop (semid, &oper, 1);
}

static int sem_unlock (int num) {
    if (__builtin_expect(nproc == 0, 0))
        return 0;
    struct sembuf oper = { 
        .sem_flg = 0,
        .sem_num = num,
        .sem_op = +1,
    };
    return semop (semid, &oper, 1);
}

void set_affinity (int core) {
    cpu_set_t cpuSet;
    pthread_t tid = pthread_self ();

    CPU_ZERO (&cpuSet);
    CPU_SET (core, &cpuSet);

    if (pthread_setaffinity_np (tid, sizeof (cpu_set_t), &cpuSet) != 0)
    {
        fprintf (stderr, "pthread_setaffinity_np, thread %d: ", core);
        perror(NULL);
    }
}

void solve (int id, int left, int right);

void *slave (void *context) {
    ThreadData *data = (ThreadData *)context;
    set_affinity (data->id);
    solve (data->id, data->left, data->right);
    pthread_exit(NULL);
}

void solve (int id, int left, int right) {

    // Рассчитываем коэффициенты схемы
#   if defined (REC)
    calc_t const aux_sum  = 1 / (2 * tau) + 1 / (2 * h);
    calc_t const aux_diff = 1 / (2 * tau) - 1 / (2 * h);
#   elif defined (LCOR)
    calc_t tauoverh = tau / h;
#   elif defined (LW)
    calc_t kf = tau;
    calc_t km = 1 - tau*tau/h/h;
    calc_t kl = tau*tau/2/h/h + tau/2/h;
    calc_t kr = tau*tau/2/h/h - tau/2/h;
#   elif defined (EXACT)
#   else
    #error Must define the scheme
#   endif
    int k, m;
    calc_t t, x;
    
    // Calculating t=0 boundary
    for (m = left, x = left * h; m < right; ++m, x += h)
        acc (0, m) = phi (x);
    
    for (k = 1, t = tau; k < sizeT; ++k, t += tau) {
        m = left, x = left * h;
        if (id == 0) 
            acc (k, m) = psi (t);
        else {
            sem_lock (id - 1);
            acc (k, m) = 
#           if defined (REC)
            rt (acc(k, m-1), acc(k-1, m), acc(k-1, m-1), 
                             aux_sum, aux_diff, f (t - tau / 2, x + h / 2));
#           elif defined (LCOR)
            corner (acc(k-1, m-1), acc(k-1, m), f(t, x) , tau, tauoverh);
#           elif defined (LW)
            lv (acc(k-1, m-1), acc(k-1, m), acc(k-1, m+1), f(t,x), 
                             kl, km, kr, kf);
#           elif defined (EXACT)
            solution (t, x);                 
#           else
            NAN;
#           endif
        }


        for (m = left + 1, x = (left + 1) * h; m < right; ++m, x += h) {
            acc (k, m) = 
#           if defined (REC)
            rt (acc(k, m-1), acc(k-1, m), acc(k-1, m-1),
                             aux_sum, aux_diff, f (t - tau / 2, x + h / 2));
#           elif defined (LCOR)
            corner (acc(k-1, m-1), acc(k-1, m), f(t, x) , tau, tauoverh);
#           elif defined (LW)
            lv (acc(k-1, m-1), acc(k-1, m), acc(k-1, m+1), f(t,x), 
                             kl, km, kr, kf);
#           elif defined (EXACT)
            solution (t, x);  
#           else
            NAN;
#           endif
        }
        if (id < nproc - 1) 
            sem_unlock (id);
    }

    return;
}

void calcBounds (int rank, int comm_size, int *left, int *right);
int maxBound (int size_global, int comm_size);

int main (int argc, char *argv[]) {

    sizeT = atoi (argv[1]);
    sizeX = atoi (argv[2]);
    tau = T / (sizeT - 1);
    h = X / (sizeX - 1);

    nproc = atoi (argv[3]);

    if (sizeX < nproc) {
        fprintf (stderr, "Error: sizeX must not be less than nproc\n");
        return -1;
    }

    int nCPUs = sysconf (_SC_NPROCESSORS_ONLN);

    if (nproc > nCPUs) {
        fprintf (stderr, "Error: nproc must not be greater than "
        "sysconf (_SC_NPROCESSORS_ONLN) = %d\n", nCPUs);
        return -1;
    }

    if (nproc > 1) { 
        semid = semget (IPC_PRIVATE, nproc - 1, 0600 | IPC_CREAT);
        if (semid == -1) {
            perror ("semget");
            exit (-1);
        }
    } else
        semid = 0;

    if ((net = (calc_t *) malloc (sizeT * sizeX * sizeof (*net))) == NULL) {
        perror ("malloc");
        return -1;
    }

    ThreadData *data = alloca (nproc * sizeof (*data));
    pthread_t *thr = alloca (nproc * sizeof (*thr));
    if (data == NULL || thr == NULL) {
        perror ("alloca");
        return -1;
    }

    for (int i = 0; i < nproc; i++) {
        calcBounds (i, nproc, &(data[i].left), &(data[i].right));
        data[i].id = i;
        int rc = pthread_create (thr + i, NULL, slave, data + i);
        if (rc != 0) {
            perror("pthread_create");
            return -1;
        }
    }

    for (int i = 0; i < nproc; i++)
        pthread_join (thr[i], NULL);

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
            printf ("%.3f ", net[k * sizeX + m]);
        printf ("\n");
    }


    cleanup:
    free (net);
}

void calcBounds (int rank, int comm_size, int *left, int *right) {
    double leftd  = (double)  rank      / comm_size * sizeX;
    double rightd = (double) (rank + 1) / comm_size * sizeX;
    *left  = -(int)(-leftd - 0.5); // round up
    *right = (int)(rightd + 0.5); // round down
}

int maxBound (int size_global, int comm_size) {
    return size_global / comm_size + 1;
}
