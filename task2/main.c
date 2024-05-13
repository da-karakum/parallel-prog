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

const int MAX_NUM_THREADS = 4;
#define GLOBAL_STACK_SIZE 128

typedef struct {
    int num_threads;
    int thread_id;
} ThreadData;

typedef struct {
    double a, b, fa, fb, sab;
} Task;

struct Global_t {
    Task list_of_tasks [GLOBAL_STACK_SIZE];
    int ntask;
    int nactive;
    double total_sum;
    int semid;
    double precision;
} Global;

double FROM = 1e-5;
double TO = 1;
double integrand (double x) {
    return sin (1 / x);
}

enum sem_t {
    SEMAPHORE_LIST, 
    SEMAPHORE_TOTAL_SUM, 
    SEMAPHORE_TASK_PRESENT, 
    SEMAPHORE_COUNT
};

const int   SEMAPHORE_LOCK = -1;
const int SEMAPHORE_UNLOCK = +1;
int semop_one (int semid, enum sem_t sem_num, int sem_flg, int sem_op)
{
    struct sembuf oper = { 
        .sem_flg = sem_flg,
        .sem_num = sem_num,
        .sem_op = sem_op,
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

void usage (const char *argv0) {
    printf ("Usage: %s [num_threads] [0 < precision < 1]\n", argv0);
}

void *slave (void *context) {

    const int SPK = 16;
    const int LOCAL_STACK_SIZE = 128;
    Task local_list [LOCAL_STACK_SIZE];
    int local_ntask = 0;
    double local_sum = 0;
    ThreadData *thread_data = (ThreadData *)context;
    Task curr;

    set_affinity (thread_data->thread_id);

    while (1) {
        // Get a task from global stack
        semop_one (Global.semid, SEMAPHORE_TASK_PRESENT, 0, SEMAPHORE_LOCK);
        semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_LOCK);

        curr = Global.list_of_tasks[--Global.ntask];
        if (Global.ntask > 0)
            semop_one (Global.semid, SEMAPHORE_TASK_PRESENT, 0, 
                        SEMAPHORE_UNLOCK);
        
        if (curr.a <= curr.b) 
            Global.nactive++;
        
        semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_UNLOCK);

        // 
        if (curr.a > curr.b)
            break;
        
        // Integrate an interval from stack
        while (1) {
            double  c = (curr.a + curr.b) / 2.0;
            double fc = integrand (c);

            double sum_ac = (c      - curr.a) * (fc      + curr.fa) / 2.0;
            double sum_cb = (curr.b -      c) * (curr.fb +      fc) / 2.0;
            double sum_acb = sum_ac + sum_cb;

            double disbalance = fabs(sum_acb - curr.sab) / fabs (curr.sab);
            if (disbalance < Global.precision) {
                local_sum += sum_acb;
                if (local_ntask == 0)
                    break;

                curr = local_list[--local_ntask];
            } else {
                if (local_ntask >= LOCAL_STACK_SIZE) {
                    printf ("error: local stack overflow\n");
                    exit (-1);
                }
                // отрезок ac
                local_list[local_ntask].a = curr.a;
                local_list[local_ntask].b = c;
                local_list[local_ntask].fa = curr.fa;
                local_list[local_ntask].fb = fc;
                local_list[local_ntask].sab = sum_ac;
                local_ntask++;
                // отрезок cb
                curr.a = c;
                curr.fa = fc;
                curr.sab = sum_cb;
            }

            if (local_ntask > SPK && Global.ntask == 0) {
                semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_LOCK);

                if (Global.ntask == 0)
                    semop_one (Global.semid, SEMAPHORE_TASK_PRESENT, 0, SEMAPHORE_UNLOCK);

                while (local_ntask > 1) {
                    curr = local_list[--local_ntask];
                    if (Global.ntask >= GLOBAL_STACK_SIZE) {
                        printf ("error: global stack overflow\n");
                        exit (-1);
                    }
                    Global.list_of_tasks[Global.ntask++] = curr;
                }

                semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_UNLOCK);
            }
        }

        semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_LOCK);

        Global.nactive--;

        if (Global.nactive == 0 && Global.ntask == 0) {
            curr.a = 1; curr.b = 0;
            for (int i = 0; i < thread_data->num_threads; i++)
                Global.list_of_tasks[Global.ntask++] = curr;
            semop_one (Global.semid, SEMAPHORE_TASK_PRESENT, 0, SEMAPHORE_UNLOCK);
        }

        semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_UNLOCK);



    }

    semop_one (Global.semid, SEMAPHORE_TOTAL_SUM, 0, SEMAPHORE_LOCK);
    Global.total_sum += local_sum;
    semop_one (Global.semid, SEMAPHORE_TOTAL_SUM, 0, SEMAPHORE_UNLOCK);



    pthread_exit(NULL);
}
/*
argv[1] - num of threads
argv[2] - precision
*/
int main (int argc, char *argv[]) {

    if (argc != 3) {
        usage (argv[0]);
        exit (-1);
    }

    char *endch;
    int num_of_threads = strtol (argv[1], &endch, 10);
    if (*endch != ' ' && *endch != '\n' && *endch != '\0') {
        usage (argv[0]);
        exit (-1);
    } else if (num_of_threads > MAX_NUM_THREADS) {
        printf ("No more than %d threads\n", MAX_NUM_THREADS);
        exit (-1);
    }

    Global.precision = strtod (argv[2], &endch);
    if (*endch != ' ' && *endch != '\n' && *endch != '\0') {
        usage (argv[0]);
        exit (-1);
    } else if (Global.precision <= 0 || Global.precision >= 1) {
        printf ("0 < precision < 1\n");
        exit (-1);
    }

    Global.ntask = 0;
    Global.nactive = 0;
    Global.total_sum = 0;

    creat ("key", 0666);
    key_t key = ftok ("key", getpid());
    if (key == -1) {
        perror ("key");
        unlink ("key");
        exit (-1);
    }

    Global.semid = semget (key, SEMAPHORE_COUNT, 0600 | IPC_CREAT);
    if (Global.semid == -1) {
        perror ("semget");
        unlink ("key");
        exit (-1);
    }

    semop_one (Global.semid, SEMAPHORE_TOTAL_SUM, 0, SEMAPHORE_UNLOCK);
    semop_one (Global.semid, SEMAPHORE_LIST, 0, SEMAPHORE_UNLOCK);



    pthread_t threads[MAX_NUM_THREADS];
    ThreadData thread_data[MAX_NUM_THREADS];

    // Create threads
    for (int i = 0; i < num_of_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].num_threads = num_of_threads;

        int rc = pthread_create(&threads[i], NULL, slave, (void *)&thread_data[i]);

        if (rc != 0) {
            perror("pthread_create");
            exit(-1);
        }

    }


    Global.list_of_tasks[0].a = FROM;
    Global.list_of_tasks[0].b = TO;
    Global.list_of_tasks[0].fa = integrand(FROM);
    Global.list_of_tasks[0].fb = integrand(TO);
    Global.list_of_tasks[0].sab = (TO - FROM) * 
            (Global.list_of_tasks[0].fa + Global.list_of_tasks[0].fb) / 2;
    Global.ntask = 1;
    semop_one (Global.semid, SEMAPHORE_TASK_PRESENT, 0, SEMAPHORE_UNLOCK);

    for (int i = 0; i < num_of_threads; i++)
        pthread_join (threads[i], NULL);

    printf ("Answer: %.17f\n", Global.total_sum);
}
