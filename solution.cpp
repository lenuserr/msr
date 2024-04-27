#include <sys/sysinfo.h>
#include "inc.h"

void print_matrix(int N, int* I, double* A) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                int m = I[i + 1] - I[i];
                int q;
                for (q = 0; q < m; ++q) {
                    if (I[I[i] + q] == j) {
                        break;
                    }
                }

                if (q >= m) {
                    printf("%lf ", 0.0);
                } else {
                    printf("%lf ", A[I[i] + q]);
                }
            } else {
                printf("%lf ", A[i]);
            }
        }
        printf("\n");
    }
}

void print_vec(int N, double* b) {
    for (int i = 0; i < N; ++i) {
        printf("%lf ", b[i]);
    }
    printf("\n");
}

void* solution(void* ptr) {
    Args* args = (Args*)ptr;
    double a = args->a; double b = args->b; double c = args->c; double d = args->d;
    double eps = args->eps; int* I = args->I; double* A = args->A; double* B = args->B;
    double* x = args->x; double* r = args->r; double* u = args->u; double* v = args->v;
    int maxit = args->maxit; int nx = args->nx; int ny = args->ny;
    int p = args->p; int k = args->k; double (*f)(double, double) = args->f;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;
    int N = (nx + 1) * (ny + 1);
    
    fill_A(nx, ny, hx, hy, I, A, p, k);
    fill_B(nx, ny, hx, hy, a, c, B, f, p, k); 

    int maxsteps = 500; // гиперпараметр
    minimal_residual_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);    
    
    reduce_sum<int>(p);
    return nullptr;
}
