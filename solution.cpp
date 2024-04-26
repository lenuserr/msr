#include "inc.h"

// 20:10. 5 hour work. должен все отдебажить и закодить 100 проц

void print_matrix(size_t N, size_t* I, double* A) {
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

void print_vec(size_t N, double* b) {
    for (int i = 0; i < N; ++i) {
        printf("%lf ", b[i]);
    }
    printf("\n");
}

void* solution(void* ptr) {
    Args* args = (Args*)ptr;
    double a = args->a; double b = args->b; double c = args->c; double d = args->d;
    double eps = args->eps; size_t* I = args->I; double* A = args->A; double* B = args->B;
    double* x = args->x; double* r = args->r; double* u = args->u; double* v = args->v;
    [[maybe_unused]] int m = args->m; int maxit = args->maxit;
    int nx = args->nx; int ny = args->ny; int p = args->p; int k = args->k;

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;
    double s = hx*hy;
    size_t N = (nx + 1) * (size_t)(ny + 1);
    //double (*f)(double, double);    
    //select_f(m, f);

    fill_A(nx, ny, hx, hy, I, A, p, k);
    fill_B(nx, ny, hx, hy, a, c, B, f_0, p, k); 
    int maxsteps = 5;
    minimal_residual_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);    

    reduce_sum<int>(p);
    return nullptr;
}
