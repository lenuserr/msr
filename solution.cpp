#include "inc.h"

// не путай, пожалуйста, m и k
void* solution(void* ptr) {
    Args* args = (Args*)ptr;
    double a = args->a; double b = args->b; double c = args->c; double d = args->d;
    double eps = args->eps; size_t* I = args->I; double* A = args->A;
    int m = args->m; int maxit = args->maxit;
    int nx = args->nx; int ny = args->ny; int p = args->p; int k = args->k;

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    fill_A(nx, ny, hx, hy, I, A, p, k);

    return nullptr;
}
