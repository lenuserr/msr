#include "inc.h"
#include <algorithm>

double r3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    int N = (nx+1) * (ny+1);
    int i1, i2;
    int i, j;
    thread_rows(N, p, k, i1, i2);
    double res = -1;
    for (int l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        res = std::max(res, fabs(f(a + i*hx, c + j*hy) - x[l]));
    }

    reduce_sum(p, &res, 1, &max);    
    return res;
}
