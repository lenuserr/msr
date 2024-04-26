#include "inc.h"

// ключевое слово "странно считает" для поиска подозрителных функций
// 20:00. Пишу код пока с Якоби. С релаксацией разобрался в теории. update: пока дебажу, что чексум запустить хотя б

// Еще одну функцию написать для проверки того , что сумма по строкам матрицы определенная, а не рандомная.
// Богачев объяснял какая должна получится в конце 8го занятия
// надеюсь из-за того, что у него v и r похожи я их нигде не перепутал)

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

#define F(I, J) (f(a + (I)*hx, с + (J)*hy))

void matrix_mult_vector_msr(int n, double* A, size_t* I, double* x, double* y, int p, int k) {
    int i, i1, i2, l, J; double s;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s = A[i] * x[i];
        l = I[i+1] - I[i];
        J = I[i];
        for (int j = 0; j < l; ++j) {
            s += A[J + j] * x[I[J + j]];
        }

        y[i] = s;
    }
}

// Ax = b, x - в начале x0, в конце решение
int minimal_residual_msr_matrix(int n, double* A, size_t* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k) {

    double prec, b_norm2, tau, c1, c2;
    int it;
    b_norm2 = scalar_product(n, b, b, p, k); // (b, b)
    prec = b_norm2 * eps * eps;
    matrix_mult_vector_msr(n, A, I, x, r, p, k); // r = Ax
    mult_sub_vector(n, r, b, 1., p, k); // r -= 1*b

    for (it = 0; it < maxit; ++it) {
        // v = M^(-1) r - метод Якоби
        apply_preconditioner_msr_matrix(n, A, I, v, r, p, k);
        matrix_mult_vector_msr(n, A, I, v, u, p, k); // u = Av
        c1 = scalar_product(n, u, r, p, k); // C1 = (u, r)
        c2 = scalar_product(n, u, u, p, k); // C2 = (u, u)
        if (c1 < prec || c2 < prec) {
            break;
        }

        tau = c1 / c2;
        mult_sub_vector(n, x, v, tau, p, k); // x -= tau * v
        mult_sub_vector(n, r, u, tau, p, k); // r -= tau * u
    }

    if (it > maxit) {
        return -1;
    }

    return it;
}

int minimal_residual_msr_matrix_full(int n, double* A, size_t* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k) {

    int step, ret, its = 0;
    for (step = 0; step < maxsteps; ++step) {
        ret = minimal_residual_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if (ret >= 0) { its += ret; break; }
        its += maxit;
    }

    if (step >= maxsteps) {
        return -1;
    } 

    return its;
}

void thread_rows(int n, int p, int k, int& i1, int& i2) {
    i1 = n*k; i1 /= p; i2 = n*(k + 1); i2 /= p; 
}

double scalar_product(int n, double* x, double* y, int p, int k) {
    int i1, i2, i; double s = 0;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s += x[i] * y[i];
    }

    reduce_sum_det(p, k, s);
    return s;
}

void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        x[i] -= tau * y[i];
    }

    reduce_sum<int>(p);
}

void apply_preconditioner_msr_matrix(int n, double* A, [[maybe_unused]] size_t* I, double* v, double* r, int p, int k) {
    // v = M^{-1} * r
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        v[i] = r[i] / A[i]; 
    }    

    reduce_sum<int>(p);
}

void ij2l(int nx, [[maybe_unused]] int ny, int i, int j, size_t& l) {
    l = i + j * (size_t)(nx + 1);
}

void l2ij(int nx, [[maybe_unused]] int ny, int& i, int& j, size_t l) {
    j = l / (nx + 1);
    i = l - j * (size_t)(nx + 1);
}

size_t get_len_msr(int nx, int ny) {
    return (nx + 1)*(ny + 1)
            + 6*(nx - 1)*(ny - 1)
            + 4*(2*(nx-1) + 2*(ny-1))
            + 2*3 + 2*2; 
}

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

// странно считает 0 1 0 1 2 2 1 1e-10 1000 1
int get_off_diag(int nx, int ny, int i, int j, size_t* I_ij) {
    int m = 0; size_t k;
    if (i < nx) { FUNC(i+1, j); }
    if (j > 0) { FUNC(i, j-1); }
    if (i > 0 && j > 0) { FUNC(i-1, j-1); }
    if (i > 0) { FUNC(i - 1, j); }
    if (j < ny) { FUNC(i, j+1); }
    if (i < nx && j < ny) { FUNC(i+1, j+1); }

    return m;
}

size_t get_len_msr_off_diag(int nx, int ny) {
    size_t m = 0; int i, j;
    for (i = 0; i <= nx; ++i) {
        for (j = 0; j <= ny; ++j) {
            m += get_off_diag(nx, ny, i, j);
        }
    }    

    return m;
}

int allocate_msr_matrix(int nx, int ny, double** p_A, size_t** p_I) {
    size_t diag_len = (nx+1)*(ny+1);
    size_t off_diag = get_len_msr_off_diag(nx, ny);
    size_t len = diag_len + off_diag + 1;

    double* A = nullptr; size_t* I = nullptr;
    A = new double[len]; if (A == nullptr) { return 1; }
    I = new size_t[len]; if (I == nullptr) { return 2; }
    *p_A = A; *p_I = I;
    return 0;  
}

/*
void l2ij(int nx, [[maybe_unused]] int ny, int& i, int& j, size_t l) {
    j = l / (nx + 1);
    i = l - j * (size_t)(nx + 1);
}
*/

void fill_I(int nx, int ny, size_t* I) { 
    size_t N = (size_t)(nx+1) * (ny+1);
    size_t l; int i, j;
    int m; size_t r = N + 1;
    for (l = 0; l < N; ++l) {
        l2ij(nx, ny, i, j, l);
        I[l] = r;
        m = get_off_diag(nx, ny, i, j, I + r);
        r += m;
    }

    I[l] = r;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag) {

    double s = hx*hy;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        *A_diag = s / 2;
        for (int l = 0; l < 6; ++l) {
            A_off_diag[l] = s / 12; 
        }
    } else if (j == 0 && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (j == ny && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 2*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 1*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (i == nx && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j == 0) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
    } else if (i == nx && j == ny) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 1*s/24;
    } else if (i == 0 && j == ny) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    } else if (i == nx && j == 0) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    }

    // abort() если сюда дошли
}

void fill_A(int nx, int ny, double hx, double hy, size_t* I, double* A, int p, int k) {
    size_t l1, l2;
    int i, j;
    size_t N = (nx + 1) * (size_t)(ny + 1);    
    l1 = N * k; l1 /= p;
    l2 = N * (k + 1); l2 /= p;

    for (size_t l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        double* A_diag = A + l;
        double* A_off_diag = A + I[l];

        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_off_diag);
    }

    reduce_sum<int>(p);
}

// обязательно проверить матрицу построенную на симметричность иначе потом долго ошибку искать
int check_symm(int nx, int ny, size_t* I, double* A, double eps, int p, int k) {
    size_t l1, l2, l;
    size_t N = (nx+1)*(size_t)(ny+1);
    l1 = k * N; l1 /= p;
    l2 = (k + 1) * N; l2 /= p;

    int err = 0;

    for (l = l1; l < l2; ++l) {
        int m = I[l + 1] - I[l];
        double* A_off_diag = A + I[l];
        for (int q = 0; q < m; ++q) {
            double a_ij = A_off_diag[q];
            int j = I[I[l] + q];
            // найти в j-ой строке внедиагональный элемент с номером l
            int m2 = I[j + 1] - I[j];
            int q2;
            for (q2 = 0; q2 < m2; ++q2) {
                if (I[I[j] + q2] == l) {
                    break;
                }
            }

            if (q2 >= m2) {
                err++;
            } else if (fabs(A[I[j] + q2] - a_ij) > eps) {
                err++;
            } 
        }
    }

    reduce_sum<int>(p, &err, 1);
    return err;
}

double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double)) {
    double w = hx*hy/192;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        return w * (36*F(i, j) + 20*(F(i+0.5, j) + F(i, j-0.5) + F(i-0.5,j-0.5) + F(i-0.5, j) + F(i, j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5, j-0.5) + F(i-0.5, j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5))
            + 2*(F(i+1,j) + F(i, j-1) + F(i-1,j-1) + F(i-1,j) + F(i,j+1) + F(i+1,j+1))
        );
    }
    
    if (i > 0 && i < nx && j == 0) {
        return w * ( 
            18*F(i,j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i,j+1) + F(i+1,j+1))
        );
    }

    if (i > 0 && i < nx && j == ny) {
        return w * (
            18*F(i, j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j-0.5) + F(i-0.5,j-0.5))
            + 4*(F(i+0.5,j-0.5) + F(i-0.5,j-1) + F(i-1,j-0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i-1,j) + F(i-1,j-1))
        );
    }

    if (i == 0 && j > 0 && j < ny) {
        return w*(
            18*F(i,j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i+0.5,j) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5,j-0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i+1,j) + F(i+1,j+1))
        );     
    }

    if (i == nx && j > 0 && j < ny) {
        return w * (
            18*F(i, j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i-0.5,j) + F(i-0.5,j-0.5))
            + 4*(F(i-0.5,j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i-1,j) + F(i-1,j-1))
        );
    }

    if (i == 0 && j == 0) {
        return w*(
            12*F(i,j) + 10*(F(i+0.5,j) + F(i,j+0.5)) + 20*F(i+0.5,j+0.5)
            + 4*(F(i+1,j+0.5) + F(i+0.5,j+1)) + 1*(F(i+1,j) + F(i,j+1)) + 2*F(i+1,j+1)
        );     
    }    

    if (i == nx && j == ny) {
        return w*(
            12*F(i,j) + 10*(F(i-0.5,j) + F(i,j-0.5)) + 20*F(i-0.5,j-0.5)
            + 4*(F(i-0.5,j+1) + F(i-1,j-0.5)) + 1*(F(i,j-1) + F(i-1,j)) + 2*F(i-1,j-1)
        );   
    }

    if (i == 0 && j == ny) {
        return w*(
            6*F(i, j) + 10*(F(i+0.5,j) + F(i,j-0.5)) + 4*F(i+0.5,j-0.5) + F(i+1,j) + F(i,j-1) 
        );
    }

    if (i == nx && j == 0) {
        return w*(
            6*F(i, j) + 10*(F(i-0.5,j) + F(i,j+0.5)) + 4*F(i-0.5,j+0.5) + F(i-1,j) + F(i,j+1) 
        );      
    }

    return 1e308;
}

void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k) {
    size_t l1, l2;
    int i, j;
    size_t N = (nx + 1) * (size_t)(ny + 1);    
    l1 = N * k; l1 /= p;
    l2 = N * (k + 1); l2 /= p;

    for (size_t l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        B[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f);
    }

    reduce_sum<int>(p);
}
