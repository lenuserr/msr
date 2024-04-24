#include <stddef.h>
#include "math.h"

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

#define F(I, J) (f(a + (I)*hx, с + (J)*hy))

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k);
int minimal_residual_msr_matrix(int n, double* A, int* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k);

int minimal_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k);

void thread_rows(int n, int p, int k, int& i1, int& i2);
double scalar_product(int n, double* x, double* y, int p, int k);
void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k);
void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v, double* r, int p, int k);
void ij2l(int nx, int ny, int i, int j, size_t& l);
void l2ij(int nx, int ny, int& i, int& j, size_t l);
size_t get_len_msr(int nx, int ny);
int get_off_diag(int nx, int ny, int i, int j, size_t* I_ij = nullptr);
size_t get_len_msr_off_diag(int nx, int ny);
int allocate_msr_matrix(int nx, int ny, double** p_A, size_t** p_I);
void fill_I(int nx, int ny, size_t* I);
void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag);
void fill_A(int nx, int ny, double hx, double hy, size_t* I, double* A, int p, int k);
int check_symm(int nx, int ny, size_t* I, double* A, int p, int k);
double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double));
void fill_b(int nx, int ny, double hx, double hy, double a, double c, double* b, double (*f)(double, double), int p, int k);
double f0(double, double);
double f1(double x, double);
double f2(double, double y);
double f3(double x, double y);
double f4(double x, double y);
double f5(double x, double y);
double f6(double x, double y);
double f7(double x, double y);

void* solution(void* ptr);
void reduce_sum(int p, int* err=nullptr, int n=0);
int init_reduce_sum(int p);
double reduce_sum_det(int p, int k, double s);

enum class Status {
    success,
    error
};

struct Args{
    double a;
    double b;
    double c;
    double d;
    double eps;
    size_t* I;
    double* A;
    int nx;
    int ny;
    int m;
    int maxit;
    int p;
    int k;
    Status status = Status::success;
};
