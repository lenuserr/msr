#include <string>
#include <cstring>
#include "inc.h"

void print_vec(int N, double* b) {
    for (int i = 0; i < N; ++i) {
        printf("%lf ", b[i]);
    }
    printf("\n");
}

int main(int argc, char* argv[]) {
    const int task = 3;    
    if (argc != 11) {
        return -1;
    }

    double a = std::stod(argv[1]);
    double b = std::stod(argv[2]);
    double c = std::stod(argv[3]);
    double d = std::stod(argv[4]);
    int nx = std::stoi(argv[5]);
    int ny = std::stoi(argv[6]);
    int m = std::stoi(argv[7]);
    double eps = std::stod(argv[8]);
    int maxit = std::stoi(argv[9]);
    int p = std::stoi(argv[10]);
    
    int* I = nullptr;
    double* A = nullptr;
    if (allocate_msr_matrix(nx, ny, &A, &I)) { return -2; }
    init_reduce_sum(p);
    int N = (nx + 1) * (ny + 1);
    double* B = new double[N];
    double* x = new double[N];
    double* r = new double[N];
    double* u = new double[N];
    double* v = new double[N];

    fill_I(nx, ny, I);

    memset(x, 0, N*sizeof(double));

    Functions F = Functions();
    F.select_f(m);
    double (*f)(double, double) = F.f;

    Args* args = new Args[p];
    pthread_t* tid = new pthread_t[p];
        
    for (int k = 1; k < p; ++k) {
        args[k].a = a; args[k].b = b; args[k].c = c; args[k].d = d;
        args[k].eps = eps; args[k].I = I; args[k].A = A; args[k].B = B; 
        args[k].x = x; args[k].r = r; args[k].u = u; args[k].v = v;
        args[k].nx = nx; args[k].ny = ny; args[k].maxit = maxit;
        args[k].p = p; args[k].k = k; args[k].f = f;

        pthread_create(tid + k, 0, &solution, args + k); 
    }

    args[0].a = a; args[0].b = b; args[0].c = c; args[0].d = d;
    args[0].eps = eps; args[0].I = I; args[0].A = A; args[0].B = B;
    args[0].x = x; args[0].r = r; args[0].u = u; args[0].v = v;
    args[0].nx = nx; args[0].ny = ny; args[0].maxit = maxit; 
    args[0].p = p; args[0].k = 0; args[0].f = f;
    solution(args + 0);

    for (int k = 1; k < p; ++k) {
        pthread_join(tid[k], 0);
    }    

    int its = args[0].its;
    double res1 = args[0].res_1; double res2 = args[0].res_2;
    double res3 = args[0].res_3; double res4 = args[0].res_4;
    double t1 = args[0].t1; double t2 = args[0].t2;

    printf (
    "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n\
    	  It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
          argv[0], task, res1, res2, res3, res4, t1, t2, its, eps, m, nx, ny, p);

    //print_vec(N, x);
    //print_vec(N, B);    

    free_results();
    delete[] I;
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] r;
    delete[] u;
    delete[] v;
    delete[] args;
    delete[] tid;

    return 0;
}
