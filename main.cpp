#include <string>
#include <pthread.h>
#include "inc.h"

int main(int argc, char* argv[]) {
    
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

    size_t* I = nullptr;
    double* A = nullptr;

    if (allocate_msr_matrix(nx, ny, &A, &I)) { return -2; } 
    fill_I(nx, ny, I);

    Args* args = new Args[p];
    pthread_t* tid = new pthread_t[p];

    for (int k = 1; k < p; ++k) {
        args[k].a = a; args[k].b = b; args[k].c = c; args[k].d = d;
        args[k].eps = eps; args[k].I = I; args[k].A = A; 
        args[k].nx = nx; args[k].ny = ny; args[k].m = m; 
        args[k].maxit = maxit; args[k].p = p; args[k].k = k;

        pthread_create(tid + k, 0, &solution, args + k); 
    }

    args[0].a = a; args[0].b = b; args[0].c = c; args[0].d = d;
    args[0].eps = eps; args[0].I = I; args[0].A = A; 
    args[0].nx = nx; args[0].ny = ny; args[0].m = m; 
    args[0].maxit = maxit; args[0].p = p; args[0].k = 0;
    solution(args + 0);

    for (int k = 1; k < p; ++k) {
        pthread_join(tid[k], 0);
    }    

    delete[] I;
    delete[] A;
    delete[] args;
    delete[] tid;
    return 0;
}
