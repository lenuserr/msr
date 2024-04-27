#include "inc.h"

double f_0(double, double) {
    return 1;
}

double f_1(double x, double) {
    return x;
}

double f_2(double, double y) {
    return y;
}

double f_3(double x, double y) {
    return x + y;
}

double f_4(double x, double y) {
    return sqrt(x*x + y*y);
}

double f_5(double x, double y) {
    return x*x + y*y;
}

double f_6(double x, double y) {
    return exp(x*x - y*y);
}

double f_7(double x, double y) {
    return 1. / (25*(x*x + y*y) + 1);
}

void Functions::select_f(int func_id) {
    switch (func_id) {
        case 0:
            f = f_0;
            break;
        case 1:
            f = f_1;
            break;
        case 2:
            f = f_2;
            break;
        case 3:
            f = f_3;
            break;
        case 4:
            f = f_4;
            break;
        case 5:
            f = f_5;
            break;
        case 6:
            f = f_6;
            break;
        case 7:
            f = f_7;
            break;
    }
}
