#include "inc.h"

double f0(double, double) {
    return 1;
}

double f1(double x, double) {
    return x;
}

double f2(double, double y) {
    return y;
}

double f3(double x, double y) {
    return x + y;
}

double f4(double x, double y) {
    return sqrt(x*x + y*y);
}

double f5(double x, double y) {
    return x*x + y*y;
}

double f6(double x, double y) {
    return exp(x*x - y*y);
}

double f7(double x, double y) {
    return 1. / (25*(x*x + y*y) + 1);
}
