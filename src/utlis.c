#include "utlis.h"
#include <math.h>

void swap(double **a, double **b){
    void *tmp = *a;
    *a = *b;
    *b = tmp;
}

int is_zero(double a){
    return fabs(a) < EPS;
}

