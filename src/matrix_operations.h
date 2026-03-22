#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
#include "graf.h"

typedef struct{
    int size;
    double *data;
} Matrix;   //N x N matrix

typedef struct{
    int size;
    double *data;
} Vector;

static inline double getM(Matrix *m, int i, int j){
    return m->data[i * m->size + j];
}

static inline void setM(Matrix *m, int i, int j, double val){
    m->data[i * m->size + j] = val;
}

static inline double getV(Vector *v, int i){
    return v->data[i];
}

static inline void setV(Vector *v, int i, int val){
    return v->data[i] = val;
}


Matrix* create_adjacency_matrix(graf g);  //stworzy macierz sasiedztwa, dodac argumenty (ustalic strukture grafu) - ozn. A
int* create_degree_vector(double** adj_matrix, int size);  //ozn. D
double** adjacency_to_laplacian_matrix(double** adj_matrix, int* deg_matrix, int size);     //L = D - A 
double vector_norm(double *v, int n);   //norma wektora
int vector_normalize(double *v, int n);    //normalizacja wektora
double scalar_product(double *a, double *b, int n);
double squared_length(double *v, int n);
void scaling(double *v, int size, int scalar);
int vector_orthagonalization(double *result, double *component, int size);
//void power_iteration(double** matrix, double* vector, int n, int iterations);
#endif
