#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
#include "graf.h"
double** create_adjacency_matrix(graf g);  //stworzy macierz sasiedztwa, dodac argumenty (ustalic strukture grafu) - ozn. A
int* create_degree_vector(double** adj_matrix, int size);  //ozn. D
double** adjacency_to_laplacian_matrix(double** adj_matrix, int* deg_matrix, int size);     //L = D - A 
double vector_norm(double *v, int n);   //norma wektora
int vector_normalize(double *v, int n);    //normalizacja wektora
double scalar_product(double *a, double *b, int n);
double squared_length(double *v, int n);
void scaling(double *v, int size, int scalar);
int vector_orthagonalization(double *result, double *component, int size);
void multiply_matrix_by_vector(double** matrix, double *vector, double *result, int n);
void power_iteration(double** matrix, double* vector, int n, int iterations);
#endif