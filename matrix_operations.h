#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
#include "graf.h"
int** create_adjacency_matrix(graf g);  //stworzy macierz sasiedztwa, dodac argumenty (ustalic strukture grafu) - ozn. A
int* create_degree_vector(int** adj_matrix, int size);  //ozn. D
int** create_laplacian_matrix(int** adj_matrix, int* deg_matrix, int size);     //L = D - A 
double vector_norm(double *v, int n);   //norma wektora
void vector_normalize(double *v, int n);    //normalizacja wektora
double vector_scalar(double *a, double *b, int n);  //mnozenie skalarne wektorow
void multiply_matrix_by_vector(double** matrix, double *vector, double *result, int n);
void power_iteration(double** matrix, double* vector, int n, int iterations);
#endif