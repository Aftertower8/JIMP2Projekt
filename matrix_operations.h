#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
double** create_adjacency_matrix();  //stworzy macierz sasiedztwa, dodac argumenty (ustalic strukture grafu)
int* create_degree_vector(double** adj_matrix, int W, int H);
double** create_laplacian_matrix(double** adj_matrix, int* deg_matrix, int W, int H);   
double vector_norm(double *v, int n);   //norma wektora
void vector_normalize(double *v, int n);    //normalizacja wektora
double vector_scalar(double *a, double *b, int n);  //mnozenie skalarne wektorow
void multiply_matrix_by_vector(double** matrix, double *vector, double *result, int n);
void power_iteration(double** matrix, double* vector, int n, int iterations);
#endif