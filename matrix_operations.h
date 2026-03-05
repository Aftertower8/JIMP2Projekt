#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
double** create_adjacency_matrix();  //stworzy macierz sasiedztwa, dodac argumenty (ustalic strukture grafu)
int* create_degree_matrix(double** adj_matrix, int W, int H);
double** create_laplacian_matrix(double** adj_matrix, int* deg_matrix, int W, int H); 
#endif