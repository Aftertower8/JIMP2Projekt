#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_
#include "graf.h"

#define MAT(m,i,j) ((m)->data[(i) * (m)->size + (j)])
#define VEC(v,i) ((v)->data[(i)])

typedef struct{
    int size;
    double *data;
} Matrix;   //N x N matrix

typedef struct{
    int size;
    double *data;
} Vector;

Matrix* allocate_matrix(int size);
Vector* allocate_vector(int size);
void    free_matrix(Matrix *M);
void    free_vec(Vector *V);
Matrix* create_adjacency_matrix(graf g);
Vector* create_degree_vector(Matrix *adj_matrix);
void adjacency_to_laplacian_matrix(Matrix *adj_matrix, Vector *deg_matrix);
void LU_decompose(Matrix *A, int *P);
int  LU_solve(Matrix *A, int *P, Vector *current, Vector *res, Vector *help_vec);
double vector_norm(Vector *v);
int    vector_normalize(Vector *v);
double scalar_product(Vector *a, Vector *b);
double squared_length(Vector *v);
void   scaling(Vector *v, double scalar);      /* zmiana: int -> double */
int    vector_orthagonalization(Vector *result, Vector *component);

#endif
