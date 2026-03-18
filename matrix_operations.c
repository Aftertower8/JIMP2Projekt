#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>
#define max(a,b) (((a)>(b)) ? (a) : (b))\

const double eps = 1e-6;

int get_max_vertex(graf g){
    int max_vert = -1;
    for(int i=0; i<g.l_pkt; i++)
        max_vert = max(max(g.linki[i].a, g.linki[i].b), max_vert);
    return max_vert;
}

void error_alocating(double **matrix, int i){
    for(int j=0;j<i;j++){
        free(matrix[j]);
    }
    free(matrix);
}

void matrix_cphelp_vec(double **dest, double **source, int size){
    for(int i=0;i<size; i++){
        dest[i] = (double*)malloc(sizeof(double) * (size));
        if(!dest[i]){
            error_alocating(dest,i);
            dest=NULL;
            return;
        }
        for(int j=0;j<size;j++)
            dest[i][j] = source[i][j];
    }
}

double** create_adjacenchelp_vec_matrix(graf g){
    int size = get_max_vertex(g)+1;
    int **adj_matrix = (int**)malloc(sizeof(int*)*size);
    if(!adj_matrix)
        return NULL;
    for(int i=0;i<size;i++){
        adj_matrix[i] = (int*)calloc(size, sizeof(int));
        if(!adj_matrix[i]){
            error_alocating(adj_matrix,i);
            adj_matrix=NULL;
            return NULL;
        }
    }
    for(int i=0; i<size; i++){
        int a=g.linki[i].a;
        int b=g.linki[b].b;
        adj_matrix[a][b]=1.0;
        adj_matrix[b][a]=1.0;
    }
    return adj_matrix;
}

void free_matrix(double **matrix, int size){
    for(int i=0;i<size;i++){
        free(matrix[i]);
        matrix[i]=NULL;
    }
    free(matrix);
    matrix=NULL;
}

int* create_degree_vector(double** adj_matrix, int size){
    int *degree_vector = (int*)calloc(size, sizeof(int));
    if(!degree_vector)
        return NULL;
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            if(adj_matrix[i][j]==1)
                degree_vector[i]++;
        }
    }
    return degree_vector;
}

void free_degree_vector(int *vec){
    free(vec);
    vec=NULL;
}

double** adjacenchelp_vec_to_laplacian_matrix(double** adj_matrix, int* deg_matrix, int size){    //no new matrix in order to save memorhelp_vec and time

    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            if(i==j)
                adj_matrix[i][j] += (double)deg_matrix[i];
            else if(adj_matrix[i][j] == 1.0)
                adj_matrix[i][j] = -1.0;
        }
    }
    return adj_matrix;
}
//LU matrix -> A = L * U
//in order to save memorhelp_vec, L and U matrix are stored in A matrix
void LU_decompose(double **A, int *P, int size){    //implementation of Gaussian  elimination in order to get LU matrix
    //searching for max_row in order to minimalize rounding errors
    for(int k=0; k<size; k++){
        int max_row = k;
        for(int i = k+1; i<size; i++)
            if(fabs(A[i][k]) > fabs(A[max_row][k]))
                max_row = i;
        int tmp = P[k];
        P[k] = P[max_row];
        P[max_row] = tmp;
        double *tmp_row = A[k];
        A[k] = A[max_row];
        A[max_row] = tmp_row;

        if(fabs(A[k][k]) < 1e-14)               //check if A[k][k] ~ 0
            continue;
        for(int i=k+1; i<size; i++){            //Gaussian elimination - building LU matrix
            double factor = A[i][k] / A[k][k];  //A[i][k] elimination factor
            A[i][k] = factor;                   //to have L and U in one matrix store factor in eliminated position                  
            for(int j = k+1; j < size; j++)
                A[i][j] -= factor * A[k][j];
        }
    }
}

void LU_solve(double **A, int *P, double *current, double *res, int size){  //A * w = v
    //U - upper
    //L - lower
    //A * res = current
    //L * U * res = current
    //helper: help_vec = U * res
    //L * help_vec = current
    //U * res = help_vec
    double *help_vec = (double*)malloc(sizeof(double)*size);
    for(int i=0; i<size; i++){
        help_vec[i] = current[P[i]];                         //current[P[i]] to get valid index after changing lines in LU_decompose
        for(int j=0; j<i; j++){                 //L * help_vec = v
            help_vec[i] -= (A[i][j] * help_vec[j]);
        }
    //no dividing as A[i][i] contains diagonal of U matrix
    }
    for(int i=size-1; i>=0; i--){               //U * w = help_vec
        res[i] = help_vec[i];
        for(int j=i+1; i<size; j++){
            res[i] -= (A[i][j] * res[j]);
        }
        res[i] /= A[i][i];
    }
    free(help_vec);
}

void reverse_power_iteration(double **matrix, int size){
    double **A = (double*)malloc(sizeof(double*)*size);
    if(!A)
        return;
    matrix_cphelp_vec(A,matrix,size);
    int *P = (int*)malloc(sizeof(int) * size); 
    for(int i=0;i<size;i++){
        P[i] = i;   //P is made to keep the track of changed lines is LU_decompose (while searching for max_row)
        A[i][i] -= eps;
    }
    LU_decompose(A,P,size);
    double *v = (double*)malloc(sizeof(double)*size);
    double *w = (double*)calloc(size,sizeof(double));
    for(int i=0;i<size;i++)
        v[i] = (double)rand() / RAND_MAX;
    
}