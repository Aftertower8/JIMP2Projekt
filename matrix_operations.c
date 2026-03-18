#include "matrix_operations.h"
#include <stdlib.h>
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

void matrix_cpy(double **dest, double **source, int size){
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

double** create_adjacency_matrix(graf g){
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

double** adjacency_to_laplacian_matrix(double** adj_matrix, int* deg_matrix, int size){    //no new matrix in order to save memory and time

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

void reverse_power_iteration(double **matrix, int size){
    double **m_cpy = (double*)malloc(sizeof(double*)*size);
    if(!m_cpy)
        return;
    matrix_cpy(m_cpy,matrix,size);
    
}