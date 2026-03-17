#include "matrix_operations.h"
#include <stdlib.h>
#define max(a,b) (((a)>(b)) ? (a) : (b))
int get_max_vertex(graf g){
    int max_vert = -1;
    for(int i=0; i<g.l_pkt; i++)
        max_vert = max(max(g.linki[i].a, g.linki[i].b), max_vert);
    return max_vert;
}

int** create_adjacency_matrix(graf g){
    int vertex_max = get_max_vertex(g);
    int **adj_matrix = (int**)malloc(sizeof(int*)*(vertex_max+1));
    if(!adj_matrix)
        return NULL;
    for(int i=0;i<vertex_max;i++){
        adj_matrix[i] = (int*)calloc(vertex_max, sizeof(int));
        if(!adj_matrix[i]){
            for(int j=0;j<i;j++){
                free(adj_matrix[j]);
                adj_matrix[j]=NULL;
            }
            free(adj_matrix);
            adj_matrix=NULL;
            return NULL;
        }
    }
    for(int i=0; i<vertex_max; i++){
        int a=g.linki[i].a;
        int b=g.linki[b].b;
        adj_matrix[a][b]=1;
        adj_matrix[b][a]=1;
    }
    return adj_matrix;
}

void free_adjacency_matrix(int **matrix, int max){
    for(int i=0;i<max;i++){
        free(matrix[i]);
        matrix[i]=NULL;
    }
    free(matrix);
    matrix=NULL;
}

int* create_degree_vector(int** adj_matrix, int max){
    int *degree_vector = (int*)calloc(max, sizeof(int));
    if(!degree_vector)
        return NULL;
    for(int i=0;i<max;i++){
        for(int j=0;j<max;j++){
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

int** create_laplacian_matrix(int** adj_matrix, int* deg_matrix, int max){

    for(int i=0;i<max;i++){
        for(int j=0;j<max;j++){
            if(i==j)
                adj_matrix[i][j] += deg_matrix[i];
            if(adj_matrix[i][j]==1)
                adj_matrix[i][j] = -1;
        }
    }
    return adj_matrix;
}