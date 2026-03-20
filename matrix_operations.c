#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>
#define max(a,b) (((a)>(b)) ? (a) : (b))
//#define MATRIX(m,i,j) (m)->data[(i)*(m)->size +(j)]

/*
    TODO:
    -poprawic komentarze
    -wprowadzic strukture na macierze i wektory
    -optymalizacja
    -byc moze zmiana metody
    -zmiana nazw zmiennych
    -sprawdzenie poprawnosci
    -przerzucenie elementow do spectral.c i .h
    -obsluga bledow
*/

const double sigma = 1e-6;
const double eps = 1e-9;

int is_zero(double a){
    return fabs(a) < eps;
}

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
        /*
        dest[i] = (double*)malloc(sizeof(double) * (size));
        if(!dest[i]){
            error_alocating(dest,i);
            dest=NULL;
            return;
        }
        */
        for(int j=0;j<size;j++)
            dest[i][j] = source[i][j];
    }
}

double** allocate_matrix(int size){
    double **matrix = (double**)malloc(sizeof(double*)*size);
    if(!matrix)
        return NULL;
    for(int i=0;i<size;i++){
        matrix[i] = (double*)calloc(size, sizeof(double));
        if(!matrix[i]){
            error_alocating(matrix,i);
            matrix=NULL;
            return NULL;
        }
    }
    return matrix;
}

double** create_adjacency_matrix(graf g){
    int size = get_max_vertex(g)+1;
    double **adj_matrix = allocate_matrix(size);
    for(int i=0; i<g.l_l; i++){
        int a=g.linki[i].a;
        int b=g.linki[i].b;
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

double** adjacency_to_laplacian_matrix(double** adj_matrix, int* deg_matrix, int size){    //no new matrix in order to save memorhelp_vec and time

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
//in order to save memory, L and U matrix are stored in A matrix
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

        if(is_zero(A[k][k]))               //check if A[k][k] ~ 0
            continue;
        for(int i=k+1; i<size; i++){            //Gaussian elimination - building LU matrix
            double factor = A[i][k] / A[k][k];  //A[i][k] elimination factor
            A[i][k] = factor;                   //to have L and U in one matrix store factor in eliminated position                  
            for(int j = k+1; j < size; j++)
                A[i][j] -= factor * A[k][j];
        }
    }
}

int LU_solve(double **A, int *P, double *current, double *res, double *help_vec, int size){  //A * w = v
    //U - upper
    //L - lower
    //A * res = current
    //L * U * res = current
    //helper: help_vec = U * res
    //L * help_vec = current
    //U * res = help_vec
    for(int i=0; i<size; i++){
        help_vec[i] = current[P[i]];                         //current[P[i]] to get valid index after changing lines in LU_decompose
        for(int j=0; j<i; j++){                 //L * help_vec = v
            help_vec[i] -= (A[i][j] * help_vec[j]);
        }
    //no dividing as A[i][i] contains diagonal of U matrix
    }
    for(int i=size-1; i>=0; i--){               //U * w = help_vec
        res[i] = help_vec[i];
        for(int j=i+1; j<size; j++){
            res[i] -= (A[i][j] * res[j]);
        }
        if(is_zero(A[i][i]))
            return -1;
        res[i] /= A[i][i];
    }
    return 0;
}

double scalar_product(double *a, double *b, int n){
    double scalar=0;
    for(int i=0;i<n;i++)
        scalar += a[i]*b[i];
    return scalar;
}

double squared_length(double *v, int n){
    double s_len=0;
    for(int i=0;i<n;i++)
        s_len += pow(v[i],2);
    return s_len;
}

void scaling(double *v, int size, int scalar){
    for(int i=0;i<size;i++)
        v[i] *= scalar;
}

int vector_orthagonalization(double *result, double *component, int size){
    double scalar = scalar_product(result,component,size);
    double squared_len = squared_length(component,size);
    if(fabs(squared_len) < 1e-14)
        return -1;
    double coeff = scalar / squared_len;
    for(int i=0;i<size; i++)
        result[i] -= coeff * component[i];
    return 0;
}

double vector_norm(double *v, int n){
    return sqrt(squared_length(v,n));
}

int vector_normalize(double *v, int n){
    double norm = vector_norm(v,n);
    if(is_zero(norm))
        return -1;
    for(int i=0;i<n;i++)
        v[i] /= norm;
    return 0;
}

//do spectral
int reverse_power_iteration(double **matrix, int size){
    int error = 0;
    int converged = 0;
    double **A = allocate_matrix(size);
    matrix_cpy(A,matrix,size);
    int *P = (int*)malloc(sizeof(int) * size); 
    for(int i=0;i<size;i++){
        P[i] = i;   //P is made to keep the track of changed lines is LU_decompose (while searching for max_row)
        A[i][i] -= sigma;
    }
    LU_decompose(A,P,size);
    double *v = (double*)malloc(sizeof(double)*size);
    double *w = (double*)calloc(size,sizeof(double));
    for(int i=0;i<size;i++)
        v[i] = (double)rand() / RAND_MAX;
    double *ones = malloc(sizeof(double) * size);   //v1 which is corresponding to lambda1 = 0
    for(int i=0;i<size;i++)
        ones[i] = 1.0;
    vector_orthagonalization(v, ones, size);
    vector_normalize(v, size);
    double lambda2 = 0.0;
    double lambda2_prev = 1.0;
    double *help_vec = (double*)malloc(sizeof(double)*size);
    double *Lap_w = malloc(sizeof(double) * size);
    for(int i=0;i<1000;i++){
        if(LU_solve(A, P, v, w, help_vec, size)==-1){
            error=1;
            break;
        }
        vector_orthagonalization(w,ones,size);  //delete component v1
        if(vector_normalize(w,size)==-1){
            error=1;
            break;
        }
        
        for(int j=0; j<size;j++){
            //lambda2 = w * (Lap * w)/(w * w)
            //w*w == 1
            Lap_w[j] = 0.0;
            for(int k=0;k<size;k++)
                Lap_w[j] += matrix[j][k] * w[j];
        }
        lambda2 = scalar_product(w, Lap_w, size);

        double *tmp = v;
        v = w;
        w = tmp;
        
        if(is_zero(lambda2-lambda2_prev)){
            converged = 1;
            break;
        }
        lambda2_prev=lambda2;
    }

    double *v2 = v;
    double sigma3 = lambda2 + sigma;

    matrix_cpy(A, matrix, size);
    for(int i=0; i<size; i++){
        P[i] = i;
        A[i][i] -= sigma3;
    }
    LU_decompose(A,P,size);
    double *current = w;
    double *nxt = malloc(sizeof(double) * size);
    if(!nxt){
        error=1;
    }
    for(int i=0;i<size;i++)
        current[i] = (double)rand() / RAND_MAX;
    vector_orthagonalization(current, ones, size);
    vector_orthagonalization(current, v2, size);
    vector_normalize(current, size);
    double lambda3 = 0.0;
    double lambda3_prev = 1.0;
    for(int i=0;i<1000;i++){
        if(LU_solve(A, P, current, nxt, help_vec, size)==-1){
            error=1;
            break;
        }
        vector_orthagonalization(nxt,ones,size);
        vector_orthagonalization(nxt, v2, size);
        if(vector_normalize(nxt,size)==-1){
            error=1;
            break;
        }
        
        for(int i=0; i<size;i++){
            //lambda2 = w * (Lap * w)/(w * w)
            //w*w == 1
            Lap_w[i] = 0.0;
            for(int j=0;j<size;j++)
                Lap_w[i] += matrix[i][j] * nxt[j];
        }
        lambda3 = scalar_product(nxt, Lap_w, size);

        double *tmp = current;
        current = nxt;
        nxt = tmp;
        
        if(is_zero(lambda3-lambda3_prev)){
            converged=1;
            break;
        }
        lambda3_prev=lambda3;
    }
    double *v3 = current;
    free_matrix(A,size);
    free(P);
    free(ones);
    free(help_vec);
    free(Lap_w);
    if(error){
        free(v);
        free(w);
        free(nxt);
        v2 = NULL;
        v3 = NULL;
        return -1;
    }
    return 1;
}