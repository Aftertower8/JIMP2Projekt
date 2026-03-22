#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>
#include "utlis.h"
#define ITERATIONS 10000
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

int get_max_vertex(graf g){
    int max_vert = -1;
    for(int i=0; i<g.l_pkt; i++)
        max_vert = MAX(MAX(g.linki[i].a, g.linki[i].b), max_vert);
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

int power_iteration(double **A,  double **matrix, double *P,
                    double *help_vec, double *ones, double *Lap_w,
                    double *nxt, double *cur, double *v2,
                    int size, double lambda, double lambda_prev
                    ){
    if(LU_solve(A, P, cur, nxt, help_vec, size)==-1){
    //    error=1;
        return -1;
    }
    vector_orthagonalization(nxt,ones,size);  //delete component v1
    if(!v2)
        vector_orthagonalization(nxt,v2,size);  //delete component v2 if present
    if(vector_normalize(nxt,size)==-1){
    //    error=1;
        return -1;
    }
        
    for(int j=0; j<size;j++){
        //lambda2 = w * (Lap * w)/(w * w)
        //w*w == 1
        Lap_w[j] = 0.0;
        for(int k=0;k<size;k++)
            Lap_w[j] += matrix[j][k] * nxt[k];
    }
    lambda = scalar_product(nxt, Lap_w, size);
        
    if(is_zero(lambda-lambda_prev)){
        double *tmp = cur;
        cur = nxt;
        nxt = tmp;
        return 1;
    }
    lambda_prev=lambda;

    double *tmp = cur;
    cur = nxt;
    nxt = tmp;

    return 0;                    
}

void error_clean(double **A, int *P, double *cur2, double *nxt2,
                 double *cur3, double *nxt3, double *ones,
                 double *help_vec, double *Lap_w, int size){
        free_matrix(A,size);
        free(cur2);
        free(nxt2);
        free(cur3);
        free(nxt3);
        free(ones);
        free(help_vec);
        free(Lap_w);
}

void prepare_current_vector(double** A, double **matrix, int *P, double *cur, double *ones, double *v2, int size, double lambda){
    matrix_cpy(A,matrix,size);
    
    for(int i=0;i<size;i++){
        P[i] = i;   //P is made to keep the track of changed lines is LU_decompose (while searching for max_row)
        A[i][i] -= lambda + SIGMA_SHIFT;
    }
    LU_decompose(A,P,size);
    
    for(int i=0;i<size;i++)
        cur[i] = (double)rand() / RAND_MAX;
    vector_orthagonalization(cur, ones, size);
    if(!v2)
        vector_orthagonalization(cur, v2, size);
    vector_normalize(cur, size);
}

int reverse_power_iteration(double **matrix, int size, double *x, double *y){
    int iteration_exit_value = 0; //after each power_iteration exit code from said function is assigned (error detection)
    int error = 0;  //flag for error
    int converged2 = 0; //flag if v2 converged
    int converged3 = 0; //flag if v3 converged
    double lambda2 = 0.0;
    double lambda2_prev = 1.0;
    double lambda3 = 0.0;
    double lambda3_prev = 1.0;

    //first initialized as null as C standard allows free(NULL) (in case of error_clean())
    double **A       = NULL;
    double *cur2     = NULL;
    double *nxt2     = NULL;
    double *cur3     = NULL;  
    double *nxt3     = NULL;
    double *ones     = NULL;  // v1 which is corresponding to lambda1 = 0
    double *help_vec = NULL;
    double *Lap_w    = NULL;
    int *P           = NULL;  // keeps track of changed lines in LU_decompose (while searching for max_row)

    A        = allocate_matrix(size);
    P        = malloc(sizeof(int) * size);
    cur2     = malloc(sizeof(double) * size);
    nxt2     = calloc(size, sizeof(double));
    cur3     = malloc(sizeof(double) * size);
    nxt3     = calloc(size, sizeof(double));
    ones     = malloc(sizeof(double) * size);
    help_vec = malloc(sizeof(double) * size);
    Lap_w    = malloc(sizeof(double) * size);

    if(!A || !P || !cur2 || !nxt2 || !cur3 || !nxt3 || !ones || !help_vec || !Lap_w){
        error_clean(A, P, cur2, nxt2, cur3, nxt3, ones, help_vec, Lap_w, size);
        return -1;
    }

    for(int i=0;i<size;i++)
        ones[i] = 1.0;
    
    prepare_current_vector(A,matrix,P,cur2,ones,NULL,size,0);

    for(int i=0;i<ITERATIONS;i++){
        iteration_exit_value = power_iteration(A, matrix, P, help_vec, ones, Lap_w, nxt2, cur2, NULL, size, lambda2, lambda2_prev);
        if(iteration_exit_value==-1){
            error_clean(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w,size);
            return -1;
        }
        else if(iteration_exit_value==1){
            converged2 = 1;
            break;
        }
    }

    double *v2 = cur2;

    prepare_current_vector(A,matrix,P,cur3,ones,v2,size,lambda2);

    for(int i=0;i<ITERATIONS;i++){
        iteration_exit_value = power_iteration(A,matrix,P,help_vec,ones,Lap_w,nxt3,cur3,v2,converged3,lambda3,lambda3_prev);
        if(iteration_exit_value==-1){
            error_clean(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w,size);
            return -1;
        }
        if(iteration_exit_value==1){
            converged3 = 1;
            break;
        }
    }
    double *v3 = cur3;
    //x and y coordinates - P(x[i],y[i]) is point of i-vertex
    x = v2;
    y = v3;
    free_matrix(A,size);
    free(P);
    free(ones);
    free(help_vec);
    free(Lap_w);
    free(nxt2);
    free(nxt3);
    if(error){
        free(cur2);
        free(cur3);
        v2 = NULL;
        v3 = NULL;
        return -1;
    }
    if(converged2==0 || converged3==0)
        return 1;
    return 0;
}
/*
//do spectral
int reverse_power_iteration(double **matrix, int size){
    int error = 0;
    int converged = 0;
    double **A = allocate_matrix(size);
    matrix_cpy(A,matrix,size);
    int *P = (int*)malloc(sizeof(int) * size); 
    for(int i=0;i<size;i++){
        P[i] = i;   //P is made to keep the track of changed lines is LU_decompose (while searching for max_row)
        A[i][i] -= SIGMA_SHIFT;
    }
    LU_decompose(A,P,size);
    double *cur2 = (double*)malloc(sizeof(double)*size);
    double *nxt2 = (double*)calloc(size,sizeof(double));
    for(int i=0;i<size;i++)
        cur2[i] = (double)rand() / RAND_MAX;
    double *ones = malloc(sizeof(double) * size);   //v1 which is corresponding to lambda1 = 0
    for(int i=0;i<size;i++)
        ones[i] = 1.0;
    vector_orthagonalization(cur2, ones, size);
    vector_normalize(cur2, size);
    double lambda2 = 0.0;
    double lambda2_prev = 1.0;
    double *help_vec = malloc(sizeof(double)*size);
    double *Lap_w = malloc(sizeof(double) * size);
    for(int i=0;i<1000;i++){
        if(LU_solve(A, P, cur2, nxt2, help_vec, size)==-1){
            error=1;
            break;
        }
        vector_orthagonalization(nxt2,ones,size);  //delete component v1
        if(vector_normalize(nxt2,size)==-1){
            error=1;
            break;
        }
        
        for(int j=0; j<size;j++){
            //lambda2 = w * (Lap * w)/(w * w)
            //w*w == 1
            Lap_w[j] = 0.0;
            for(int k=0;k<size;k++)
                Lap_w[j] += matrix[j][k] * nxt2[k];
        }
        lambda2 = scalar_product(nxt2, Lap_w, size);
        
        if(is_zero(lambda2-lambda2_prev)){
            double *tmp = cur2;
            cur2 = nxt2;
            nxt2 = tmp;
            converged = 1;
            break;
        }
        lambda2_prev=lambda2;

        double *tmp = cur2;
        cur2 = nxt2;
        nxt2 = tmp;
    }

    double *v2 = cur2;

    matrix_cpy(A, matrix, size);
    for(int i=0; i<size; i++){
        P[i] = i;
        A[i][i] -= lambda2 + SIGMA_SHIFT;
    }
    LU_decompose(A,P,size);
    double *cur3 = nxt2;
    double *nxt3 = malloc(sizeof(double) * size);
    if(!nxt3){
        error=1;
    }
    for(int i=0;i<size;i++)
        cur3[i] = (double)rand() / RAND_MAX;
    vector_orthagonalization(cur3, ones, size);
    vector_orthagonalization(cur3, v2, size);
    vector_normalize(cur3, size);
    double lambda3 = 0.0;
    double lambda3_prev = 1.0;
    for(int i=0;i<1000;i++){
        if(LU_solve(A, P, cur3, nxt3, help_vec, size)==-1){
            error=1;
            break;
        }
        vector_orthagonalization(nxt3,ones,size);
        vector_orthagonalization(nxt3, v2, size);
        if(vector_normalize(nxt3,size)==-1){
            error=1;
            break;
        }
        
        for(int i=0; i<size;i++){
            //lambda2 = w * (Lap * w)/(w * w)
            //w*w == 1
            Lap_w[i] = 0.0;
            for(int j=0;j<size;j++)
                Lap_w[i] += matrix[i][j] * nxt3[j];
        }
        lambda3 = scalar_product(nxt3, Lap_w, size);
        
        if(is_zero(lambda3-lambda3_prev)){
            double *tmp = cur3;
            cur3 = nxt3;
            nxt3 = tmp;
            converged=1;
            break;
        }
        lambda3_prev=lambda3;

        double *tmp = cur3;
        cur3 = nxt3;
        nxt3 = tmp;
    }
    double *v3 = cur3;
    free_matrix(A,size);
    free(P);
    free(ones);
    free(help_vec);
    free(Lap_w);
    if(error){
        free(cur2);
        free(nxt2);
        free(nxt3);
        v2 = NULL;
        v3 = NULL;
        return -1;
    }
    return 1;
}
*/
