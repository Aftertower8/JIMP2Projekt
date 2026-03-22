#include "matrix_operations.h"
#include <stdlib.h>
#include <math.h>
#include "utlis.h"
#define ITERATIONS 10000
#define MAT(m,i,j) ((m)->data[(i) * (m)->size + (j)])
#define VEC(v,i) ((v)->data[(i)])
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

void error_alocating(double **matrix, int i){   //niepotrzebne
    for(int j=0;j<i;j++){
        free(matrix[j]);
    }
    free(matrix);
}

int matrix_cpy(Matrix *dest, Matrix *source){
    if(!dest || !source)
        return -1;
    for(int i=0;i<dest->size * dest->size; i++)
        source->data[i] = dest->data[i];
    return 0;
}

Matrix *allocate_matrix(int size){
    Matrix *M = malloc(sizeof(Matrix));
    if(!M)
        return NULL;
    M->size = size;
    M->data = calloc(size * size,sizeof(double));
    if(!M->data){
        free(M);
        return NULL;
    }
    return M;
}

Vector *allocate_vector(int size){
    Vector *V = malloc(sizeof(Vector));
    if(!V)
        return NULL;
    V->size = size;
    V->data = calloc(size,sizeof(double));
    if(!V->data){
        free(V);
        return NULL;
    }
    return V;
}

Matrix *create_adjacency_matrix(graf g){
    int size = get_max_vertex(g)+1;
    Matrix *adj = allocate_matrix(size);
    for(int i=0; i<g.l_l; i++){
        int a=g.linki[i].a;
        int b=g.linki[i].b;
        MAT(adj,a,b) = 1.0;
        MAT(adj,b,a) = 1.0;
    }
    return adj;
}

void free_matrix(Matrix *M){
    free(M->data);
    free(M);
}

void free_vec(Vector *V){
    free(V->data);
    fre(V);
}

Vector* create_degree_vector(Matrix *M){
    Vector *degree_vector = malloc(sizeof(Vector));
    if(!degree_vector)
        return NULL;
    degree_vector->size = M->size;
    degree_vector->data = calloc(degree_vector->size, sizeof(int));
    if(!degree_vector->data){
        free(degree_vector);
        return NULL;
    }
    
    for(int i=0;i<M->size;i++){
        for(int j=0;j<M->size;j++){
            if(MAT(M,i,j)==1)
                degree_vector->data[i]++;
        }
    }
    return degree_vector;
}

void free_degree_vector(int *vec){  //niepotrzebne
    free(vec);
    vec=NULL;
}

Matrix* adjacency_to_laplacian_matrix(Matrix* adj_matrix, Vector* deg_matrix){    //no new matrix in order to save memorhelp_vec and time

    for(int i=0;i<adj_matrix->size;i++){
        for(int j=0;j<adj_matrix->size;j++){
            if(i==j)
                MAT(adj_matrix,i,j) += deg_matrix->data[i];
            else if(MAT(adj_matrix,i,j) == 1.0)
                MAT(adj_matrix,i,j) = -1.0;
        }
    }
    return adj_matrix;
}
//LU matrix -> A = L * U
//in order to save memory, L and U matrix are stored in A matrix
void LU_decompose(Matrix *A, int *P){    //implementation of Gaussian  elimination in order to get LU matrix
    double *tmp_row = malloc(sizeof(double) * A->size);
    //searching for max_row in order to minimalize rounding errors
    for(int k=0; k<A->size; k++){
        int max_row = k;
        for(int i = k+1; i<A->size; i++)
            if(fabs(MAT(A,i,k)) > fabs(MAT(A,max_row,k)))
                max_row = i;
        if(max_row != k){
            int tmp = P[k];
            P[k] = P[max_row];
            P[max_row] = tmp;

            for(int j=0; j<A->size; j++){
                tmp_row[j] = MAT(A, k,j);
                MAT(A,k,j) = MAT(A, max_row,j);
                MAT(A,max_row,j) = tmp_row[j];
            }

        }
        if(is_zero(MAT(A,k,k)))               //check if A[k][k] ~ 0
            continue;
        for(int i=k+1; i<A->size; i++){            //Gaussian elimination - building LU matrix
            double factor = MAT(A,i,k) / MAT(A,k,k);  //A[i][k] elimination factor
            MAT(A,i,k) = factor;                   //to have L and U in one matrix store factor in eliminated position                  
            for(int j = k+1; j < A->size; j++)
                MAT(A,i,j) -= factor * MAT(A,k,j);
        }
    }
}


int LU_solve(Matrix *A, int *P, Vector *current, Vector *res, Vector *help_vec){  //A * w = v
    //U - upper
    //L - lower
    //A * res = current
    //L * U * res = current
    //helper: help_vec = U * res
    //L * help_vec = current
    //U * res = help_vec
    for(int i=0; i<A->size; i++){
        VEC(help_vec, i) = VEC(current, P[i]); //VEC(current, P[i]) to get valid index after changing lines in LU_decompose
        for(int j=0; j<i; j++){                 //L * help_vec = v
            VEC(help_vec, i) -= MAT(A, i, j) * VEC(help_vec, j);
        }
    //no dividing as A[i][i] contains diagonal of U matrix
    }
    for(int i=A->size-1; i>=0; i--){               //U * w = help_vec
        VEC(res, i) = VEC(help_vec, i);
        for(int j=i+1; j<A->size; j++){
            VEC(res, i) -= MAT(A, i, j) * VEC(res, j);
        }
        if(is_zero(MAT(A, i, i)))
            return -1;
        VEC(res, i) /= MAT(A, i, i);
    }
    return 0;
}

double scalar_product(Vector *a, Vector *b){
    double scalar=0;
    for(int i=0;i<a->size;i++)
        scalar += VEC(a,i)*VEC(b,i);
    return scalar;
}

double squared_length(Vector *v){
    double s_len=0;
    for(int i=0;i<v->size;i++)
        s_len += pow(VEC(v,i),2);
    return s_len;
}

void scaling(Vector *v, double scalar){
    for(int i=0;i<v->size;i++)
        VEC(v,i) *= scalar;
}

int vector_orthagonalization(Vector *result, Vector *component){
    double scalar = scalar_product(result,component);
    double squared_len = squared_length(component);
    if(fabs(squared_len) < EPS)
        return -1;
    double coeff = scalar / squared_len;
    for(int i=0;i<result->size; i++)
        VEC(result,i) -= coeff * VEC(component,i);
    return 0;
}

double vector_norm(Vector *v){
    return sqrt(squared_length(v));
}

int vector_normalize(Vector *v){
    double norm = vector_norm(v);
    if(is_zero(norm))
        return -1;
    for(int i=0;i<v->size;i++)
        VEC(v,i) /= norm;
    return 0;
}

int power_iteration(Matrix *A,  Matrix *matrix, Vector *P,
                    Vector *help_vec, Vector *ones, Vector *Lap_w,
                    Vector *nxt, Vector *cur, Vector *v2,
                    double lambda, double lambda_prev
                    ){
    if(LU_solve(A, P, cur, nxt, help_vec)==-1){
        return -1;
    }
    vector_orthagonalization(nxt,ones);  //delete component v1
    if(!v2)
        vector_orthagonalization(nxt,v2);  //delete component v2 if present
    if(vector_normalize(nxt)==-1){
        return -1;
    }
        
    for(int j=0; j<matrix->size;j++){
        //lambda2 = w * (Lap * w)/(w * w)
        //w*w == 1
        VEC(Lap_w,j) = 0.0;
        for(int k=0;k<matrix->size;k++)
            VEC(Lap_w,j) += MAT(matrix,j,k) * VEC(nxt,k);
    }
    lambda = scalar_product(nxt, Lap_w);
        
    if(is_zero(lambda-lambda_prev)){
        Vector *tmp = cur;
        cur = nxt;
        nxt = tmp;
        return 1;
    }
    lambda_prev=lambda;

    Vector *tmp = cur;
    cur = nxt;
    nxt = tmp;

    return 0;                    
}

void error_clean(Matrix *A, int *P, Vector *cur2, Vector *nxt2,
                 Vector *cur3, Vector *nxt3, Vector *ones,
                 Vector *help_vec, Vector *Lap_w){
        free_matrix(A);
        free_vec(cur2);
        free_vec(nxt2);
        free_vec(cur3);
        free_vec(nxt3);
        free_vec(ones);
        free_vec(help_vec);
        free_vec(Lap_w);
}

void prepare_current_vector(Matrix* A, Matrix *matrix, int *P, Vector *cur, Vector *ones, Vector *v2, double lambda){
    matrix_cpy(A,matrix);
    
    for(int i=0;i<A->size;i++){
        P[i] = i;   //P is made to keep the track of changed lines is LU_decompose (while searching for max_row)
        MAT(A,i,i) -= lambda + SIGMA_SHIFT;
    }
    LU_decompose(A,P);
    
    for(int i=0;i<A->size;i++)
        VEC(cur,i) = (double)rand() / RAND_MAX;
    vector_orthagonalization(cur, ones);
    if(!v2)
        vector_orthagonalization(cur, v2);
    vector_normalize(cur);
}

int reverse_power_iteration(Matrix *matrix, Vector *x, Vector *y){
    int iteration_exit_value = 0; //after each power_iteration exit code from said function is assigned (error detection)
    int error = 0;  //flag for error
    int converged2 = 0; //flag if v2 converged
    int converged3 = 0; //flag if v3 converged
    double lambda2 = 0.0;
    double lambda2_prev = 1.0;
    double lambda3 = 0.0;
    double lambda3_prev = 1.0;
    int size = matrix->size;

    //first initialized as null as C standard allows free(NULL) (in case of error_clean())
    Matrix *A        = NULL;
    Vector *cur2     = NULL;
    Vector *nxt2     = NULL;
    Vector *cur3     = NULL;  
    Vector *nxt3     = NULL;
    Vector *ones     = NULL;  // v1 which is corresponding to lambda1 = 0
    Vector *help_vec = NULL;
    Vector *Lap_w    = NULL;
    int *P           = NULL;  // keeps track of changed lines in LU_decompose (while searching for max_row)

    A        = allocate_matrix(size);
    P        = malloc(sizeof(int) * size);
    cur2     = allocate_vector(size);
    nxt2     = allocate_vector(size);
    cur3     = allocate_vector(size);
    nxt3     = allocate_vector(size);
    ones     = allocate_vector(size);
    help_vec = allocate_vector(size);
    Lap_w    = allocate_vector(size);

    if(!A || !P || !cur2 || !nxt2 || !cur3 || !nxt3 || !ones || !help_vec || !Lap_w){
        error_clean(A, P, cur2, nxt2, cur3, nxt3, ones, help_vec, Lap_w);
        return -1;
    }

    for(int i=0;i<size;i++)
        VEC(ones,i) = 1.0;
    
    prepare_current_vector(A,matrix,P,cur2,ones,NULL,0);

    for(int i=0;i<ITERATIONS;i++){
        iteration_exit_value = power_iteration(A, matrix, P, help_vec, ones, Lap_w, nxt2, cur2, NULL, lambda2, lambda2_prev);
        if(iteration_exit_value==-1){
            error_clean(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w);
            return -1;
        }
        else if(iteration_exit_value==1){
            converged2 = 1;
            break;
        }
    }

    Vector *v2 = cur2;

    prepare_current_vector(A,matrix,P,cur3,ones,v2,lambda2);

    for(int i=0;i<ITERATIONS;i++){
        iteration_exit_value = power_iteration(A,matrix,P,help_vec,ones,Lap_w,nxt3,cur3,v2,lambda3,lambda3_prev);
        if(iteration_exit_value==-1){
            error_clean(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w);
            return -1;
        }
        if(iteration_exit_value==1){
            converged3 = 1;
            break;
        }
    }
    Vector *v3 = cur3;
    //x and y coordinates - P(x[i],y[i]) is point of i-vertex
    x = v2;
    y = v3;
    free_matrix(A);
    free(P);
    free_vec(ones);
    free_vec(help_vec);
    free_vec(Lap_w);
    free_vec(nxt2);
    free_vec(nxt3);
    if(error){
        free_vec(cur2);
        free_vec(cur3);
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
