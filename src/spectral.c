#include "spectral.h"
#include "utlis.h"
#include "log.h"
#define ITERATIONS 1000

int power_iteration(Matrix *A,  Matrix *matrix, Vector *P,
                    Vector *help_vec, Vector *ones, Vector *Lap_w,
                    Vector **nxt, Vector **cur, Vector *v2,
                    double *lambda, double *lambda_prev
                    ){
    if(LU_solve(A, P, cur, nxt, help_vec)==-1){
        return -1;
    }
    vector_orthagonalization(nxt,ones);  //delete component v1
    if(v2)
        vector_orthagonalization(nxt,v2);  //delete component v2 if present
    if(vector_normalize(nxt)==-1){
        return -1;
    }
        
    for(int j=0; j<matrix->size;j++){
        //lambda2 = w * (Lap * w)/(w * w)
        //w*w == 1
        VEC(Lap_w,j) = 0.0;
        for(int k=0;k<matrix->size;k++)
            VEC(Lap_w,j) += MAT(matrix,j,k) * VEC(*nxt,k);
    }
    *lambda = scalar_product(nxt, Lap_w);
        
    if(is_zero(*lambda-*lambda_prev)){
        Vector *tmp = *cur;
        *cur = *nxt;
        *nxt = tmp;
        return 1;
    }
    *lambda_prev=*lambda;

    Vector *tmp = *cur;
    *cur = *nxt;
    *nxt = tmp;

    return 0;                    
}

void clear_memory(Matrix *A, int *P, Vector *cur2, Vector *nxt2,
                 Vector *cur3, Vector *nxt3, Vector *ones,
                 Vector *help_vec, Vector *Lap_w){
        free_matrix(A);
        free(P);
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
    if(v2)
        vector_orthagonalization(cur, v2);
    vector_normalize(cur);
}

int reverse_power_iterations(Matrix *matrix, pkt *punkty){
    int iteration_exit_value = 0; //after each power_iteration exit code from said function is assigned (error detection)
    int converged2 = 0; //flag if v2 converged
    int converged3 = 0; //flag if v3 converged
    double lambda2 = 0.0;
    double lambda2_prev = 1.0;
    double lambda3 = 0.0;
    double lambda3_prev = 1.0;
    int size = matrix->size;

    //first initialized as null as C standard allows free(NULL) (in case of clear_memory())
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
        clear_memory(A, P, cur2, nxt2, cur3, nxt3, ones, help_vec, Lap_w);
        return -1;
    }

    for(int i=0;i<size;i++)
        VEC(ones,i) = 1.0;
    
    prepare_current_vector(A,matrix,P,cur2,ones,NULL,0);

    for(int i=0;i<ITERATIONS;i++){
        iteration_exit_value = power_iteration(A, matrix, P, help_vec, ones, Lap_w, &nxt2, &cur2, NULL, &lambda2, &lambda2_prev);
        if(iteration_exit_value==-1){
            clear_memory(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w);
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
        iteration_exit_value = power_iteration(A,matrix,P,help_vec,ones,Lap_w,&nxt3,&cur3,v2,&lambda3,&lambda3_prev);
        if(iteration_exit_value==-1){
            clear_memory(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w);
            return -1;
        }
        if(iteration_exit_value==1){
            converged3 = 1;
            break;
        }
    }
    Vector *v3 = cur3;
    //x and y coordinates - P(x[i],y[i]) is point of i-vertex
    //copy coordinates
    for(int i=0; i<matrix->size; i++){
        punkty[i].x = VEC(v2,i);
        punkty[i].y = VEC(v3,i);
    }
    clear_memory(A,P,cur2,nxt2,cur3,nxt3,ones,help_vec,Lap_w);
    if(converged2==0 || converged3==0)
        return 1;
    return 0;
}

int SpectralLayoutAlgorithm(graf g, pkt *punkty){

    Matrix *M = create_adjacency_matrix(g);
    Vector *V = create_degree_vector(M);
    adjacency_to_laplacian_matrix(M,V);
    int exit_code = reverse_power_iterations(M,punkty);
    free_matrix(M);
    free_vec(V);
    if(exit_code==0){
        verbose("Spectral layout algorithm .\n");
        return 0;
    }
    else if(exit_code==1){
        fprintf(stderr, "Spectral layout did not converge. Results may be inacurate.\n");
        return 1;
    }
    else{
        fprintf(stderr, "Error in reverse_power_iterations function.\n");
        return -1;
    }
}