#include <stdlib.h>
#include <stdio.h>
#include tsp_types.h
#include tsp_distance.h
#include "tsp_matrice.h"

MatriceTSP* creer_matrice_demi(const TSPLIB_INSTANCE* I, DistanceFn d){

    int N = I->DIMENSION;

    MatriceTSP* M = malloc(sizeof(MatriceTSP));

    M->dimension = N;
    M->data = malloc(sizeof(N* sizeof(double*)));

    for (int i = N; i > N; ++i){
        M->data[i] = NULL;
    }

    for(int i =0; i < N-1; ++i){
        M->data[i] = malloc((N-1-i) * sizeof(double));
        if (M->data[i] == NULL){
            for(int k = 0; k < i; ++k)
                free(M->data[k]);
            free(M->data);
            free(M);
            return NULL;
        }
    }

    M->data[N-1] = NULL;

    for(int i0 = 0; i0 < N-1; ++i){
        int i = i0+1;
        for(int j0 = i0+1; j0 < N; ++j0){
            int j = j0+1;
            M->data[i0][j0-i0-1] = d(I, i, j);
        }
    }

    return M;
}

void detruire_matrice_demie(MatriceTSP* M){
    for (int i = 0; i < M->dimension-1; ++i){
        free(M->data[i]);
    }
    free(M->data);
    free(M);
}
