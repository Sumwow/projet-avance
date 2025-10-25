#include <stdlib.h>
#include <stdio.h>
#include "tsp_matrice.h"

MatriceTSP* creer_matrice_demie(const TSPLIB_INSTANCE* I, DistanceFn d){
    if (I == NULL || d == NULL) return NULL;

    int N = I->DIMENSION;
    if (N <= 0) return NULL;

    MatriceTSP* M = (MatriceTSP*)malloc(sizeof(MatriceTSP));
    if (M == NULL) return NULL;

    M->dimension = N;
    M->data = (double**)malloc((size_t)N * sizeof(double*));
    if (M->data == NULL) {
        free(M);
        return NULL;
    }

    /* Initialisation à NULL */
    for (int i = 0; i < N; ++i) M->data[i] = NULL;

    /* Allocation des lignes triangulaires */
    for (int i = 0; i < N - 1; ++i) {
        M->data[i] = (double*)malloc((size_t)(N - 1 - i) * sizeof(double));
        if (M->data[i] == NULL) {
            for (int k = 0; k < i; ++k) free(M->data[k]);
            free(M->data);
            free(M);
            return NULL;
        }
    }
    M->data[N - 1] = NULL;

    /* Remplissage */
    for (int i0 = 0; i0 < N - 1; ++i0) {
        int i = i0 + 1;
        for (int j0 = i0 + 1; j0 < N; ++j0) {
            int j = j0 + 1;
            M->data[i0][j0 - i0 - 1] = d(I, i, j);
        }
    }

    return M;
}

void detruire_matrice_demie(MatriceTSP* M){
    if (M == NULL) return;
    for (int i = 0; i < M->dimension - 1; ++i)
        free(M->data[i]);
    free(M->data);
    free(M);
}

