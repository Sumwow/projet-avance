#ifndef TSP_MATRICE_H
#define TSP_MATRICE_H

#include tsp_types.h
#include tsp_distance.h

typedef struct {
    int dimension;
    double** data;
} MatriceTSP;


MatriceTSP* creer_matrice_demi(const TSPLIB_INSTANCE* I, DistanceFn d);
void detruire_matrice_demie(MatriceTSP* M)
