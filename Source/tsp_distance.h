#ifndef TSP_DISTANCE_H
#define TSP_DISTANCE_H
#include <stddef.h>
#include "tsp_types.h"

/* Type unique de pointeur de fonction pour les distances */
typedef double (*DistanceFn)(const TSPLIB_INSTANCE* I, int i, int j);

/* Sélecteur automatique selon type*/
DistanceFn distance_pour(const TSPLIB_INSTANCE* I);

/* Fonctions de distance */
double dist_euc2d(const TSPLIB_INSTANCE* I, int i, int j);
double dist_geo  (const TSPLIB_INSTANCE* I, int i, int j);
double dist_att  (const TSPLIB_INSTANCE* I, int i, int j);

/* Calcul de la longueur d'une tournée */
double longueur_tour(const TSPLIB_INSTANCE* I, const TOUR_TSP* T, DistanceFn d);

#endif