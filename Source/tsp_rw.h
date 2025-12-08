#ifndef TSP_RW_H
#define TSP_RW_H

#include "tsp_types.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"

double marche_aleatoire(const TSPLIB_INSTANCE* I, DistanceFn d, TOUR_TSP* tour);

#endif