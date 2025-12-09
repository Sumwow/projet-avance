#ifndef TSP_2OPT_H
#define TSP_2OPT_H

#include "tsp_types.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"

double two_opt(const TSPLIB_INSTANCE* I, DistanceFn d, TOUR_TSP* tour);

#endif