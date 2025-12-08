#ifndef TSP_FORCE_BRUTE_H
#define TSP_FORCE_BRUTE_H

#include "tsp_distance.h"

double force_brute(const TSPLIB_INSTANCE* instance,
                   DistanceFn calculDistance,
                   TOUR_TSP* meilleureTournee,
                   TOUR_TSP* pireTournee);

#endif
