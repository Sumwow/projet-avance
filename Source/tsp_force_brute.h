#ifndef TSP_FORCE_BRUTE_H
#define TSP_FORCE_BRUTE_H

#include "tsp_distance.h"  /* TSPLIB_INSTANCE, TOUR_TSP, DistanceFn */

/* Recherche exhaustive (DistanceFn).
   Retourne la longueur minimale (double) ou -1.0 si erreur. */
double force_brute(const TSPLIB_INSTANCE* instance,
                   DistanceFn calculDistance,
                   TOUR_TSP* meilleureTournee,
                   TOUR_TSP* pireTournee);

#endif /* TSP_FORCE_BRUTE_H */
