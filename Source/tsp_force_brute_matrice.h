#ifndef TSP_FORCE_BRUTE_MATRICE_H
#define TSP_FORCE_BRUTE_MATRICE_H

#include "tsp_matrice.h"

double force_brute_matrice(const MatriceTSP* matrice,
                           TOUR_TSP* meilleureTournee,
                           TOUR_TSP* pireTournee);

#endif