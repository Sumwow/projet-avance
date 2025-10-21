#ifndef TSP_FORCE_BRUTE_MATRICE_H
#define TSP_FORCE_BRUTE_MATRICE_H

#include "tsp_matrice.h"   /* MatriceTSP, TOUR_TSP */

/* Recherche exhaustive (demi-matrice).
   Retourne la longueur minimale (double) ou -1.0 si erreur. */
double force_brute_matrice(const MatriceTSP* matrice,
                           TOUR_TSP* meilleureTournee,
                           TOUR_TSP* pireTournee);

#endif /* TSP_FORCE_BRUTE_MATRICE_H */
