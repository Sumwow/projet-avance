#include <math.h>
#include <string.h>
#include "tsp_distance.h"

/*acces au nœud (ville) numéro idx1 d’une instance TSPLIB*/
static inline const NODE* node_at(const TSPLIB_INSTANCE* I, int idx1) {
    return &I->NODE_COORD_SECTION[idx1 - 1];
}


/* Distance euclidienne arrondie à l'entier le plus proche */
double dist_euc2d(const TSPLIB_INSTANCE* I, int i, int j) {
    if (i == j) return 0.0;
    const NODE* a = node_at(I, i);
    const NODE* b = node_at(I, j);
    double dx = a->X - b->X;
    double dy = a->Y - b->Y;
    double r = sqrt(dx * dx + dy * dy);
    return floor(r + 0.5);
}


/* convertit une coordonnée géographique en radiant */
static inline double geo_to_rad(double dddmm) {
    int deg = (int)dddmm;
    double min = dddmm - (double)deg;
    double val = (double)deg + 5.0 * min / 3.0;
    return M_PI * val / 180.0;
}

/* Distance géographique TSPLIB */
double dist_geo(const TSPLIB_INSTANCE* I, int i, int j) {
    if (i == j) return 0.0;

    const NODE* a = node_at(I, i);
    const NODE* b = node_at(I, j);

    double lati = geo_to_rad(a->X);
    double loni = geo_to_rad(a->Y);
    double latj = geo_to_rad(b->X);
    double lonj = geo_to_rad(b->Y);

    double q1 = cos(loni - lonj);
    double q2 = cos(lati - latj);
    double q3 = cos(lati + latj);

    double R = 6378.388;
    double r = acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3));
    return floor(R * r + 1.0);
}


/* Didienne (TSPLIB) */
double dist_att(const TSPLIB_INSTANCE* I, int i, int j) {
    if (i == j) return 0.0;
    const NODE* a = node_at(I, i);
    const NODE* b = node_at(I, j);
    double dx = a->X - b->X;
    double dy = a->Y - b->Y;
    double rij = sqrt((dx * dx + dy * dy) / 10.0);
    double tij = floor(rij + 0.5);
    if (tij < rij)
        tij += 1.0;
    return tij;
}

/* Choisi en la fonction a appliquer en fonction du type*/
DistanceFn distance_pour(const TSPLIB_INSTANCE* I) {
    if (strcmp(I->EDGE_WEIGHT_TYPE, "ATT") == 0)
        return dist_att;
    if (strcmp(I->EDGE_WEIGHT_TYPE, "EUC_2D") == 0 ||
        strcmp(I->EDGE_WEIGHT_TYPE, "EUCL_2D") == 0)
        return dist_euc2d;
    if (strcmp(I->EDGE_WEIGHT_TYPE, "GEO") == 0)
        return dist_geo;
    return dist_euc2d; /* par défaut */
}

/* Calcul longueur Tourner */
double longueur_tour(const TSPLIB_INSTANCE* I, const TOUR_TSP* T, DistanceFn d) {
    if (!I || !T || T->DIMENSION <= 0) return 0.0;

    double L = 0.0;

    for (int k = 0; k < T->DIMENSION - 1; ++k) {
        int a = T->SECTION_TOUR[k];
        int b = T->SECTION_TOUR[k + 1];
        L += d(I, a, b);
    }

    if (T->FERMEE && T->DIMENSION > 1) {
        int last = T->SECTION_TOUR[T->DIMENSION - 1];
        int first = T->SECTION_TOUR[0];
        L += d(I, last, first);
    }

    return L;
}
