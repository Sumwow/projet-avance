#ifndef TSP_MATRICE_H
#define TSP_MATRICE_H

#include "tsp_types.h"
#include "tsp_distance.h"

/* Structure représentant une matrice triangulaire compacte */
typedef struct {
    int dimension;   /* nombre de villes */
    double** data;   /* ligne i (0-based) contient N-1-i distances */
} MatriceTSP;

/* Création et remplissage de la matrice triangulaire */
MatriceTSP* creer_matrice_demie(const TSPLIB_INSTANCE* I, DistanceFn d);

/* Libération mémoire */
void detruire_matrice_demie(MatriceTSP* M);

/* Accès symétrique : indices TSPLIB 1..N */
static inline double matrice_distance(const MatriceTSP* M, int i, int j) {
    if (!M || !M->data || i == j) return 0.0;
    /* Si j > i, on est dans la partie stockée directement */
    if (j > i)  return M->data[i - 1][j - i - 1];
    /* Sinon, on inverse les indices (symétrie du TSP) */
    else        return M->data[j - 1][i - j - 1];
}

#endif

