#ifndef TSP_BRUTE_GENERIQUE_H
#define TSP_BRUTE_GENERIQUE_H

#include "tsp_matrice.h"

void brute_set_matrice(MatriceTSP *m);

double brute(int nb_nodes,
             int nb_ressources,
             int *best_perm,
             unsigned long long *count_best,
             void * (*cout)(void *, int *));

void *cout_tsp_matrice(void *unused, int *perm);

#endif