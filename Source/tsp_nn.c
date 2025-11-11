#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include "tsp_nn.h"

double plus_proche_voisin(const TSPLIB_INSTANCE* I,
                          DistanceFn d,
                          TOUR_TSP* tour,
                          int depart)
{
    if (!I || I->DIMENSION <= 0 || !d || !tour) {
        fprintf(stderr, "NN(coord): arguments invalides\n");
        return -1.0;
    }
    const int N = I->DIMENSION;
    if (!tour->SECTION_TOUR || tour->DIMENSION != N) {
        fprintf(stderr, "NN(coord): tour non allouee ou dimension incoherente\n");
        return -1.0;
    }
    if (depart <= 0 || depart > N) depart = 1;

    int* visite = (int*)calloc((size_t)(N + 1), sizeof(int));
    if (!visite) {
        fprintf(stderr, "NN(coord): allocation visite echouee\n");
        return -1.0;
    }

    int courant = depart;
    tour->SECTION_TOUR[0] = courant;
    visite[courant] = 1;

    double longueur = 0.0;

    for (int k = 1; k < N; ++k) {
        int meilleur = -1;
        int bestDist = INT_MAX;

        for (int j = 1; j <= N; ++j) {
            if (!visite[j]) {
                int dij = d(I, courant, j);
                if (dij < bestDist) {
                    bestDist = dij;
                    meilleur = j;
                }
            }
        }

        if (meilleur <= 0) {
            free(visite);
            fprintf(stderr, "NN(coord): impossible de trouver un voisin\n");
            return -1.0;
        }

        tour->SECTION_TOUR[k] = meilleur;
        visite[meilleur] = 1;
        longueur += (double)bestDist;
        courant = meilleur;
    }

    if (tour->FERMEE && N > 1) {
        int retour = d(I, courant, depart);
        longueur += (double)retour;
    }

    tour->LONGUEUR = longueur;
    free(visite);
    return longueur;
}

/*Version DemiMatrice (pas obligatoire je vous laisse le choix de la garder ou non) */

double plus_proche_voisin_matrice(const MatriceTSP* M,
                                  TOUR_TSP* tour,
                                  int depart)
{
    if (!M || !M->data || M->dimension <= 0 || !tour) {
        fprintf(stderr, "NN(matrice): arguments invalides\n");
        return -1.0;
    }
    const int N = M->dimension;
    if (!tour->SECTION_TOUR || tour->DIMENSION != N) {
        fprintf(stderr, "NN(matrice): tour non allouee ou dimension incoherente\n");
        return -1.0;
    }
    if (depart <= 0 || depart > N) depart = 1;

    int* visite = (int*)calloc((size_t)(N + 1), sizeof(int));
    if (!visite) {
        fprintf(stderr, "NN(matrice): allocation visite echouee\n");
        return -1.0;
    }

    int courant = depart;
    tour->SECTION_TOUR[0] = courant;
    visite[courant] = 1;

    double longueur = 0.0;

    for (int k = 1; k < N; ++k) {
        int meilleur = -1;
        double bestDist = DBL_MAX;

        for (int j = 1; j <= N; ++j) {
            if (!visite[j]) {
                double dij = matrice_distance(M, courant, j);
                if (dij < bestDist) {
                    bestDist = dij;
                    meilleur = j;
                }
            }
        }

        if (meilleur <= 0) {
            free(visite);
            fprintf(stderr, "NN(matrice): impossible de trouver un voisin\n");
            return -1.0;
        }

        tour->SECTION_TOUR[k] = meilleur;
        visite[meilleur] = 1;
        longueur += bestDist;
        courant = meilleur;
    }

    if (tour->FERMEE && N > 1) {
        longueur += matrice_distance(M, courant, depart);
    }

    tour->LONGUEUR = longueur;
    free(visite);
    return longueur;
}
