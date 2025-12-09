#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "tsp_2opt.h"

static void inverser_segment(int* tour, int debut, int fin) {
    while (debut < fin) {
        int tempo = tour[debut];
        tour[debut] = tour[fin];
        tour[fin] = tempo;
        debut++;
        fin--;
    }
}

double two_opt(const TSPLIB_INSTANCE* instance, DistanceFn distance, TOUR_TSP* tour)
{
    if (instance == NULL || distance == NULL || tour == NULL || tour->SECTION_TOUR == NULL) {
        printf("Erreur : arguments invalides dans two_opt.\n");
        return -1.0;
    }
    int nbVilles = instance->DIMENSION;
    int* parcours = tour->SECTION_TOUR;
    bool ameliorationTrouvee = true;

    while (ameliorationTrouvee) {
        ameliorationTrouvee = false;
        for (int i = 0; i < nbVilles - 2; i++) {

            int villeA = parcours[i];
            int villeB = parcours[i + 1];
            for (int k = i + 2; k < nbVilles - 1; k++) {

                int villeC = parcours[k];
                int villeD = parcours[k + 1];
                
                int distAB = distance(instance, villeA, villeB);
                int distCD = distance(instance, villeC, villeD);

                int distAC = distance(instance, villeA, villeC);
                int distBD = distance(instance, villeB, villeD);

                if (distAC + distBD < distAB + distCD) {
                    inverser_segment(parcours, i + 1, k);
                    ameliorationTrouvee = true;
                }
            }
        }
    }

    double longueurTotale = longueur_tour(instance, tour, distance);

    tour->LONGUEUR = longueurTotale;
    return longueurTotale;
}

