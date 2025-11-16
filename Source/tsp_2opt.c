#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "tsp_2opt.h"

/*On inverse positions du point debut et fin */
static void inverser_segment(int* tour, int debut, int fin) {
    while (debut < fin) {
        int tempo = tour[debut];
        tour[debut] = tour[fin];
        tour[fin] = tempo;
        debut++;
        fin--;
    }
}

double two_opt(const TSPLIB_INSTANCE* instance,
               DistanceFn distance,
               TOUR_TSP* tour)
{
    if (instance == NULL || distance == NULL || tour == NULL || tour->SECTION_TOUR == NULL) {
        printf("Erreur : arguments invalides dans two_opt.\n");
        return -1.0;
    }

    int nbVilles = instance->DIMENSION;
    int* parcours = tour->SECTION_TOUR;

    bool ameliorationTrouvee = true;

    /*On boucle tant qu'une amélioration du temps est trouvée*/
    while (ameliorationTrouvee) {

        ameliorationTrouvee = false;

        /*i correspond au premier segment*/
        for (int i = 0; i < nbVilles - 2; i++) {

            int villeA = parcours[i];
            int villeB = parcours[i + 1];

            /*k correspond au deuxième segment*/
            for (int k = i + 2; k < nbVilles - 1; k++) {

                int villeC = parcours[k];
                int villeD = parcours[k + 1];

                /*On récupère la distance avant 2-opt soit: (A-B) + (C-D)*/
                int distAB = distance(instance, villeA, villeB);
                int distCD = distance(instance, villeC, villeD);

                /*On récupère la distance après 2-opt soit: (A-C) + (B-D) */
                int distAC = distance(instance, villeA, villeC);
                int distBD = distance(instance, villeB, villeD);

                /*Si la modification améliore la tournée, alors on inverse*/
                if (distAC + distBD < distAB + distCD) {
                    inverser_segment(parcours, i + 1, k);
                    ameliorationTrouvee = true;
                }
            }
        }
    }

    /*On calcul la nouvelle longueur totale de la tournée*/
    double longueurTotale = 0.0;

    for (int i = 0; i < nbVilles - 1; i++) {
        longueurTotale += distance(instance, parcours[i], parcours[i + 1]);
    }

    /*Si la tournée est fermée, alors on rajoute le segment de retour*/
    if (tour->FERMEE) {
        longueurTotale += distance(instance, parcours[nbVilles - 1], parcours[0]);
    }

    tour->LONGUEUR = longueurTotale;
    return longueurTotale;
}