#include <stdio.h>
#include <stdlib.h>
#include "tsp_distance.h"

/* On permute dans l'ordre lexicographique.
   Ca renvoie 1 si y'a une permutation et 0 dans le cas contraireµ*/
int prochaine_permutation(int *villes, int nbVilles) {
    int i = nbVilles - 2;

    // On cherche le premier élément qui est plus petit que le suivant
    while (i >= 0 && villes[i] >= villes[i + 1]) {
        i--;
    }
    if (i < 0) {
        return 0;
    }

    // On cherche le plus petit élément > à villes[i]
    int j = nbVilles - 1;
    while (villes[j] <= villes[i]) {
        j--;
    }

    // On échange les deux éléments trouvés
    int temp = villes[i];
    villes[i] = villes[j];
    villes[j] = temp;

    // On renverse la partie du tableau après i
    int debut = i + 1;
    int fin = nbVilles - 1;
    while (debut < fin) {
        temp = villes[debut];
        villes[debut] = villes[fin];
        villes[fin] = temp;
        debut++;
        fin--;
    }

    return 1;
}

double force_brute(const TSPLIB_INSTANCE* instance, DistanceFn calculDistance, TOUR_TSP* meilleureTournee, TOUR_TSP* pireTournee) {

    int nbVilles = instance->DIMENSION;

    // On stocke l'ordre des villes dans un tableau
    int *ordreVilles = malloc(nbVilles * sizeof(int));
    if (ordreVilles == NULL) {
        printf("Erreur force brute\n");
        return -1;
    }

    // On crée la tournée de départ : [1, 2, 3, ..., N]
    for (int i = 0; i < nbVilles; i++) {
        ordreVilles[i] = i + 1;
    }

    TOUR_TSP tourActuelle;
    tourActuelle.DIMENSION = nbVilles;
    tourActuelle.SECTION_TOUR = ordreVilles;
    tourActuelle.FERMEE = 1; // on revient à la ville de départ

    // On calcul la longueur de la tournée de départ
    double longueurMin = longueur_tour(instance, &tourActuelle, calculDistance);
    double longueurMax = longueurMin;

    // On la considère comme la pire et la meilleure tournée
    for (int i = 0; i < nbVilles; i++) {
        meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
        pireTournee->SECTION_TOUR[i] = ordreVilles[i];
    }
    meilleureTournee->LONGUEUR = longueurMin;
    pireTournee->LONGUEUR = longueurMax;

    // On parcourt toutes les permutations suivantes
    while (prochaine_permutation(ordreVilles, nbVilles)) {
        double longueurActuelle = longueur_tour(instance, &tourActuelle, calculDistance);

        // On considère la tournée actuelle comme étant la meilleure si elle est plus courte
        if (longueurActuelle < longueurMin) {
            longueurMin = longueurActuelle;
            for (int i = 0; i < nbVilles; i++)
                meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
            meilleureTournee->LONGUEUR = longueurActuelle;
        }

        // On considère la tournée actuelle comme étant la pire si elle est plus longue
        if (longueurActuelle > longueurMax) {
            longueurMax = longueurActuelle;
            for (int i = 0; i < nbVilles; i++)
                pireTournee->SECTION_TOUR[i] = ordreVilles[i];
            pireTournee->LONGUEUR = longueurActuelle;
        }
    }

    free(ordreVilles);
    return longueurMin;
}