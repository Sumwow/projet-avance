#include <stdio.h>
#include <stdlib.h>
#include "tsp_matrice.h"

double lire_distance_matrice(const MatriceTSP* matrice, int i, int j) {
    if (i == j) return 0.0;

    // On s'assure que i < j pour accéder correctement à la matrice
    if (i > j) {
        int temp = i;
        i = j;
        j = temp;
    }

    int ligne = i - 1;           // ligne correspondant à la ville i
    int colonne = j - i - 1;     // colonne correspondant à la ville j
    return matrice->data[ligne][colonne];
}


double longueur_tour_matrice(const MatriceTSP* matrice, const TOUR_TSP* tour) {
    double longueur = 0.0;
    int nbVilles = tour->DIMENSION;

    // On additionne les distances entre chaque paire consécutive
    for (int i = 0; i < nbVilles - 1; i++) {
        int ville1 = tour->SECTION_TOUR[i];
        int ville2 = tour->SECTION_TOUR[i + 1];
        longueur += lire_distance_matrice(matrice, ville1, ville2);
    }

    // Si la tournée est fermée, on ajoute la distance retour
    if (tour->FERMEE) {
        int derniereVille = tour->SECTION_TOUR[nbVilles - 1];
        int premiereVille = tour->SECTION_TOUR[0];
        longueur += lire_distance_matrice(matrice, derniereVille, premiereVille);
    }

    return longueur;
}

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

double force_brute_matrice(const MatriceTSP* matrice, TOUR_TSP* meilleureTournee, TOUR_TSP* pireTournee) {

    int nbVilles = matrice->dimension;

    // On stocke l'ordre des villes dans un tableau
    int *ordreVilles = malloc(nbVilles * sizeof(int));
    if (ordreVilles == NULL) {
        printf("Erreur force brute (allocation)\n");
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

    // On calcule la longueur de la tournée de départ
    double longueurMin = longueur_tour_matrice(matrice, &tourActuelle);
    double longueurMax = longueurMin;

    // On la considère comme la meilleure et la pire tournée
    for (int i = 0; i < nbVilles; i++) {
        meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
        pireTournee->SECTION_TOUR[i] = ordreVilles[i];
    }
    meilleureTournee->LONGUEUR = longueurMin;
    pireTournee->LONGUEUR = longueurMax;

    // On parcourt toutes les permutations suivantes
    while (prochaine_permutation(ordreVilles, nbVilles)) {

        double longueurActuelle = longueur_tour_matrice(matrice, &tourActuelle);

        // Si la tournée actuelle est plus courte → nouvelle meilleure tournée
        if (longueurActuelle < longueurMin) {
            longueurMin = longueurActuelle;
            for (int i = 0; i < nbVilles; i++)
                meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
            meilleureTournee->LONGUEUR = longueurActuelle;
        }

        // Si la tournée actuelle est plus longue → nouvelle pire tournée
        if (longueurActuelle > longueurMax) {
            longueurMax = longueurActuelle;
            for (int i = 0; i < nbVilles; i++)
                pireTournee->SECTION_TOUR[i] = ordreVilles[i];
            pireTournee->LONGUEUR = longueurActuelle;
        }
    }

    // On libère la mémoire allouée pour le tableau temporaire
    free(ordreVilles);

    // On renvoie la longueur de la meilleure tournée
    return longueurMin;
}