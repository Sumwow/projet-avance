#include <stdio.h>
#include <stdlib.h>
#include "tsp_matrice.h"

/* Locales : helpers matrice + permutation */
static double lire_distance_matrice(const MatriceTSP* matrice, int i, int j) {
    if (i == j) return 0.0;
    if (i > j) { int t = i; i = j; j = t; }
    int ligne = i - 1;
    int colonne = j - i - 1;
    return matrice->data[ligne][colonne];
}

static double longueur_tour_matrice(const MatriceTSP* matrice, const TOUR_TSP* tour) {
    double L = 0.0;
    int n = tour->DIMENSION;
    for (int k = 0; k < n - 1; k++) {
        int a = tour->SECTION_TOUR[k];
        int b = tour->SECTION_TOUR[k + 1];
        L += lire_distance_matrice(matrice, a, b);
    }
    if (tour->FERMEE && n > 1) {
        int a = tour->SECTION_TOUR[n - 1];
        int b = tour->SECTION_TOUR[0];
        L += lire_distance_matrice(matrice, a, b);
    }
    return L;
}

static int prochaine_permutation(int *villes, int nbVilles) {
    int i = nbVilles - 2;
    while (i >= 0 && villes[i] >= villes[i + 1]) i--;
    if (i < 0) return 0;

    int j = nbVilles - 1;
    while (villes[j] <= villes[i]) j--;

    int tmp = villes[i]; villes[i] = villes[j]; villes[j] = tmp;

    int d = i + 1, f = nbVilles - 1;
    while (d < f) {
        tmp = villes[d]; villes[d] = villes[f]; villes[f] = tmp;
        d++; f--;
    }
    return 1;
}

/* Force brute via demi-matrice (patch allocations inclus) */
double force_brute_matrice(const MatriceTSP* matrice,
                           TOUR_TSP* meilleureTournee,
                           TOUR_TSP* pireTournee) {

    int nbVilles = matrice->dimension;
    int *ordreVilles = (int*)malloc((size_t)nbVilles * sizeof(int));
    if (!ordreVilles) {
        fprintf(stderr, "Erreur force brute (allocation)\n");
        return -1.0;
    }
    for (int i = 0; i < nbVilles; i++) ordreVilles[i] = i + 1;

    TOUR_TSP tourActuelle;
    tourActuelle.DIMENSION = nbVilles;
    tourActuelle.SECTION_TOUR = ordreVilles;
    tourActuelle.FERMEE = 1;

    double longueurMin = longueur_tour_matrice(matrice, &tourActuelle);
    double longueurMax = longueurMin;

    /* Patch minimal : allouer sorties si besoin */
    if (meilleureTournee->SECTION_TOUR == NULL || meilleureTournee->DIMENSION != nbVilles) {
        if (meilleureTournee->SECTION_TOUR) free(meilleureTournee->SECTION_TOUR);
        meilleureTournee->SECTION_TOUR = (int*)malloc((size_t)nbVilles * sizeof(int));
        if (!meilleureTournee->SECTION_TOUR) { free(ordreVilles); return -1.0; }
        meilleureTournee->DIMENSION = nbVilles;
        meilleureTournee->FERMEE = 1;
    }
    if (pireTournee->SECTION_TOUR == NULL || pireTournee->DIMENSION != nbVilles) {
        if (pireTournee->SECTION_TOUR) free(pireTournee->SECTION_TOUR);
        pireTournee->SECTION_TOUR = (int*)malloc((size_t)nbVilles * sizeof(int));
        if (!pireTournee->SECTION_TOUR) { free(ordreVilles); free(meilleureTournee->SECTION_TOUR); return -1.0; }
        pireTournee->DIMENSION = nbVilles;
        pireTournee->FERMEE = 1;
    }

    for (int i = 0; i < nbVilles; i++) {
        meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
        pireTournee->SECTION_TOUR[i]     = ordreVilles[i];
    }
    meilleureTournee->LONGUEUR = longueurMin;
    pireTournee->LONGUEUR      = longueurMax;

    while (prochaine_permutation(ordreVilles, nbVilles)) {
        double L = longueur_tour_matrice(matrice, &tourActuelle);

        if (L < longueurMin) {
            longueurMin = L;
            for (int i = 0; i < nbVilles; i++)
                meilleureTournee->SECTION_TOUR[i] = ordreVilles[i];
            meilleureTournee->LONGUEUR = L;
        }
        if (L > longueurMax) {
            longueurMax = L;
            for (int i = 0; i < nbVilles; i++)
                pireTournee->SECTION_TOUR[i] = ordreVilles[i];
            pireTournee->LONGUEUR = L;
        }
    }

    free(ordreVilles);
    return longueurMin;
}

