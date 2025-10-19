#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"

/* Affiche l’utilisation du programme */
static void usage(const char* prog) {
    fprintf(stderr, "Usage: %s -f <fichier.tsp> [-c]\n", prog);
}

/* Crée la tournée canonique : 1, 2, 3, ..., N, retour au 1 */
static TOUR_TSP tour_canonique(int n) {
    TOUR_TSP T;
    T.DIMENSION = n;
    T.FERMEE = 1;
    T.LONGUEUR = -1.0;
    T.SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
    if (!T.SECTION_TOUR) {
        fprintf(stderr, "Erreur d’allocation mémoire pour la tournée.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++)
        T.SECTION_TOUR[i] = i + 1;
    return T;
}

int main(int argc, char** argv) {
    const char* file = NULL;
    int check = 0;

    /* Lecture des arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            file = argv[++i];
        } else if (strcmp(argv[i], "-c") == 0) {
            check = 1;
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (!file) {
        usage(argv[0]);
        return 1;
    }

    /* Lecture du fichier TSPLIB */
    TSPLIB_INSTANCE I;
    int rc = lire_tsplib(file, &I);
    if (rc != 0) {
        fprintf(stderr, "Erreur de lecture du fichier TSPLIB (%d)\n", rc);
        return 2;
    }

    /* Création de la tournée canonique */
    TOUR_TSP T = tour_canonique(I.DIMENSION);

    /* Sélection de la fonction de distance */
    DistanceFn d = distance_pour(&I);

    /* Calcul de la longueur de la tournée canonique */
    double L = longueur_tour(&I, &T, d);

    /* Calcul de la demi-matrice optimisée */
    MatriceTSP* M = creer_matrice_demie(&I, d);

    /* Affichage des informations de l’instance */
    if (check) {
        printf("NAME: %s\n", I.NAME);
        printf("TYPE: %s\n", I.TYPE);
        printf("DIMENSION: %d\n", I.DIMENSION);
        printf("EDGE_WEIGHT_TYPE: %s\n", I.EDGE_WEIGHT_TYPE);
    }

    /* Affiche toujours la longueur canonique */
    printf("CANONICAL_LENGTH=%.0f\n", L);

    /* Affiche toujours la demi-matrice triangulaire inférieure */
    printf("\n=== DEMI-MATRICE DES DISTANCES (INFÉRIEURE) ===\n\n");
    for (int i = 1; i <= I.DIMENSION; i++) {
        for (int j = 1; j <= i; j++) {
            double dist = matrice_distance(M, i, j);
            printf("%6.0f ", dist);
        }
        printf("\n");
    }
    printf("\n");

    /* Libération de la mémoire */
    free(T.SECTION_TOUR);
    detruire_matrice_demie(M);
    liberer_instance(&I);

    return 0;
}

