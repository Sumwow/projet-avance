#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdbool.h>
#include "tsp_distance.h"

bool signal_recu = false;

void handler(int sig){
  (void) sig;
  signal_recu = true;
}

static void print_force_brute(double meilleur, double pire){
  printf("Meilleur longueur : %f\n",meilleur);
  printf("Pire longueur : %f\n",pire);
  printf("Entrez un caractère pour reprendre\n");
  
  getchar();
  printf("Reprise de force brute\n");
  signal_recu = false;
}

/* Locale : prochaine permutation lexicographique */
static int prochaine_permutation(int *villes, int nbVilles) {
    int i = nbVilles - 2;
    while (i >= 0 && villes[i] >= villes[i + 1]) i--;
    if (i < 0) return 0;

    int j = nbVilles - 1;
    while (villes[j] <= villes[i]) j--;

    int temp = villes[i]; villes[i] = villes[j]; villes[j] = temp;

    int debut = i + 1, fin = nbVilles - 1;
    while (debut < fin) {
        temp = villes[debut]; villes[debut] = villes[fin]; villes[fin] = temp;
        debut++; fin--;
    }
    return 1;
}

/* Force brute via DistanceFn (patch allocations inclus) */
double force_brute(const TSPLIB_INSTANCE* instance, DistanceFn calculDistance,
                   TOUR_TSP* meilleureTournee, TOUR_TSP* pireTournee) {

    int nbVilles = instance->DIMENSION;
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
    
    
    
    double longueurMin = longueur_tour(instance, &tourActuelle, calculDistance);
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
        
        signal(SIGINT,handler);
        
        double L = longueur_tour(instance, &tourActuelle, calculDistance);

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
        
        if(signal_recu){
          print_force_brute(meilleureTournee->LONGUEUR,pireTournee->LONGUEUR);
        }
    }

    free(ordreVilles);
    return longueurMin;
}

