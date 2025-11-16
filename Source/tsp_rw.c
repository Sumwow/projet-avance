#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "tsp_rw.h"

double marche_aleatoire(const TSPLIB_INSTANCE* I, DistanceFn d, TOUR_TSP* tour){
    
    if (!I || I->DIMENSION <= 0 || !d || !tour) {
        fprintf(stderr, "NN(coord): arguments invalides\n");
        return -1.0;
    }
    const int N = I->DIMENSION;
    if (!tour->SECTION_TOUR || tour->DIMENSION != N) {
        fprintf(stderr, "NN(coord): tour non allouee ou dimension incoherente\n");
        return -1.0;
    }

    //Tableau de boolean pour savoir quels points on a déjà visisté, initialisé à false
    //Tableau à n+1 cases comme ça : point n = indice n et non pas indice n-1.
    //On se sert pas d'indice 0
    bool pointsVisite[N+1];
    memset(pointsVisite, false, sizeof(pointsVisite));

    //init random
    srand(time(NULL));

    //Tirage aléatoire du départ
    //Tire un nombre aléatoire parmis nos points, le tout + 1 pour éviter de tomber sur 0.
    int depart = (rand() % N) + 1; 
    pointsVisite[depart] = true;

    double longueur = 0.0;

    int courant = depart;
    tour->SECTION_TOUR[0] = courant;

    for (int i = 1; i < N; i++) {
        int next;
        do {
            next = (rand() % N) + 1;
        } while (pointsVisite[next]);

        longueur += d(I, courant, next);
        pointsVisite[next] = true;
        tour->SECTION_TOUR[i] = next;
        courant = next;
    }
    longueur += d(I, courant, depart);
    tour->FERMEE = 0;
    tour->LONGUEUR = longueur;
    
    return longueur;
}