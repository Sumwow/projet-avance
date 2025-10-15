#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"

static void usage(const char* prog){
    fprintf(stderr, "Usage: %s -f <fichier.tsp> [-c]\n", prog);
}

static TOUR_TSP tour_canonique(int n) {
    TOUR_TSP T; T.DIMENSION = n; T.FERMEE = 1; T.LONGUEUR = -1.0;
    T.SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
    for (int i=0;i<n;i++) T.SECTION_TOUR[i] = i+1;
    return T;
}

int main(int argc, char** argv){
    const char* file = NULL; int check = 0;
    for (int i=1;i<argc;i++){
        if (strcmp(argv[i], "-f")==0 && i+1<argc){ file = argv[++i]; }
        else if (strcmp(argv[i], "-c")==0){ check = 1; }
        else { usage(argv[0]); return 1; }
    }
    if (!file){ usage(argv[0]); return 1; }

    TSPLIB_INSTANCE I; int rc = lire_tsplib(file, &I);
    if (rc!=0){ fprintf(stderr, "Erreur lecture TSPLIB (%d)\n", rc); return 2; }

    TOUR_TSP T = tour_canonique(I.DIMENSION);
    DistanceFn d = distance_pour(&I);
    double L = longueur_tour(&I, &T, d);

    if (check){
        printf("NAME: %s\n", I.NAME);
        printf("TYPE: %s\n", I.TYPE);
        printf("DIMENSION: %d\n", I.DIMENSION);
        printf("EDGE_WEIGHT_TYPE: %s\n", I.EDGE_WEIGHT_TYPE);
    }
    /* Ligne machine-lisible utile aux tests */
    printf("CANONICAL_LENGTH=%.0f\n", L);

    free(T.SECTION_TOUR);
    liberer_instance(&I);
    return 0;
}
