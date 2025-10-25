#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"
#include "tsp_force_brute.h"
#include "tsp_force_brute_matrice.h"

/* Affiche l’utilisation du programme */
static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s -f <fichier.tsp> [-c] [-b] [-M] [-F]\n"
        "  -f <fichier.tsp> : fichier TSPLIB (obligatoire)\n"
        "  -c               : afficher metadonnees (NAME/TYPE/DIM/EDGE_WEIGHT_TYPE)\n"
        "  -b               : lancer la force brute (meilleure & pire tournee)\n"
        "  -M               : (avec -b) utiliser la demi-matrice pour la force brute\n"
        "  -F               : forcer le bruteforce meme si DIMENSION>12 (tres lent)\n",
        prog);
}

/* Crée la tournée canonique : 1..N (fermée) */
static TOUR_TSP tour_canonique(int n) {
    TOUR_TSP T;
    T.DIMENSION = n;
    T.FERMEE = 1;
    T.LONGUEUR = -1.0;
    T.SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
    if (!T.SECTION_TOUR) {
        fprintf(stderr, "Erreur d’allocation memoire pour la tournee.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++) T.SECTION_TOUR[i] = i + 1;
    return T;
}

static void print_tour_inline(const TOUR_TSP* T) {
    for (int i = 0; i < T->DIMENSION; ++i) {
        printf("%d", T->SECTION_TOUR[i]);
        if (i + 1 < T->DIMENSION) printf(" ");
    }
    printf("\n");
}

int main(int argc, char** argv) {
    const char* file = NULL;
    FILE* output = stdout;
    char* outfile = NULL;  /* -o */
    int check = 0;         /* -c */
    int do_bruteforce = 0; /* -b */
    int use_matrix_bf = 0; /* -M */
    int force_large = 0;   /* -F */
    
    //Si -h alors le programme affiche juste le usage.
    if(strcmp(argv[1],"-h") == 0){
      usage(argv[0]); return 1;
    }
    
    /* Args */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            file = argv[++i];
        } else if(strcmp(argv[i], "-o") == 0 && i + 1 < argc){
            outfile = argv[i+1];
        } else if (strcmp(argv[i], "-c") == 0) {
            check = 1;
        } else if (strcmp(argv[i], "-b") == 0) {
            do_bruteforce = 1;
        } else if (strcmp(argv[i], "-M") == 0) {
            use_matrix_bf = 1;
        } else if (strcmp(argv[i], "-F") == 0) {
            force_large = 1;
        } else {
            usage(argv[0]); return 1;
        }
    }
    if (!file) { usage(argv[0]); return 1; }
    
    if(outfile){
      output = fopen(outfile,"a");
      if(!output){
        printf("Erreur ouverture fichier sortie");
        return 1;
      }
      if(dup2(fileno(output), fileno(stdout)) == 1){
        printf("Erreur dup2 redirection sortie");
        return 1;
      }
    }
    
    
    /* Lecture TSPLIB */
    TSPLIB_INSTANCE I;
    int rc = lire_tsplib(file, &I);
    if (rc != 0) {
        fprintf(stderr, "Erreur de lecture TSPLIB (%d)\n", rc);
        return 2;
    }

    /* Canonique + distance fn + demi-matrice */
    TOUR_TSP T = tour_canonique(I.DIMENSION);
    DistanceFn d = distance_pour(&I);
    double L = longueur_tour(&I, &T, d);

    MatriceTSP* M = creer_matrice_demie(&I, d);
    if (!M) {
        fprintf(stderr, "Erreur: creation de la demi-matrice echouee.\n");
        free(T.SECTION_TOUR);
        liberer_instance(&I);
        return 3;
    }

    if (check) {
        printf("NAME: %s\n", I.NAME);
        printf("TYPE: %s\n", I.TYPE);
        printf("DIMENSION: %d\n", I.DIMENSION);
        printf("EDGE_WEIGHT_TYPE: %s\n", I.EDGE_WEIGHT_TYPE);
    }

    /* Longueur canonique */
    printf("CANONICAL_LENGTH=%.0f\n", L);

    /* Demi-matrice toujours affichée */
    printf("\n=== DEMI-MATRICE DES DISTANCES (INFÉRIEURE) ===\n\n");
    for (int i = 1; i <= I.DIMENSION; i++) {
        for (int j = 1; j <= i; j++) {
            double dist = matrice_distance(M, i, j);
            printf("%6.0f ", dist);
        }
        printf("\n");
    }
    printf("\n");

    /* Force brute (optionnel) */
    if (do_bruteforce) {
        if (I.DIMENSION > 12 && !force_large) {
            fprintf(stderr, "Attention: DIMENSION=%d trop grande pour brute force (max 12). Utilisez -F pour forcer.\n",
                    I.DIMENSION);
        } else {
            TOUR_TSP best = {0}, worst = {0};
            double best_len = -1.0, worst_len = -1.0;

            if (use_matrix_bf) {
                best.SECTION_TOUR = NULL; worst.SECTION_TOUR = NULL;
                best_len = force_brute_matrice(M, &best, &worst);
                worst_len = worst.LONGUEUR;

                printf("=== FORCE BRUTE (Matrice) ===\n");
                printf("Meilleure longueur  = %.0f\n", best_len);
                printf("Meilleure tournee    : ");
                print_tour_inline(&best);
                printf("Pire longueur       = %.0f\n", worst_len);
                printf("Pire tournee         : ");
                print_tour_inline(&worst);

                if (best.SECTION_TOUR)  free(best.SECTION_TOUR);
                if (worst.SECTION_TOUR) free(worst.SECTION_TOUR);
            } else {
                best.SECTION_TOUR = NULL; worst.SECTION_TOUR = NULL;
                best_len = force_brute(&I, d, &best, &worst);
                worst_len = worst.LONGUEUR;

                printf("=== FORCE BRUTE (Distance) ===\n");
                printf("Meilleure longueur  = %.0f\n", best_len);
                printf("Meilleure tournee    : ");
                print_tour_inline(&best);
                printf("Pire longueur       = %.0f\n", worst_len);
                printf("Pire tournee         : ");
                print_tour_inline(&worst);

                if (best.SECTION_TOUR)  free(best.SECTION_TOUR);
                if (worst.SECTION_TOUR) free(worst.SECTION_TOUR);
            }
            printf("\n");
        }
    }
  
    if(output != stdout) fclose(output);
  
    /* Free */
    free(T.SECTION_TOUR);
    detruire_matrice_demie(M);
    liberer_instance(&I);
    return 0;
}

