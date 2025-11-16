#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"
#include "tsp_force_brute.h"
#include "tsp_force_brute_matrice.h"
#include "tsp_nn.h"
#include "tsp_rw.h"

static void usage(const char* prog){
    fprintf(stderr,
        "Usage: %s -f <fichier.tsp> [-c] [-m bf|nn] [-M] [-F] [-o <fichier_sortie.txt>]\n"
        "  -h       : aide usage\n"
        "  -c       : tour canonique\n"
        "  -m bf    : force brute\n"
        "  -m nn    : plus proche voisin\n"
        "  -m rw    : marche aléatoire\n"
        "  -M       : avec -m bf ou -m nn : version demi-matrice\n"
        "  -F       : (avec -m bf) forcer si N>12\n"
        "  -o file  : met les resultats dans un fichier (append)\n",
        prog);
}

/* [1,2,3,...,N] sans espaces */
static void print_tour_compact_no_space(const TOUR_TSP* T){
    putchar('[');
    for (int i = 0; i < T->DIMENSION; ++i) {
        if (i) putchar(',');
        printf("%d", T->SECTION_TOUR[i]);
    }
    putchar(']');
}

/* Crée la tournée canonique 1..n fermée */
static TOUR_TSP tour_canonique(int n){
    TOUR_TSP T;
    T.DIMENSION = n;
    T.FERMEE = 1;
    T.LONGUEUR = -1.0;
    T.SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
    if (!T.SECTION_TOUR) { fprintf(stderr, "alloc tour\n"); exit(1); }
    for (int i = 0; i < n; ++i) T.SECTION_TOUR[i] = i + 1;
    return T;
}

/* Ligne normalisée unique */
static void print_line(const char* name, const char* algo, double secs, double length, const TOUR_TSP* T){
    printf("Tour %s %s %.6f %.0f ", name, algo, secs, length);
    print_tour_compact_no_space(T);
    putchar('\n');
}

int main(int argc, char** argv){
    int opt;
    FILE* output = stdout;
    const char* outfile = NULL;
    const char* filename = NULL;
    int want_canonical = 0;    /* -c */
    const char* method = NULL; /* -m bf|nn */
    int use_matrix = 0;        /* -M (pour bf ou nn) */
    int force_large = 0;       /* -F (pour bf) */

    while ((opt = getopt(argc, argv, "hf:cm:MFo:")) != -1) {
        switch (opt) {
            case 'h': usage(argv[0]); return 0;
            case 'f': filename = optarg; break;
            case 'c': want_canonical = 1; break;
            case 'm': method = optarg; break;
            case 'o': outfile = optarg; break;
            case 'M': use_matrix = 1; break;
            case 'F': force_large = 1; break;
            default: usage(argv[0]); return 1;
        }
    }
    if (!filename) { usage(argv[0]); return 1; }

    if (outfile){
        output = fopen(outfile,"a");
        if(!output){
            fprintf(stderr, "Erreur ouverture fichier sortie\n");
            return 1;
        }
        if(dup2(fileno(output), fileno(stdout)) == -1){
            fprintf(stderr, "Erreur dup2 redirection sortie\n");
            return 1;
        }
    }

    /* Si -m est fourni, il prime sur -c (on ignore -c) */
    if (method) want_canonical = 0;

    /* Lecture TSPLIB */
    TSPLIB_INSTANCE I;
    if (lire_tsplib(filename, &I) != 0) {
        /* ligne "safe" pour ne pas casser le parseur Python */
        printf("Tour unknown canonical 0.000000 0 []\n");
        return 2;
    }
    DistanceFn d = distance_pour(&I);
    if (!d) {
        printf("Tour %s canonical 0.000000 0 []\n", I.NAME[0] ? I.NAME : "unknown");
        liberer_instance(&I);
        return 3;
    }

    /* Tournée canonique (toujours dispo) */
    TOUR_TSP canon = tour_canonique(I.DIMENSION);

    /* Méthode canonique (si -c OU aucune méthode) */
    if (want_canonical || !method) {
        clock_t t0 = clock();
        double L = longueur_tour(&I, &canon, d);
        clock_t t1 = clock();
        double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

        print_line(I.NAME, "canonical", secs, L, &canon);

        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 0;
    }

    /* -m nn : plus proche voisin */
    if (strcmp(method, "nn") == 0) {
        TOUR_TSP T;
        T.DIMENSION = I.DIMENSION;
        T.FERMEE = 1;
        T.LONGUEUR = -1.0;
        T.SECTION_TOUR = (int*)malloc((size_t)I.DIMENSION * sizeof(int));
        if (!T.SECTION_TOUR) { fprintf(stderr, "alloc tour nn\n"); free(canon.SECTION_TOUR); liberer_instance(&I); return 4; }

        clock_t t0 = clock();
        double L;
        if (use_matrix) {
            /* Variante demi-matrice */
            MatriceTSP* M = creer_matrice_demie(&I, d);
            if (!M) {
                L = plus_proche_voisin(&I, d, &T, 1);
            } else {
                L = plus_proche_voisin_matrice(M, &T, 1);
                detruire_matrice_demie(M);
            }
        } else {
            L = plus_proche_voisin(&I, d, &T, 1);
        }
        clock_t t1 = clock();
        double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

        print_line(I.NAME, "nn", secs, L, &T);

        free(T.SECTION_TOUR);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 0;
    }

    if (strcmp(method, "rw") == 0) {
        TOUR_TSP T;
        T.DIMENSION = I.DIMENSION;
        T.FERMEE = 1;
        T.LONGUEUR = -1.0;
        T.SECTION_TOUR = (int*)malloc((size_t)I.DIMENSION * sizeof(int));
        if (!T.SECTION_TOUR) { fprintf(stderr, "alloc tour rw\n"); free(canon.SECTION_TOUR); liberer_instance(&I); return 4; }

        clock_t t0 = clock();
        double L;
        L = marche_aleatoire(&I,d,&T);
        clock_t t1 = clock();
        double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

        print_line(I.NAME, "rw", secs, L, &T);

        free(T.SECTION_TOUR);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 0;
    }

    /* -m bf : force brute  */
    if (strcmp(method, "bf") != 0) {
        fprintf(stderr, "Erreur: methode non supportee: %s (seuls 'canonical' via -c, 'bf' ou 'nn')\n", method);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 5;
    }

    /* Brute force demandée */
    if (I.DIMENSION > 12 && !force_large) {
        fprintf(stderr, "Erreur: DIMENSION=%d > 12. Utilisez -F pour forcer la brute force.\n", I.DIMENSION);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 6;
    }

    TOUR_TSP best = (TOUR_TSP){0}, worst = (TOUR_TSP){0};
    best.SECTION_TOUR = NULL; worst.SECTION_TOUR = NULL;

    double best_len;
    clock_t t0 = clock();
    if (use_matrix) {
        MatriceTSP* M = creer_matrice_demie(&I, d);
        if (!M) {
            best_len = force_brute(&I, d, &best, &worst);
        } else {
            best_len = force_brute_matrice(M, &best, &worst);
            detruire_matrice_demie(M);
        }
    } else {
        best_len = force_brute(&I, d, &best, &worst);
    }
    clock_t t1 = clock();
    double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

    if (best_len < 0.0) {
        fprintf(stderr, "Erreur: force brute a echoue.\n");
        if (best.SECTION_TOUR) free(best.SECTION_TOUR);
        if (worst.SECTION_TOUR) free(worst.SECTION_TOUR);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 7;
    }

    print_line(I.NAME, "bf", secs, best_len, &best);

    if (output != stdout) fclose(output);

    if (best.SECTION_TOUR) free(best.SECTION_TOUR);
    if (worst.SECTION_TOUR) free(worst.SECTION_TOUR);
    free(canon.SECTION_TOUR);
    liberer_instance(&I);
    return 0;
}
