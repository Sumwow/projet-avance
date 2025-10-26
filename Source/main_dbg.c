#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>   /* getopt */
#include <time.h>

#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"
#include "tsp_force_brute.h"
#include "tsp_force_brute_matrice.h"

static void usage(const char* prog){
    fprintf(stderr,
        "Usage: %s -f <fichier.tsp> [-c] [-m bf] [-M] [-F] [-X]\n"
        "  -h       : aide usage\n"
        "  -c       : canonique\n"
        "  -m bf    : force brute\n"
        "  -M       : (bf) version matrice\n"
        "  -F       : (bf) forcer si N>12\n"
        "  -X       : affichages detailles (stderr)\n", prog);
}

/* [1,2,3,...,N] sans espaces — pour stdout */
static void print_tour_compact_no_space_stdout(const TOUR_TSP* T){
    putchar('[');
    for (int i = 0; i < T->DIMENSION; ++i) {
        if (i) putchar(',');
        printf("%d", T->SECTION_TOUR[i]);
    }
    putchar(']');
}

/* [1,2,3,...,N] sans espaces — pour stderr */
static void fprint_tour_compact_no_space_stderr(const TOUR_TSP* T){
    fprintf(stderr, "[");
    for (int i = 0; i < T->DIMENSION; ++i) {
        if (i) fprintf(stderr, ",");
        fprintf(stderr, "%d", T->SECTION_TOUR[i]);
    }
    fprintf(stderr, "]\n");
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

/* Ligne normalisée unique sur stdout */
static void print_line_stdout(const char* name, const char* algo, double secs, double length, const TOUR_TSP* T){
    printf("Tour %s %s %.6f %.0f ", name, algo, secs, length);
    print_tour_compact_no_space_stdout(T);
    putchar('\n');
}

int main(int argc, char** argv){
    int opt;
    const char* filename = NULL;
    int want_canonical = 0;    /* -c */
    const char* method = NULL; /* -m bf */
    int use_matrix = 0;        /* -M */
    int force_large = 0;       /* -F */
    int verbose = 0;           /* -X : affiche sur stderr */

    while ((opt = getopt(argc, argv, "hf:cm:MFX")) != -1) {
        switch (opt) {
            case 'h': usage(argv[0]); return 0;
            case 'f': filename = optarg; break;
            case 'c': want_canonical = 1; break;
            case 'm': method = optarg; break;
            case 'M': use_matrix = 1; break;
            case 'F': force_large = 1; break;
            case 'X': verbose = 1; break;
            default: usage(argv[0]); return 1;
        }
    }
    if (!filename) { usage(argv[0]); return 1; }
    if (method) want_canonical = 0; /* -m prime sur -c */

    /* Lecture instance */
    TSPLIB_INSTANCE I;
    if (lire_tsplib(filename, &I) != 0) {
        printf("Tour unknown canonical 0.000000 0 []\n"); /* safe pour le parser */
        return 2;
    }
    DistanceFn d = distance_pour(&I);
    if (!d) {
        printf("Tour %s canonical 0.000000 0 []\n", I.NAME[0] ? I.NAME : "unknown");
        liberer_instance(&I);
        return 3;
    }

    /* Canonique dispo partout */
    TOUR_TSP canon = tour_canonique(I.DIMENSION);

    /* Canonique (ou pas de -m) */
    if (want_canonical || !method) {
        clock_t t0 = clock();
        double L = longueur_tour(&I, &canon, d);
        clock_t t1 = clock();
        double secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

        print_line_stdout(I.NAME, "canonical", secs, L, &canon);

        if (verbose) {
            fprintf(stderr, "[DBG] Instance=%s, DIM=%d, EDGE_WEIGHT_TYPE=%s\n",
                    I.NAME, I.DIMENSION, I.EDGE_WEIGHT_TYPE);
            fprintf(stderr, "[DBG] Canonical length=%.0f, time=%.6f s\n", L, secs);

            MatriceTSP* Mdbg = creer_matrice_demie(&I, d);
            if (Mdbg) {
                fprintf(stderr, "\n=== DEMI-MATRICE (INFÉRIEURE) ===\n\n");
                for (int i = 1; i <= I.DIMENSION; ++i) {
                    for (int j = 1; j <= i; ++j) {
                        fprintf(stderr, "%6.0f ", matrice_distance(Mdbg, i, j));
                    }
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "\n");
                detruire_matrice_demie(Mdbg);
            }
        }

        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 0;
    }

    /* Force brute uniquement si -m bf */
    if (strcmp(method, "bf") != 0) {
        fprintf(stderr, "Erreur: methode non supportee: %s (seuls -c ou -m bf)\n", method);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 4;
    }
    if (I.DIMENSION > 12 && !force_large) {
        fprintf(stderr, "Erreur: DIMENSION=%d > 12. Utilisez -F pour forcer.\n", I.DIMENSION);
        free(canon.SECTION_TOUR);
        liberer_instance(&I);
        return 5;
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
        return 6;
    }

    /* stdout — ligne normalisée */
    print_line_stdout(I.NAME, "bf", secs, best_len, &best);

    /* stderr — détails avec meilleure longueur */
    if (verbose) {
        fprintf(stderr, "[DBG] FORCE BRUTE RESULTATS\n");
        fprintf(stderr, "  -> Meilleure longueur : %.0f\n", best_len);
        fprintf(stderr, "  -> Temps CPU          : %.6f s\n", secs);
        fprintf(stderr, "  -> Meilleure tournée  : ");
        fprint_tour_compact_no_space_stderr(&best);

        if (worst.SECTION_TOUR && worst.DIMENSION > 0) {
            fprintf(stderr, "  -> Pire longueur      : %.0f\n", worst.LONGUEUR);
            fprintf(stderr, "  -> Pire tournée       : ");
            fprint_tour_compact_no_space_stderr(&worst);
        }

        MatriceTSP* Mdbg = creer_matrice_demie(&I, d);
        if (Mdbg) {
            fprintf(stderr, "\n=== DEMI-MATRICE (INFÉRIEURE) ===\n\n");
            for (int i = 1; i <= I.DIMENSION; ++i) {
                for (int j = 1; j <= i; ++j) {
                    fprintf(stderr, "%6.0f ", matrice_distance(Mdbg, i, j));
                }
                fprintf(stderr, "\n");
            }
            fprintf(stderr, "\n");
            detruire_matrice_demie(Mdbg);
        }
    }

    if (best.SECTION_TOUR) free(best.SECTION_TOUR);
    if (worst.SECTION_TOUR) free(worst.SECTION_TOUR);
    free(canon.SECTION_TOUR);
    liberer_instance(&I);
    return 0;
}
