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
#include "tsp_2opt.h"
#include "tsp_ga.h"
#include "tsp_dpx.h"


static void usage(const char* prog){
    fprintf(stderr,
        "Usage: %s -f <fichier.tsp> [-c] [-m bf|nn|rw|2optnn|2optrw|ga|gadpx] [-M] [-F] [-o <fichier_sortie.txt>]\n"
        "  -h       : aide usage\n"
        "  -c       : tour canonique\n"
        "  -m bf    : force brute\n"
        "  -m nn    : plus proche voisin\n"
        "  -m rw    : marche aleatoire\n"
        "  -m 2optnn: 2-opt sur plus proche voisin\n"
        "  -m 2optrw: 2-opt sur marche aleatoire\n"
        "  -m ga <Ngens> <taux_mut> <Nindiv> : algorithme genetique\n"
        "  -m gadpx <Ngens> <taux_mut> <Nindiv> : GA + DPX\n"
        "  -M       : avec -m bf ou -m nn : version demi-matrice\n"
        "  -F       : (avec -m bf) forcer si N>12\n"
        "  -o file  : append des resultats dans un fichier\n",
        prog);
}

static void print_tour_compact_no_space(const TOUR_TSP* T){
    putchar('[');
    for (int i = 0; i < T->DIMENSION; ++i) {
        if (i) putchar(',');
        printf("%d", T->SECTION_TOUR[i]);
    }
    putchar(']');
}


static TOUR_TSP tour_canonique(int n){
    TOUR_TSP T;
    T.DIMENSION = n;
    T.FERMEE    = 1;
    T.LONGUEUR  = -1.0;
    T.SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
    if (!T.SECTION_TOUR) { fprintf(stderr, "alloc tour\n"); exit(1); }
    for (int i = 0; i < n; ++i) T.SECTION_TOUR[i] = i + 1;
    return T;
}

static void print_line(const char* name, const char* algo,
                       double secs, double length, const TOUR_TSP* T){
    printf("Tour %s %s %.6f %.0f ", name, algo, secs, length);
    print_tour_compact_no_space(T);
    putchar('\n');
}

typedef enum {
    ALGO_CANONICAL,
    ALGO_BF,
    ALGO_NN,
    ALGO_RW,
    ALGO_2OPTNN,
    ALGO_2OPTRW,
    ALGO_GA,
    ALGO_GADPX,
    ALGO_UNKNOWN
} AlgoType;

static AlgoType parse_method(const char* m){
    if (!m)              return ALGO_CANONICAL;
    if (!strcmp(m,"bf")) return ALGO_BF;
    if (!strcmp(m,"nn")) return ALGO_NN;
    if (!strcmp(m,"rw")) return ALGO_RW;
    if (!strcmp(m,"2optnn")) return ALGO_2OPTNN;
    if (!strcmp(m,"2optrw")) return ALGO_2OPTRW;
    if (!strcmp(m,"ga"))     return ALGO_GA;
    if (!strcmp(m,"gadpx"))  return ALGO_GADPX;
    return ALGO_UNKNOWN;
}

/* ============================== main =============================== */

int main(int argc, char** argv){
    int opt;
    FILE* output = stdout;
    const char* outfile   = NULL;
    const char* filename  = NULL;
    int want_canonical    = 0;
    const char* methodStr = NULL;
    int use_matrix        = 0;
    int force_large       = 0;

    while ((opt = getopt(argc, argv, "hf:cm:MFo:")) != -1) {
        switch (opt) {
            case 'h': usage(argv[0]); return 0;
            case 'f': filename       = optarg; break;
            case 'c': want_canonical = 1;      break;
            case 'm': methodStr      = optarg; break;
            case 'o': outfile        = optarg; break;
            case 'M': use_matrix     = 1;      break;
            case 'F': force_large    = 1;      break;
            default : usage(argv[0]); return 1;
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

    AlgoType algo;
    if (methodStr){
        algo = parse_method(methodStr);
    } else if (want_canonical){
        algo = ALGO_CANONICAL;
    } else {
        algo = ALGO_CANONICAL;   /* défaut */
    }

    if (algo == ALGO_UNKNOWN){
        fprintf(stderr, "Erreur: methode non supportee: %s\n", methodStr);
        return 5;
    }

    TSPLIB_INSTANCE I;
    if (lire_tsplib(filename, &I) != 0) {
        printf("Tour unknown canonical 0.000000 0 []\n");
        return 2;
    }
    DistanceFn d = distance_pour(&I);
    if (!d) {
        printf("Tour %s canonical 0.000000 0 []\n", I.NAME[0] ? I.NAME : "unknown");
        liberer_instance(&I);
        return 3;
    }

    TOUR_TSP canon = tour_canonique(I.DIMENSION);

    TOUR_TSP work = {0};
    int need_work =
        (algo == ALGO_NN || algo == ALGO_RW ||
         algo == ALGO_2OPTNN || algo == ALGO_2OPTRW);

    if (need_work){
        work.DIMENSION = I.DIMENSION;
        work.FERMEE    = 1;
        work.LONGUEUR  = -1.0;
        work.SECTION_TOUR = (int*)malloc((size_t)I.DIMENSION * sizeof(int));
        if (!work.SECTION_TOUR){
            fprintf(stderr, "alloc tour travail\n");
            liberer_instance(&I);
            free(canon.SECTION_TOUR);
            return 4;
        }
    }

    TOUR_TSP best = {0}, worst = {0};
    TOUR_TSP best_ga = {0};

    double L = -1.0;
    double secs = 0.0;
    int status = 0;

    clock_t t0 = clock();

    switch (algo){

        case ALGO_CANONICAL: {
            L = longueur_tour(&I, &canon, d);
            break;
        }

        case ALGO_NN: {
            if (use_matrix){
                MatriceTSP* M = creer_matrice_demie(&I, d);
                if (!M) L = plus_proche_voisin(&I, d, &work, 1);
                else {
                    L = plus_proche_voisin_matrice(M, &work, 1);
                    detruire_matrice_demie(M);
                }
            } else {
                L = plus_proche_voisin(&I, d, &work, 1);
            }
            break;
        }

        case ALGO_RW: {
            L = marche_aleatoire(&I, d, &work);
            break;
        }

        case ALGO_2OPTRW: {
            marche_aleatoire(&I, d, &work);
            work.FERMEE = 1;
            L = two_opt(&I, d, &work);
            break;
        }

        case ALGO_2OPTNN: {
            plus_proche_voisin(&I, d, &work, 1);
            L = two_opt(&I, d, &work);
            break;
        }

        case ALGO_GA: {
            int pop_size    = 30;
            int generations = 1000;
            double mutation_rate = 0.10;

            if (optind + 2 < argc) {
                int tmp_gen = atoi(argv[optind]);
                double tmp_mut = atof(argv[optind+1]);
                int tmp_pop = atoi(argv[optind+2]);
                if (tmp_pop > 0) pop_size = tmp_pop;
                if (tmp_gen > 0) generations = tmp_gen;
                if (tmp_mut >= 0.0 && tmp_mut <= 1.0) mutation_rate = tmp_mut;
            }

            best_ga.SECTION_TOUR = NULL;
            L = tsp_ga_light(&I, d, pop_size, generations, mutation_rate, &best_ga);
            if (L < 0.0){
                fprintf(stderr, "Erreur: GA a echoue.\n");
                status = 9;
            }
            break;
        }

        case ALGO_GADPX: {
            int pop_size    = 30;
            int generations = 1000;
            double mutation_rate = 0.10;

            if (optind + 2 < argc) {
                int tmp_gen = atoi(argv[optind]);
                double tmp_mut = atof(argv[optind+1]);
                int tmp_pop = atoi(argv[optind+2]);
                if (tmp_pop > 0) pop_size = tmp_pop;
                if (tmp_gen > 0) generations = tmp_gen;
                if (tmp_mut >= 0.0 && tmp_mut <= 1.0) mutation_rate = tmp_mut;
            }

            best_ga.SECTION_TOUR = NULL;
            L = tsp_ga_dpx(&I, d, pop_size, generations, mutation_rate, &best_ga);
            if (L < 0.0){
                fprintf(stderr, "Erreur: GA+DPX a echoue.\n");
                status = 9;
            }
            break;
        }

        case ALGO_BF: {
            if (I.DIMENSION > 12 && !force_large){
                fprintf(stderr, "Erreur: DIMENSION=%d > 12. Utilisez -F pour forcer la brute force.\n",
                        I.DIMENSION);
                status = 6;
                break;
            }
            best.SECTION_TOUR = NULL;
            worst.SECTION_TOUR = NULL;

            if (use_matrix){
                MatriceTSP* M = creer_matrice_demie(&I, d);
                if (!M) {
                    L = force_brute(&I, d, &best, &worst);
                } else {
                    L = force_brute_matrice(M, &best, &worst);
                    detruire_matrice_demie(M);
                }
            } else {
                L = force_brute(&I, d, &best, &worst);
            }
            if (L < 0.0){
                fprintf(stderr, "Erreur: force brute a echoue.\n");
                status = 7;
            }
            break;
        }

        default:
            status = 5;
            break;
    }

    clock_t t1 = clock();
    secs = (double)(t1 - t0) / CLOCKS_PER_SEC;

    if (status != 0) goto fin;

    switch (algo){
        case ALGO_CANONICAL:
            print_line(I.NAME, "canonical", secs, L, &canon);
            break;
        case ALGO_NN:
            print_line(I.NAME, "nn", secs, L, &work);
            break;
        case ALGO_RW:
            print_line(I.NAME, "rw", secs, L, &work);
            break;
        case ALGO_2OPTRW:
            print_line(I.NAME, "2optrw", secs, L, &work);
            break;
        case ALGO_2OPTNN:
            print_line(I.NAME, "2optnn", secs, L, &work);
            break;
        case ALGO_GA:
            print_line(I.NAME, "ga", secs, L, &best_ga);
            break;
        case ALGO_GADPX:
            print_line(I.NAME, "gadpx", secs, L, &best_ga);
            break;
        case ALGO_BF:
            print_line(I.NAME, "bf", secs, L, &best);
            break;
        default:
            break;
    }

fin:
    if (output != stdout) fclose(output);

    if (work.SECTION_TOUR)     free(work.SECTION_TOUR);
    if (best.SECTION_TOUR)     free(best.SECTION_TOUR);
    if (worst.SECTION_TOUR)    free(worst.SECTION_TOUR);
    if (best_ga.SECTION_TOUR)  free(best_ga.SECTION_TOUR);
    free(canon.SECTION_TOUR);
    liberer_instance(&I);

    return status;
}

