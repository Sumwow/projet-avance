#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "tsp_types.h"
#include "tsp_io.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"
//#include "tsp_force_brute.h"
//#include "tsp_force_brute_matrice.h"
#include "tsp_force_brute_generique.h"
#include "tsp_nn.h"
#include "tsp_rw.h"
#include "tsp_2opt.h"
#include "tsp_ga.h"
#include "tsp_dpx.h"

static const TSPLIB_INSTANCE* gI_coord = NULL;
static DistanceFn gD_coord = NULL;


static void brute_set_instance(const TSPLIB_INSTANCE* I, DistanceFn d) {
    gI_coord = I;
    gD_coord = d;
}

static void* cout_tsp_coord(void* unused, int* permutation_villes) {
    (void)unused;

    if (!gI_coord || !gD_coord || !permutation_villes)
        return NULL;

    int n = gI_coord->DIMENSION;
    double longueur_totale = 0.0;

    for (int i = 0; i < n - 1; ++i) {
        int a = permutation_villes[i];
        int b = permutation_villes[i + 1];
        longueur_totale += gD_coord(gI_coord, a, b);
    }

    if (n > 1) {
        int last  = permutation_villes[n - 1];
        int first = permutation_villes[0];
        longueur_totale += gD_coord(gI_coord, last, first);
    }

    double* res = (double*)malloc(sizeof(double));
    if (!res) return NULL;
    *res = longueur_totale;
    return res;
}

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
    ALGO_CANONIQUE,
    ALGO_BF,
    ALGO_NN,
    ALGO_RW,
    ALGO_2OPTNN,
    ALGO_2OPTRW,
    ALGO_GA,
    ALGO_GADPX,
    ALGO_INCONNU
} AlgoType;

static AlgoType parse_method(const char* m){
    if (!m)              return ALGO_CANONIQUE;
    if (!strcmp(m,"bf")) return ALGO_BF;
    if (!strcmp(m,"nn")) return ALGO_NN;
    if (!strcmp(m,"rw")) return ALGO_RW;
    if (!strcmp(m,"2optnn")) return ALGO_2OPTNN;
    if (!strcmp(m,"2optrw")) return ALGO_2OPTRW;
    if (!strcmp(m,"ga"))     return ALGO_GA;
    if (!strcmp(m,"gadpx"))  return ALGO_GADPX;
    return ALGO_INCONNU;
}

/* ============================== main =============================== */

int main(int argc, char** argv){
    int opt;
    FILE* output = stdout;
    const char* outfile = NULL;
    const char* filename = NULL;
    const char* methode = NULL;
    int canonique = 0;
    int matrice = 0;
    int grande_instance = 0;

    while ((opt = getopt(argc, argv, "hf:cm:MFo:")) != -1) {
        switch (opt) {
            case 'h': usage(argv[0]); return 0;
            case 'f': filename = optarg; break;
            case 'c': canonique = 1; break;
            case 'm': methode = optarg; break;
            case 'o': outfile = optarg; break;
            case 'M': matrice = 1; break;
            case 'F': grande_instance = 1; break;
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
    if (methode){
        algo = parse_method(methode);
    } else if (canonique){
        algo = ALGO_CANONIQUE;
    } else {
        algo = ALGO_CANONIQUE;   /* par défaut */
    }

    if (algo == ALGO_INCONNU){
        fprintf(stderr, "Erreur: methode non supportee: %s\n", methode);
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

    TOUR_TSP tour = {0};
    int algo_choisi =
        (algo == ALGO_NN || algo == ALGO_RW ||
         algo == ALGO_2OPTNN || algo == ALGO_2OPTRW);

    if (algo_choisi){
        tour.DIMENSION = I.DIMENSION;
        tour.FERMEE    = 1;
        tour.LONGUEUR  = -1.0;
        tour.SECTION_TOUR = (int*)malloc((size_t)I.DIMENSION * sizeof(int));
        if (!tour.SECTION_TOUR){
            fprintf(stderr, "alloc tour travail\n");
            liberer_instance(&I);
            free(canon.SECTION_TOUR);
            return 4;
        }
    }

    TOUR_TSP meilleur = {0}, pire = {0};
    TOUR_TSP ga = {0};

    double L = -1.0;
    double secs = 0.0;
    int status = 0;

    clock_t t0 = clock();

    switch (algo){

        case ALGO_CANONIQUE: {
            L = longueur_tour(&I, &canon, d);
            break;
        }

        case ALGO_NN: {
            if (matrice){
                MatriceTSP* M = creer_matrice_demie(&I, d);
                if (!M) L = plus_proche_voisin(&I, d, &tour, 1);
                else {
                    L = plus_proche_voisin_matrice(M, &tour, 1);
                    detruire_matrice_demie(M);
                }
            } else {
                L = plus_proche_voisin(&I, d, &tour, 1);
            }
            break;
        }

        case ALGO_RW: {
            L = marche_aleatoire(&I, d, &tour);
            break;
        }

        case ALGO_2OPTRW: {
            marche_aleatoire(&I, d, &tour);
            tour.FERMEE = 1;
            L = two_opt(&I, d, &tour);
            break;
        }

        case ALGO_2OPTNN: {
            plus_proche_voisin(&I, d, &tour, 1);
            L = two_opt(&I, d, &tour);
            break;
        }

        case ALGO_GA: {
            int taille_pop = 10;
            int generations = 1000;
            double taux_mutation = 0.10;

            if (optind + 2 < argc) {
                int tmp_gen = atoi(argv[optind]);
                double tmp_mut = atof(argv[optind+1]);
                int tmp_pop = atoi(argv[optind+2]);
                if (tmp_pop > 0) taille_pop = tmp_pop;
                if (tmp_gen > 0) generations = tmp_gen;
                if (tmp_mut >= 0.0 && tmp_mut <= 1.0) taux_mutation = tmp_mut;
            }

            ga.SECTION_TOUR = NULL;
            L = tsp_ga_light(&I, d, taille_pop, generations, taux_mutation, &ga);
            if (L < 0.0){
                fprintf(stderr, "Erreur: GA a echoue.\n");
                status = 9;
            }
            break;
        }

        case ALGO_GADPX: {
            int taille_pop = 10;
            int generations = 1000;
            double taux_mutation = 0.10;

            if (optind + 2 < argc) {
                int tmp_gen = atoi(argv[optind]);
                double tmp_mut = atof(argv[optind+1]);
                int tmp_pop = atoi(argv[optind+2]);
                if (tmp_pop > 0) taille_pop = tmp_pop;
                if (tmp_gen > 0) generations = tmp_gen;
                if (tmp_mut >= 0.0 && tmp_mut <= 1.0) taux_mutation = tmp_mut;
            }

            ga.SECTION_TOUR = NULL;
            L = tsp_ga_dpx(&I, d, taille_pop, generations, taux_mutation, &ga);
            if (L < 0.0){
                fprintf(stderr, "Erreur: GA+DPX a echoue.\n");
                status = 9;
            }
            break;
        }

        case ALGO_BF: {
            if (I.DIMENSION > 12 && !grande_instance){
                fprintf(stderr, "Erreur: DIMENSION=%d > 12. Utilisez -F pour forcer la brute force.\n",
                        I.DIMENSION);
                status = 6;
                break;
            }
            int N = I.DIMENSION;
            meilleur.DIMENSION = N;
            meilleur.SECTION_TOUR = (int*) malloc((size_t)N*sizeof(int));
            pire.SECTION_TOUR = NULL;
            unsigned long long cout_meilleur_long = 0ULL;
            if (matrice){
                MatriceTSP* M = creer_matrice_demie(&I, d);
                if (M) {
                    brute_set_matrice(M);  /* définie dans tsp_force_brute_generique.c */
                    L = brute(N, 0, meilleur.SECTION_TOUR, &cout_meilleur_long, cout_tsp_matrice);
                    detruire_matrice_demie(M);
                } else {
                    brute_set_instance(&I, d);
                    L = brute(N, 0, meilleur.SECTION_TOUR, &cout_meilleur_long, cout_tsp_coord);
                }
            } else {
                brute_set_instance(&I, d);
                L = brute(N, 0, meilleur.SECTION_TOUR, &cout_meilleur_long, cout_tsp_coord);
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
        case ALGO_CANONIQUE:
            print_line(I.NAME, "canonical", secs, L, &canon);
            break;
        case ALGO_NN:
            print_line(I.NAME, "nn", secs, L, &tour);
            break;
        case ALGO_RW:
            print_line(I.NAME, "rw", secs, L, &tour);
            break;
        case ALGO_2OPTRW:
            print_line(I.NAME, "2optrw", secs, L, &tour);
            break;
        case ALGO_2OPTNN:
            print_line(I.NAME, "2optnn", secs, L, &tour);
            break;
        case ALGO_GA:
            print_line(I.NAME, "ga", secs, L, &ga);
            break;
        case ALGO_GADPX:
            print_line(I.NAME, "gadpx", secs, L, &ga);
            break;
        case ALGO_BF:
            print_line(I.NAME, "bf", secs, L, &meilleur);
            break;
        default:
            break;
    }

fin:
    if (output != stdout) fclose(output);

    if (tour.SECTION_TOUR) free(tour.SECTION_TOUR);
    if (meilleur.SECTION_TOUR) free(meilleur.SECTION_TOUR);
    if (pire.SECTION_TOUR) free(pire.SECTION_TOUR);
    if (ga.SECTION_TOUR) free(ga.SECTION_TOUR);
    free(canon.SECTION_TOUR);
    liberer_instance(&I);

    return status;
}
