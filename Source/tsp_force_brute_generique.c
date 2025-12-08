#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <signal.h>
#include "tsp_matrice.h"


static MatriceTSP *global_matrice_tsp = NULL;

void brute_set_matrice(MatriceTSP *matrice) {
    global_matrice_tsp = matrice;
}

/*SIGNAL*/

static bool signal_recu_forcebrute = false;

void handler_forcebrute(int sig) {
    (void)sig;
    signal_recu_forcebrute = true;
}

static void print_force_brute(double meilleur, double pire) {
    printf("Meilleur longueur : %f\n", meilleur);
    printf("Pire longueur : %f\n", pire);
    printf("Entrez un caractère pour reprendre\n");
    getchar();
    printf("Reprise de force brute\n");
    signal_recu_forcebrute = false;
}

/*PERMUTATION*/

static int prochaine_permutation(int *villes, int nbVilles) {
    int i = nbVilles - 2;

    while (i >= 0 && villes[i] >= villes[i + 1]) i--;
    if (i < 0) return 0;   

    int j = nbVilles - 1;
    while (villes[j] <= villes[i]) j--;

    int tmp = villes[i];
    villes[i] = villes[j];
    villes[j] = tmp;

    int d = i + 1, f = nbVilles - 1;
    while (d < f) {
        tmp = villes[d];
        villes[d] = villes[f];
        villes[f] = tmp;
        d++;
        f--;
    }
    return 1;
}

/*DISTANCE MATRICE*/

static double lire_distance_matrice(const MatriceTSP* matrice, int i, int j) {
    if (i == j) return 0.0;

    if (i > j) { 
        int t = i; 
        i = j; 
        j = t; 
    }

    int ligne = i - 1;
    int colonne = j - i - 1;

    return matrice->data[ligne][colonne];
}

/*FONCTION COUT*/

void *cout_tsp_matrice(void *param_inutile, int *permutation_villes) {
    (void)param_inutile; 

    
    if (!global_matrice_tsp)
        return NULL;

    int nombre_villes = global_matrice_tsp->dimension;
    double longueur_totale = 0.0;

    for (int indice = 0; indice < nombre_villes - 1; indice++) {
        int ville_depart = permutation_villes[indice];
        int ville_arrivee = permutation_villes[indice + 1];
        longueur_totale += lire_distance_matrice(global_matrice_tsp, ville_depart, ville_arrivee);
    }

    if (nombre_villes > 1) {
        int derniere_ville = permutation_villes[nombre_villes - 1];
        int premiere_ville = permutation_villes[0];
        longueur_totale += lire_distance_matrice(global_matrice_tsp, derniere_ville, premiere_ville);
    }

    double *resultat = (double*)malloc(sizeof(double));
    if (!resultat) return NULL;

    *resultat = longueur_totale;

    return resultat;
}

/*BRUTE FORCE GENERIQUE*/

double brute(int nombre_noeuds, int nombre_ressources, int *meilleure_permutation, unsigned long long *cout_meilleur, void * (*fonction_cout)(void *, int *)) {

    int *ordre_courant = (int*)malloc((size_t)nombre_noeuds * sizeof(int));
    if (!ordre_courant) {
        fprintf(stderr, "Erreur malloc brute\n");
        return -1.0;
    }

    for (int indice = 0; indice < nombre_noeuds; indice++)
        ordre_courant[indice] = indice + 1;

    signal(SIGINT, handler_forcebrute);

    void *valeur_cout_brute = fonction_cout(NULL, ordre_courant);
    if (!valeur_cout_brute) {
        free(ordre_courant);
        return -1.0;
    }

    double meilleur_cout_double;
    unsigned long long meilleur_cout_long;

    if (nombre_ressources == 0) {
        double *cout_temporaire = (double*)valeur_cout_brute;
        meilleur_cout_double = *cout_temporaire;
        meilleur_cout_long = (unsigned long long)(*cout_temporaire);
    } else {
        unsigned long long *cout_temporaire = (unsigned long long*)valeur_cout_brute;
        meilleur_cout_long = *cout_temporaire;
        meilleur_cout_double = (double)meilleur_cout_long;
    }
    free(valeur_cout_brute);

    for (int indice = 0; indice < nombre_noeuds; indice++)
        meilleure_permutation[indice] = ordre_courant[indice];

    *cout_meilleur = meilleur_cout_long;

    double pire_cout_double = meilleur_cout_double;
    unsigned long long pire_cout_long = meilleur_cout_long;

    while (prochaine_permutation(ordre_courant, nombre_noeuds)) {

        signal(SIGINT, handler_forcebrute);

        void *valeur_cout_courant = fonction_cout(NULL, ordre_courant);
        if (!valeur_cout_courant) {
            free(ordre_courant);
            return -1.0;
        }

        double cout_courant_double;
        unsigned long long cout_courant_long;

        if (nombre_ressources == 0) {
            double *cout_temporaire = (double*)valeur_cout_courant;
            cout_courant_double = *cout_temporaire;
            cout_courant_long = (unsigned long long)(*cout_temporaire);
        } else {
            unsigned long long *cout_temporaire = (unsigned long long*)valeur_cout_courant;
            cout_courant_long = *cout_temporaire;
            cout_courant_double = (double)cout_courant_long;
        }
        free(valeur_cout_courant);

        if (cout_courant_double < meilleur_cout_double) {
            meilleur_cout_double = cout_courant_double;
            meilleur_cout_long = cout_courant_long;

            for (int indice = 0; indice < nombre_noeuds; indice++)
                meilleure_permutation[indice] = ordre_courant[indice];

            *cout_meilleur = meilleur_cout_long;
        }

        if (cout_courant_double > pire_cout_double) {
            pire_cout_double = cout_courant_double;
            pire_cout_long = cout_courant_long;
        }

        if (signal_recu_forcebrute) {
            print_force_brute(meilleur_cout_double, pire_cout_double);
        }
    }

    free(ordre_courant);
    return meilleur_cout_double;
}
