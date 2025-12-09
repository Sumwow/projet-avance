#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "tsp_2opt.h"
#include "tsp_ga.h"
#include <float.h> // Pour DBL_MAX

static double calculLongueur(const TSPLIB_INSTANCE* I, DistanceFn d, const int* perm, int n)
{
    TOUR_TSP T;
    T.DIMENSION     = n;
    T.SECTION_TOUR  = (int*)perm;
    T.FERMEE        = 1;
    T.LONGUEUR      = -1.0;
    return longueur_tour(I, &T, d);
}

static void randomPermutation(int* perm, int n)
{
    for (int i = 0; i < n; ++i)
        perm[i] = i + 1;

    for (int i = n - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        int tmp = perm[i];
        perm[i] = perm[j];
        perm[j] = tmp;
    }
}

static void swap_mutation(int* perm, int n, double mutation_rate)
{
    if (n < 2) return;

    double r = (double)rand() / (double)RAND_MAX;
    if (r >= mutation_rate) return;

    int i = rand() % n;
    int j = rand() % n;
    while (j == i)
        j = rand() % n;

    int tmp = perm[i];
    perm[i] = perm[j];
    perm[j] = tmp;
}

static void tournament_selection(const TSPLIB_INSTANCE* I, 
                                 DistanceFn d, 
                                 int** population,
                                 int pop_size,
                                 int** selected,
                                 int tournament_size,
                                 int n)
{
    if (tournament_size < 1) tournament_size = 1;

    for (int i = 0; i < pop_size; ++i) {
        int best_idx = rand() % pop_size;
        double best_fit = calculLongueur(I, d, population[best_idx], n);

        for (int t = 1; t < tournament_size; ++t) {
            int idx = rand() % pop_size;
            double f = calculLongueur(I, d, population[idx], n);
            if (f < best_fit) {
                best_fit = f;
                best_idx = idx;
            }
        }
        for (int k = 0; k < n; ++k)
            selected[i][k] = population[best_idx][k];
    }
}

static void triLongueur(const TSPLIB_INSTANCE* I,
                            DistanceFn d,
                            int** offspring,
                            int pop_size,
                            int n)
{
    for (int i = 0; i < pop_size - 1; ++i) {
        int best_idx = i;
        double best_fit = calculLongueur(I, d, offspring[i], n);
        for (int j = i + 1; j < pop_size; ++j) {
            double f = calculLongueur(I, d, offspring[j], n);
            if (f < best_fit) {
                best_fit = f;
                best_idx = j;
            }
        }
        if (best_idx != i) {
            int* tmp = offspring[i];
            offspring[i] = offspring[best_idx];
            offspring[best_idx] = tmp;
        }
    }
}


// Structure utilitaire pour gérer les fragments
// next[i] = ville suivante après i, ou -1 si i est une fin de fragment
// prev[i] = ville précédente avant i, ou -1 si i est un début de fragment
static void dpx_crossover(const TSPLIB_INSTANCE* I,
                          DistanceFn d,
                          const int* p1,
                          const int* p2,
                          int* child,
                          int n)
{
    int* next = (int*)malloc(n * sizeof(int));
    int* prev = (int*)malloc(n * sizeof(int));

    int* visited = (int*)calloc(n, sizeof(int)); 

    for(int i=0; i<n; ++i) {
        next[i] = -1;
        prev[i] = -1;
    }

    int (*neighbors_p2)[2] = malloc(n * sizeof(*neighbors_p2));
    for (int i = 0; i < n; ++i) {
        neighbors_p2[i][0] = -1;
        neighbors_p2[i][1] = -1;
    }
    
    for (int i = 0; i < n; ++i) {
        int u = p2[i];
        int v_prev = p2[(i - 1 + n) % n];
        int v_next = p2[(i + 1) % n];
        neighbors_p2[u-1][0] = v_prev;
        neighbors_p2[u-1][1] = v_next;
    }

    for (int i = 0; i < n; ++i) {
        int u = p1[i];
        int v = p1[(i + 1) % n];

        int u_idx = u - 1;
        int is_common = (neighbors_p2[u_idx][0] == v) || (neighbors_p2[u_idx][1] == v);

        if (is_common) {
            next[u_idx] = v;
            prev[v - 1] = u;
        }
    }

    int start_node = -1;
    for(int i=0; i<n; ++i) {
        if (prev[i] == -1) {
            start_node = i + 1; // Ville 1..N
            break;
        }
    }

    if (start_node == -1) start_node = p1[0];

    int current = start_node;
    int count = 0;

    while (count < n) {
        child[count++] = current;
        visited[current - 1] = 1;

        if (next[current - 1] != -1) {
            current = next[current - 1];
        } else {
            
            double min_dist = DBL_MAX;
            int best_k = -1;

            for (int k_idx = 0; k_idx < n; ++k_idx) {
                int k_city = k_idx + 1;
              
                if (!visited[k_idx] && prev[k_idx] == -1) {

                    double dist = d(I, current, k_city);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_k = k_city;
                    }
                }
            }

            if (best_k == -1) {
                // Normalement cela n'arrive qu'à la toute fin
                break; 
            }
            current = best_k;
        }
    }

    free(next);
    free(prev);
    free(visited);
    free(neighbors_p2);
}

double tsp_ga_dpx(const TSPLIB_INSTANCE* I,
                    DistanceFn d,
                    int population_size,
                    int generations,
                    double mutation_rate,
                    TOUR_TSP* best_tour)
{
if (!I || !d || !best_tour) return -1.0;

    int n = I->DIMENSION;
    if (n <= 0 || population_size <= 0 || generations <= 0)
        return -1.0;
    static int seeded = 0;
    if (!seeded) {
        srand((unsigned int)time(NULL));
        seeded = 1;
    }
    
    int** population = (int**)malloc((size_t)population_size * sizeof(int*));
    int** selected   = (int**)malloc((size_t)population_size * sizeof(int*));
    int** offspring  = (int**)malloc((size_t)population_size * sizeof(int*));
    if (!population || !selected || !offspring) {
        free(population); free(selected); free(offspring);
        return -1.0;
    }

    for (int i = 0; i < population_size; ++i) {
        population[i] = (int*)malloc((size_t)n * sizeof(int));
        selected[i]   = (int*)malloc((size_t)n * sizeof(int));
        offspring[i]  = (int*)malloc((size_t)n * sizeof(int));
        if (!population[i] || !selected[i] || !offspring[i]) {
            for (int k = 0; k <= i; ++k) {
                free(population[k]);
                free(selected[k]);
                free(offspring[k]);
            }
            free(population); free(selected); free(offspring);
            return -1.0;
        }
    }
    for (int i = 0; i < population_size; ++i)
        randomPermutation(population[i], n);
        
    int* best_perm = (int*)malloc((size_t)n * sizeof(int));
    if (!best_perm) {
        for (int i = 0; i < population_size; ++i) {
            free(population[i]); free(selected[i]); free(offspring[i]);
        }
        free(population); free(selected); free(offspring);
        return -1.0;
    }

    double best_len = calculLongueur(I, d, population[0], n);
    for (int k = 0; k < n; ++k)
        best_perm[k] = population[0][k];

    int tournament_size = (int)(0.5 * population_size);
    if (tournament_size < 2) tournament_size = 2;

    for (int g = 0; g < generations; ++g) {
        tournament_selection(I, d, population, population_size, selected, tournament_size, n);
    tournament_selection(I, d, population, population_size, selected, tournament_size, n);

    for (int i = 0; i < population_size; i += 2) {
        int next = (i + 1 < population_size) ? (i + 1) : 0;

        dpx_crossover(I, d, selected[i], selected[next], offspring[i], n);

        {
            TOUR_TSP temp_tour;
            temp_tour.DIMENSION = n;
            temp_tour.FERMEE = 1;
            temp_tour.SECTION_TOUR = offspring[i];
            temp_tour.LONGUEUR = -1.0;

            two_opt(I, d, &temp_tour);

        }

        if (i + 1 < population_size) {

            dpx_crossover(I, d, selected[next], selected[i], offspring[i + 1], n);

            {
                TOUR_TSP temp_tour;
                temp_tour.DIMENSION = n;
                temp_tour.FERMEE = 1;
                temp_tour.SECTION_TOUR = offspring[i + 1];
                temp_tour.LONGUEUR = -1.0;
                
                two_opt(I, d, &temp_tour);
            }
        }
    }
    
        for (int i = 0; i < population_size; ++i)
            swap_mutation(offspring[i], n, mutation_rate);
            
        triLongueur(I, d, offspring, population_size, n);

        for (int i = 0; i < population_size; ++i)
            memcpy(population[i], offspring[i], (size_t)n * sizeof(int));

        int worst_idx = population_size - 1;
        randomPermutation(population[worst_idx], n);

        int best_gen_idx = 0;
        int worst_gen_idx = 0;
        double best_gen_fit = calculLongueur(I, d, population[0], n);
        double worst_gen_fit = best_gen_fit;

        for (int i = 1; i < population_size; ++i) {
            double f = calculLongueur(I, d, population[i], n);
            if (f < best_gen_fit) {
                best_gen_fit = f;
                best_gen_idx = i;
            }
            if (f > worst_gen_fit) {
                worst_gen_fit = f;
                worst_gen_idx = i;
            }
        }
        if (best_gen_fit < best_len) {
            best_len = best_gen_fit;
            for (int k = 0; k < n; ++k)
                best_perm[k] = population[best_gen_idx][k];
        }

        for (int k = 0; k < n; ++k)
            population[worst_gen_idx][k] = best_perm[k];
    }
    
    if (!best_tour->SECTION_TOUR) {
        best_tour->SECTION_TOUR = (int*)malloc((size_t)n * sizeof(int));
        if (!best_tour->SECTION_TOUR) {
            best_len = -1.0;
        }
    }
    if (best_tour->SECTION_TOUR) {
        best_tour->DIMENSION = n;
        best_tour->FERMEE    = 1;
        best_tour->LONGUEUR  = best_len;
        for (int i = 0; i < n; ++i)
            best_tour->SECTION_TOUR[i] = best_perm[i];
    }

    free(best_perm);
    for (int i = 0; i < population_size; ++i) {
        free(population[i]);
        free(selected[i]);
        free(offspring[i]);
    }
    free(population);
    free(selected);
    free(offspring);

    return best_len;
}