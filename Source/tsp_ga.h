#ifndef TSP_GA_H
#define TSP_GA_H

#include "tsp_types.h"
#include "tsp_distance.h"

double tsp_ga_light(const TSPLIB_INSTANCE* I,
                    DistanceFn d,
                    int population_size,
                    int generations,
                    double mutation_rate,
                    TOUR_TSP* best_tour);

#endif

