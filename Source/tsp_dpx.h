#ifndef TSP_DPX_H
#define TSP_DPX_H

#include "tsp_types.h"
#include "tsp_distance.h"
#include "tsp_2opt.h"

double tsp_ga_dpx(const TSPLIB_INSTANCE* I,
                    DistanceFn d,
                    int population_size,
                    int generations,
                    double mutation_rate,
                    TOUR_TSP* best_tour);
#endif