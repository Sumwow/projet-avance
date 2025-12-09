#ifndef TSP_NN_H
#define TSP_NN_H

#include "tsp_types.h"
#include "tsp_distance.h"
#include "tsp_matrice.h"

double plus_proche_voisin(const TSPLIB_INSTANCE* I,
                          DistanceFn d,
                          TOUR_TSP* tour,
                          int depart);

double plus_proche_voisin_matrice(const MatriceTSP* M,
                                  TOUR_TSP* tour,
                                  int depart);

#endif
