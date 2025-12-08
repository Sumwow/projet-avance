#ifndef TSP_IO_H
#define TSP_IO_H

#include <stdio.h>
#include "tsp_types.h"


int lire_tsplib(const char* path, TSPLIB_INSTANCE* I);


void liberer_instance(TSPLIB_INSTANCE* I);

#endif