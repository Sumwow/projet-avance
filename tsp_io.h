#pragma once
#include <stdio.h>
#include "tsp_types.h"


int lire_tsplib(const char* path, TSPLIB_INSTANCE* I);


void liberer_instance(TSPLIB_INSTANCE* I);
