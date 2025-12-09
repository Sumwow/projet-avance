#ifndef TSP_TYPES_H
#define TSP_TYPES_H

#ifndef TSPLIB_MAX_NAME
#define TSPLIB_MAX_NAME 128
#endif
#ifndef TSPLIB_MAX_COMMENT
#define TSPLIB_MAX_COMMENT 512
#endif

#include <stddef.h>
#include <stdint.h>

/* Structure représentant un nœud (ville) avec son identifiant et ses coordonnées */
typedef struct {
    int ID;
    double X;
    double Y;
} NODE;

/* Structure représentant une instance TSPLIB avec ses métadonnées et ses coordonnées de nœuds */
typedef struct TSPLIB_INSTANCE {
    char NAME[TSPLIB_MAX_NAME];
    char COMMENT[TSPLIB_MAX_COMMENT];
    char TYPE[8];
    int DIMENSION;
    char EDGE_WEIGHT_TYPE[16];
    NODE* NODE_COORD_SECTION;
    void* USER_DATA;
} TSPLIB_INSTANCE;

/* Structure représentant une tournée TSP avec son ordre de visite et sa longueur totale */
typedef struct TOUR_TSP {
    int DIMENSION;
    int* SECTION_TOUR;
    int FERMEE;
    double LONGUEUR;
} TOUR_TSP;

#endif