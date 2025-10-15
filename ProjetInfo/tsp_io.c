#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "tsp_io.h"

static char* trim(char* s) {
    while (*s && isspace((unsigned char)*s)) ++s;
    size_t n = strlen(s);
    while (n && isspace((unsigned char)s[n-1])) s[--n] = '\0';
    return s;
}

int lire_tsplib(const char* path, TSPLIB_INSTANCE* I) {
    memset(I, 0, sizeof(*I));

    FILE* f = fopen(path, "r");
    if (!f) { perror("fopen"); return -1; }

    char line[512];
    int in_coords = 0;
    int coord_count = 0;

    while (fgets(line, sizeof(line), f)) {
        char* s = trim(line);
        if (*s == '\0') continue;
        if (strncmp(s, "EOF", 3) == 0) break;

        if (!in_coords) {
            if (strncmp(s, "NAME", 4) == 0) {
                char* p = strchr(s, ':'); if (p) strncpy(I->NAME, trim(p+1), sizeof(I->NAME)-1);
            } else if (strncmp(s, "COMMENT", 7) == 0) {
                char* p = strchr(s, ':'); if (p) strncpy(I->COMMENT, trim(p+1), sizeof(I->COMMENT)-1);
            } else if (strncmp(s, "TYPE", 4) == 0) {
                char* p = strchr(s, ':'); if (p) strncpy(I->TYPE, trim(p+1), sizeof(I->TYPE)-1);
            } else if (strncmp(s, "DIMENSION", 9) == 0) {
                char* p = strchr(s, ':'); if (p) I->DIMENSION = (int)strtol(trim(p+1), NULL, 10);
            } else if (strncmp(s, "EDGE_WEIGHT_TYPE", 16) == 0) {
                char* p = strchr(s, ':'); if (p) strncpy(I->EDGE_WEIGHT_TYPE, trim(p+1), sizeof(I->EDGE_WEIGHT_TYPE)-1);
            } else if (strncmp(s, "NODE_COORD_SECTION", 19) == 0) {
                if (I->DIMENSION <= 0) { fclose(f); return -2; }
                I->NODE_COORD_SECTION = (NODE*)calloc((size_t)I->DIMENSION, sizeof(NODE));
                if (!I->NODE_COORD_SECTION) { fclose(f); return -3; }
                in_coords = 1;
            }
        } else {
            if (strncmp(s, "EOF", 3) == 0) break;
            if (strstr(s, "SECTION")) continue; /* sécurité */
            int id; double x,y;
            if (sscanf(s, "%d %lf %lf", &id, &x, &y) == 3) {
                if (id < 1 || id > I->DIMENSION) { fclose(f); return -4; }
                I->NODE_COORD_SECTION[id-1].ID = id;
                I->NODE_COORD_SECTION[id-1].X = x;
                I->NODE_COORD_SECTION[id-1].Y = y;
                coord_count++;
            }
        }
    }

    fclose(f);

    if (coord_count != I->DIMENSION) return -5;
    if (strcmp(I->TYPE, "TSP") != 0) return -6;
    return 0;
}

void liberer_instance(TSPLIB_INSTANCE* I) {
    if (!I) return;
    free(I->NODE_COORD_SECTION);
    I->NODE_COORD_SECTION = NULL;
}
