#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    char *filename = NULL;
    int calcule_canonique = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            filename = argv[i + 1];
            i++;

        } else if (strcmp(argv[i], "-c") == 0) {
            calcule_canonique = 1;
        }
    }

    if (filename == NULL) {
        printf("Usage : %s -f <fichier> [-c]\n", argv[0]);
        return 1;
    }

    printf("Fichier à lire : %s\n", filename);

    if (calcule_canonique) {
        printf("Calcul de la tournée canonique...\n");
        // plus tard : lecture du fichier + calcul
        
    } else {
        printf("Aucun calcul demandé.\n");
    }

    return 0;
}