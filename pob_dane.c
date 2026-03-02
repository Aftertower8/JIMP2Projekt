#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pob_dane.h"

int pobierz_dane(int argc, char** argv, char a[], FILE** out){
for (int j = 2; j < argc; j++) { //wczytywanie z flagami
    if (strcmp(argv[j], "-v") == 0) { //tryb verbose
        verbose(); //<-do zrobienia
    }
    if (strcmp(argv[j], "-h") == 0) { //tryb help
        help(); //<-do zrobienia
    }
    if (strcmp(argv[j], "-a") == 0) { //zmiana na inny algorytm
        j++;
        if (j < argc)
            if(sscanf(argv[j], "%49s", a) != 1) {
              fprintf(stderr,"Bledny argument dla -a\n");
              return 1;
            }
    }
    if (strcmp(argv[j], "-w") == 0) { //otworzenie pliku z wynikiem
        j++;
        if (j < argc){
            out = fopen(argv[j], "a");
          if (!out) {
            fprintf(stderr, "Nie mozna otworzyc pliku");
            return 1;
          }
        }
    }
    if (strcmp(argv[j], "-w") != 0 &&
        strcmp(argv[j], "-a") != 0 &&
        strcmp(argv[j], "-h") != 0 &&
        strcmp(argv[j], "-v") != 0) {
        fprintf(stderr,"Nieznana flaga: %s\n", argv[j]);
        return 1;
    }
  }
  return 0;
}
