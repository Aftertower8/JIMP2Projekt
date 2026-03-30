#include <stdio.h>
#include <stdlib.h>
#include "zapis.h"

void zapisz(char* file, double x, double y){
    FILE* out = fopen(file, "a");
    fprintf(out, "pozycje: %d, %d", x, y);
}
