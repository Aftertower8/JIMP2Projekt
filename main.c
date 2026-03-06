#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pob_dane.h"

int main(int argc, char** argv)
{
  if(argc == 1){ //Sprawdzanie ilosci argumentow
      fprintf(stderr,"Nie podano zadnych argumentow.\nProgram wymaga podania przynajmniej sciezki do pliku z grafem.\n");
      return 1;
  }

  wczyt_graf(argv[1]); //wczytanie pliku z grafem
  char a[50] = "triangulacja"; //Bazowo bez podania flagi argumentu bedzie triangulacja
  FILE* out = NULL;
  if (pobierz_dane(argc, argv, a, &out) != 0) {
        help();
        return 2;
  }

  //wlaczanie funkcji posrednich
  //do zrobienia
  //zapisywanie wyniku
  //do zrobienia
  
  if(out) fclose(out);
  return 0;
}
