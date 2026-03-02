#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "triangulacja.h"

double odl (double x1, double y1, double x2, double y2) { //odleglosc punktow od siebie
  double dl = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  return dl;
}
void oblicz(graf* g, int a, int b, int c, int ac, int bc){ //obliczenie punktow przeciecia sie dwoch okregow o srodkach sasiadow punktu c oraz promieniu odleglosci pomiedzy nimi a punktem c

  double d = odl(g->punkty[b].x, g->punkty[b].y,
               g->punkty[a].x, g->punkty[a].y);

  if(d == 0) return;

  double ra = g->linki[ac].waga;
  double rb = g->linki[bc].waga;

  double k = (ra*ra - rb*rb + d*d) / (2*d);

  double h = sqrt(ra*ra - k*k);

  double x = g->punkty[a].x + k*(g->punkty[b].x - g->punkty[a].x)/d;
  double y = g->punkty[a].y + k*(g->punkty[b].y - g->punkty[a].y)/d;

  double x1f = x + h*(g->punkty[b].y - g->punkty[a].y)/d;
  double y1f = y - h*(g->punkty[b].x - g->punkty[a].x)/d;

  double x2f = x - h*(g->punkty[b].y - g->punkty[a].y)/d;
  double y2f = y + h*(g->punkty[b].x - g->punkty[a].x)/d;

  if(sprawdz(x1f,y1f,g)){
    g->punkty[c].x = x1f;
    g->punkty[c].y = y1f;
    g->punkty[c].op_x = x2f;
    g->punkty[c].op_y = y2f;
  }else{
    g->punkty[c].x = x2f;
    g->punkty[c].y = y2f;
    g->punkty[c].op_x = x1f;
    g->punkty[c].op_y = y1f;
  }
}
