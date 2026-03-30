#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "triangulation.h"

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

  if(find_pkt_pos(g, x1f,y1f)){
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

void triangulation(graf* g){
  g->punkty[g->linki[0].b].x == g->linki[0].waga;
  int k = 0;
  for(int i = 1; i<g->l_l; i++){
    if(find_link(g, g->linki[k].a, g->linki[i].a) != -1 )
      for(int j = 1; j<g->l_l; j++){
        if(find_link(g, g->linki[k].b, g->linki[j].a) != -1){
          oblicz(g, g->linki[k]a, g->linki[i-1].b, g->linki[i].a, g->linki[i], g->linki[j]); //od indexu linku do którego znalazło ten sam punkt w obu
          k++  
          break;
        }
        else if(find_link(g, g->linki[k].b, g->linki[j].a) != -1){
          oblicz(g, g->linki[k]a, g->linki[i-1].b, g->linki[i].a, g->linki[i], g->linki[j]);
          k++  
          break;
        }
      }

    else if(find_link(g, g->linki[k].a, g->linki[i].b) != -1 )
      for(int j = 1; j<g->l_l; j++){
        if(find_link(g, g->linki[k].b, g->linki[j].b) != -1){
          oblicz(g, g->linki[k].a, g->linki[i-1].b, g->linki[i].b, g->linki[i], g->linki[j]); //g->linki[i].b/a to index obecnego punktu trzeciego który znaleźliśmym do niego wsadza sie wynik
          k++;  
          break;
        }
        else if(find_link(g, g->linki[k].b, g->linki[j].b) != -1){
          oblicz(g, g->linki[k].a, g->linki[i-1].b, g->linki[i].b, g->linki[i], g->linki[j]); 
          k++;  
          break;
        }
      }
    }
      if(find_pkt_pos(g, g->punkty[i].x, g->punkty[i].y) == -1)
        //backtrack?
}

int find_pkt_pos(graf* g, double x, double y){
  for(int i = 0; i < g->l_pkt; i++){
    if(g->punkty[i].x == x && g->punkty[i].y == y)
      return g->punkty[i].n;
  return -1;
  }
}

int find_link(graf* g, a, n){
  for(int i = 0; i < g->l_l; i++){
    if((g->linki[i].a == a && g->linki[i].b == n) || (g->linki[i].a == n && g->linki[i].b == a))
      return i;
  return -1;
  }
}

//trzeba dodac funkcje iterujaca po punktach funkcją oblicz i w razie potrzeby backtrackujaca sie jesli punkty zaczynaja nachodzic na siebie
