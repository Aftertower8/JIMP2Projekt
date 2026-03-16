#ifndef _GRAF_H_
#define _GRAF_H_

#include "triangulacja.h"

typedef struct pkt { 
 double x; 
 double op_x;
 double y;
 double op_y;
 int n; 
} pkt; 

typedef struct link { 
  int a; 
  int b; 
  char name;
  double waga;
} link; 

typedef struct graf { 
  pkt** punkty; 
  int l_pkt;
  link* linki;
  int l_l;
} graf;

graf* load_graf (graf* g);
graf* create_graf (pkt* p, int pl, link* l, int ll);
pkt* create_pkt (double x, double op_x, double y, double op_y, int n);
link* create_link (int a, int b, char name, double waga);

#endif
