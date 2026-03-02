#ifndef _TRIANGULACJA_H_
#define _TRIANGULACJA_H_

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
  pkt* punkty; 
  int l_pkt;
  link* linki;
  int l_link;
} graf;

//trzeba dodac jeszcze alokowanie pamieci

double odl (double x1, double y1, double x2, double y2);
void oblicz(graf* g, int a, int b, int c, int ac, int bc);

#endif
