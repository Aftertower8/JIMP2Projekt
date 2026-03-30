#ifndef _TRIANGULACJA_H_
#define _TRIANGULACJA_H_

//trzeba dodac jeszcze alokowanie pamieci

double odl (double x1, double y1, double x2, double y2);
void oblicz(graf* g, int a, int b, int c, int ac, int bc);
int find_pkt_pos(graf* g);
int find_link(graf* g, a, n);
void triangulation(graf* g);
#endif
