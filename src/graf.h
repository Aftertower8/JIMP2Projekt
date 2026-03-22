#ifndef _GRAF_H_ 

#define _GRAF_H_ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct pkt {  
    double x;  
    double op_x; 
    double y; 
    double op_y; 
    int n;  
} pkt;  

typedef struct link {  
    char name[50]; 
    int a;  
    int b; 
    double waga; 
} link;  

typedef struct graf {  
    pkt* punkty;  
    int l_pkt; 
    link* linki; 
    int l_l; 

} graf; 

graf* load_graf (char* argv); 
graf* create_graf (); 
int add_link(graf* g, const char* name, int a, int b, double waga);
int add_pkt(graf* g, int n);
int find_pkt(graf* g, int n);
void free_graf(graf* g);

#endif 
