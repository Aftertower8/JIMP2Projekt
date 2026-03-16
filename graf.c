#include "graf.h"
#include <stdio.h>
#include <stdlib.h>

pkt* create_pkt (double x, double op_x, double y, double op_y, int n){
    pkt* p = (pkt*)malloc(sizeof(pkt));

    if( p = NULL){
        p->x = (double)malloc(sizeof(double));
        p->op_x = (double)malloc(sizeof(double));
        p->y = (double)malloc(sizeof(double));
        p->op_y = (double)malloc(sizeof(double));
        p->n = (int)malloc(sizeof(int));
    }

    if( !p || !p->x || !p->op_x || !p->y || !p->op_y || !p->n){
        free(p->x);
        free(p->op_x);
        free(p->y);
        free(p->op_y);
        free(p->n);
        free(p);
    }
}

graf* create_graf (pkt* p, int pl, link* l, int ll){
    graf* g = (graf*)malloc(sizeof(graf));

    if( g = NULL){
        fprintf(stderr, "Blad alokacji pamieci dla roju.\n");
        g->l_pkt = (int)malloc(sizeof(int));
        g->l_l = (int)malloc(sizeof(int));
    }

    if(!g->linki || !g->l_l){
        free(g->linki);
        free(g->l_l);
        free(g);
        return NULL;
    }
}

void create_graf(graf* g) {
    if (!g) return;
    int i = 0;
    while(g != NULL) {
        free_pkt(g, i);
        free_link(g, i);
        i++;
    }
    free(g->l_l);
    free(g->l_pkt);
    free(g);
}

void free_link(graf* g, int i){
        free(g->linki[i]->a);
        free(g->linki[i]->b);
        free(g->linki[i]->name);
        free(g->linki[i]->waga);
        free(g->linki[i]);
}

void free_pkt(graf* g, int i){
    free(g->pkt[i]->x);
    free(g->pkt[i]->op_x);
    free(g->pkt[i]->y);
    free(g->pkt[i]->op_y);
    free(g->pkt[i]->n);
    free(g->pkt[i]);
}

link* create_link (int a, int b, char name, double waga){
    link* l = (link*)malloc(sizeof(link));

    if( l = NULL){
        l-> = (double)malloc(sizeof(double));
        l->op_x = (double)malloc(sizeof(double));
        l->y = (double)malloc(sizeof(double));
        l->op_y = (double)malloc(sizeof(double));
        l-> = (int)malloc(sizeof(int));
    }

    if( !p || !p->x || !p->op_x || !p->y || !p->op_y || !p->n){
        free(p->x);
        free(p->op_x);
        free(p->y);
        free(p->op_y);
        free(p->n);
        free(p);
    }

}
