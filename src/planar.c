#include "planar.h"
#define min(a, b) (((a)<(b)) ? (a) : (b))

int checkPlanar(graf g){
    if(g.l_l > 3 * g.l_pkt - 6)
        return -1;
}

