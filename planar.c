#include "planar.h"
int checkPlanar(graf g){
    if(g.l_link > 3 * g.l_pkt - 6)
        return -1;
}