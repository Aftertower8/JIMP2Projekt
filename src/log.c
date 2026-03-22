#include <stdio.h>
#include "log.h"
int ver=0;
void verbose(char *message){
    if(ver==1)
        printf("%s",message);
}
void verbose_activate(){
    ver=1;
}
