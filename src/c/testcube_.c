#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "array.h"
// -----------------------------------------------------------------------------
// gcc -Wall -pedantic -std=c99 -O2 testcube_.c array.h xmalloc.c -o testcube_
// ./testcube_
// -----------------------------------------------------------------------------
int main(void){
int ***cubi;
make_cube(cubi,2,3,4);

// double **cubi;
// make_matrix(cubi,2,3);

cubi[0][0][0]=1;
cubi[1][0][0]=2;

cubi[0][1][0]=3;
cubi[1][1][0]=4;

cubi[0][2][0]=5;
cubi[1][2][0]=6;
// .
cubi[0][0][1]=7;
cubi[1][0][1]=8;

cubi[0][1][1]=9;
cubi[1][1][1]=10;

cubi[0][2][1]=11;
cubi[1][2][1]=12;
// .
cubi[0][0][2]=13;
cubi[1][0][2]=14;

cubi[0][1][2]=15;
cubi[1][1][2]=16;

cubi[0][2][2]=17;
cubi[1][2][2]=18;
// .
cubi[0][0][3]=19;
cubi[1][0][3]=20;

cubi[0][1][3]=21;
cubi[1][1][3]=22;

cubi[0][2][3]=23;
cubi[1][2][3]=24;

print_cube(" %i ",cubi,2,3,4);
free_cube(cubi);

// print_matrix(" %2.2f ",cubi,2,3);
// free_matrix(cubi);

return 0;
}
// -----------------------------------------------------------------------------

