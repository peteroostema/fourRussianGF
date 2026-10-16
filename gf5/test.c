/* this is just a sample row computation I did to varify my testing routine
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */


#include <stdio.h>
#include <stdlib.h>

int main(){
   int result = 1*2 + 2*4 + 2*2 + 4*2 + 3*2 + 2*4 + 3*0 + 4*2 + 2*1 + 3*0 + 0*3 + 0*0 + 4*1 + 3*2 + 2*3 + 1*4;
   int resul2 = 0*1 + 4*1 + 0*4 + 2*4 + 2*2 + 3*1 + 0*3 + 4*3 + 2*4 + 2*0 + 3*1 + 0*0 + 1*1 + 1*4 + 1*1 + 1*0;
   printf("%d\n%d\n", result, resul2);
}