/* timeAproximations.c gives an order of magnitude approximation of
 * how long each field needs for a four Russian matrix multiplication
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){
   long int m = 16384;
   long int k = 16384;
   long int n = 16384;
   long int t = 4;
   
   long int r4 = (k/64)*((n/64) * (64/t) + (1 << t)*(0.5*t + (n/64) * (0.5*t) * (64/t)) + m*(4 + (n/64)*4));
   long int reg = m*n*k/64;
   printf("GF2\n4Rus:    %ld \nregular: %ld\n", r4, reg);
   
   long int r45 = (k/16)*((n/16) * (16/t) + pow(5, t)*(0.5*t + (n/16) * (0.5*t) * (16/t)*3) + m*(4*4 + (n/16)*4*3));
   long int reg5 = (m*n*k/16)*3;
   printf("GF5\n4Rus:    %ld \nregular: %ld\n", r45, reg5);
   
   long int r47 = (k/16)*((n/16) * (16/t) + pow(7, t)*(0.5*t + (n/16) * (0.5*t) * (16/t)*3) + m*(4*4 + (n/16)*4*3));
   long int reg7 = (m*n*k/16)*3;
   printf("GF7\n4Rus:    %ld \nregular: %ld\n", r47, reg7);
}