/* negation.c has a few different operations for GF5 commented out
 * examples: doubling, negation, substraction, addition
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printBinary(long unsigned int number){
   long unsigned int tmp = number;
   for (int k = 0; k < 64; k++){
      if (tmp & 0x8000000000000000){
         printf("1");
      } else {
         printf("0");
      }
      tmp <<= 1;
   }
   printf("\n");
}

void printGF5(long unsigned int r2, long unsigned int r1, long unsigned int r0){
   for (int i = 63; i >= 0; i--){
      long unsigned int tmp0 = (r0 >> i) & 1;
      long unsigned int tmp1 = (r1 >> i) & 1;
      long unsigned int tmp2 = (r2 >> i) & 1;
//      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("%lu", tmp2*4 + tmp1*2 + tmp0);
   }
   printf("\n");
}

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

int main(){
   srand(time(NULL));
   long unsigned int a0, a1, a2, b0, b1, b2, r0, r1, r2, s0, s1, s2, s3, t;
   
   printf("a\n");
   for (int i = 0; i < 64; i++){
      long unsigned int tmp = rand() % 5;
      printf("%lu", tmp);
      a0 <<= 1;
      a1 <<= 1;
      a2 <<= 1;
      if (tmp & 0b100){
         a2 |= 1;
      }
      if (tmp & 0b010){
         a1 |= 1;
      }
      if (tmp & 0b001){
         a0 |= 1;
      }
   }
   printf("\n");
   printBinary(a2);
   printBinary(a1);
   printBinary(a0);
   printf("\n");
   
   printf("b\n");
   for (int i = 0; i < 64; i++){
      long unsigned int tmp = rand() % 5;
      printf("%lu", tmp);
      b0 <<= 1;
      b1 <<= 1;
      b2 <<= 1;
      if (tmp & 0b100){
         b2 |= 1;
      }
      if (tmp & 0b010){
         b1 |= 1;
      }
      if (tmp & 0b001){
         b0 |= 1;
      }
   }
//   printf("\n");
//   printBinary(b2);
//   printBinary(b1);
//   printBinary(b0);
//   printf("\n");
   
   unsigned long long startTime = rdtsc();
//   s0 = a0 ^ b0;// bit packed addition
//   s1 = a1 ^ b1;
//   s2 = a2 ^ b2;
//   s2 ^= s1 & (a0 & b0);
////   s3 ^= s2 & (a1 & b1);
//   s1 ^= a0 & b0;
//   s2 ^= a1 & b1;
//   s3 = a2 & b2;
   
//   t = s2 | s1;// paper fast method (incorrect)
//   r2 = s0 ^ t;
//   r1 = (r2 & s0) ^ (s3 ^ s1);
//   r0 = (t ^ s2) | (r1 & s3);
//   t = s1 | s0;
   
//   unsigned long int n0, n1, n2, e3, e2, e1, e0;// paper slower meathod
//   n0 = ~s2; n1 = s3; n2 = ~s3;
//   e0 = n0 ^ s0;
//   e1 = n1 ^ s1;
//   e2 = n2;
//   e2 ^= (e1 & (n0 & s0)) | (n1 & s1);
//   e3 = n2 & n0 & s1 & s0;
//   e1 ^= n0 & s0;
//   r0 = e0 & e3;
//   r1 = e1 & e3;
//   r2 = e2;
   
//   r2 = (s2 & ~s1 & ~s0); //my meathod
//   r1 = (a1) | (s1 & s0) | (~s2 & s1);
//   r0 = (~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2);
   
   r2 = (~a1 & a0); //my meathod
   r1 = a1;
   r0 = (a2 | a1) & ~a0;
   
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
//   printf("s\n");
//   printGF5(s2, s1, s0);
//   printBinary(s3);
//   printBinary(s2);
//   printBinary(s1);
//   printBinary(s0);
//   printf("\n");
   
//      printf("n\n");
//      printGF5(n2, n1, n0);
//      printBinary(n2);
//      printBinary(n1);
//      printBinary(n0);
//      printf("\n");
//      printf("e\n");
//      printGF5(e2, e1, e0);
//      printBinary(e3);
//      printBinary(e2);
//      printBinary(e1);
//      printBinary(e0);
//      printf("\n");
   
   printf("r\n");
   printGF5(r2, r1, r0);
   printBinary(r2);
   printBinary(r1);
   printBinary(r0);
   printf("\n");
   

}//wrong 6 > 6, 2 > 7, 1 > 6

