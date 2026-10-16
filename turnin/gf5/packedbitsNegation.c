/* packedbitsNegation.c gives bitwise negations and doubling for packed
 * elements of GF5
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

void printGF5(long unsigned int r){
   for (int i = 15; i >= 0; i--){
      long unsigned int tmp = (r >> i*4) & 0xf;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("   %lu", tmp);
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
   long unsigned int a, b, d, a0, a1, a2, b0, b1, b2, r, s0, s1, s2, s3, t = 0;
   
   printf("a\n");
   for (int i = 0; i < 16; i++){
      long unsigned int tmp = rand() % 5;
      printf("   %lu", tmp);
      a <<= 4;
      a |= tmp;
   }
   printf("\n");
   printBinary(a);
   printf("\n");
   
   printf("b\n");
   for (int i = 0; i < 16; i++){
      long unsigned int tmp = rand() % 5;
      printf("   %lu", tmp);
      b <<= 4;
      b |= tmp;
   }
//   printf("\n");
//   printBinary(b2);
//   printBinary(b1);
//   printBinary(b0);
//   printf("\n");
   r = 0; d = 0;
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
   a2 = (a & 0x4444444444444444) >> 2;
   a1 = (a & 0x2222222222222222) >> 1;
   a0 = a & 0x1111111111111111;
   r |= (~a1 & a0) << 2; //my meathod
   r |= a1 << 1;
   r |= (a2 | a1) & ~a0;
   
   d |= (a1 & ~a0) << 2; //my meathod
   d |= ((~a1 & a0) | a2) << 1;
   d |= a2 | (a1 & a0);
   
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
//   printf("r2\n");
//   printGF5((~a1 & a0) << 2);
//   printBinary((~a1 & a0) << 2);
//   printf("\n");
//   printf("r2\n");
//   printGF5(a1 << 1);
//   printBinary(a1 << 1);
//   printf("\n");
//   printf("r2\n");
//   printGF5((a2 | a1) & ~a0 << 0);
//   printBinary((a2 | a1) & ~a0 << 0);
//   printf("\n");
   
//      printf("d2\n");
//      printGF5((a1 & ~a0) << 2);
//      printBinary((a1 & ~a0) << 2);
//      printf("\n");
//      printf("d1\n");
//      printGF5(((~a1 & a0) | a2) << 1);
//      printBinary(((~a1 & a0) | a2) << 1);
//      printf("\n");
//      printf("d0\n");
//      printGF5(a2 | (a1 & a0) << 0);
//      printBinary(a2 | (a1 & a0) << 0);
//      printf("\n");
//   
//   printf("a2\n");
//   printGF5(a2);
//   printBinary(a2);
//   printf("\n");
//   printf("a1\n");
//   printGF5(a1);
//   printBinary(a1);
//   printf("\n");
//   printf("a0\n");
//   printGF5(a0);
//   printBinary(a0);
//   printf("\n");
   
   printf("r\n");
   printGF5(r);
   printBinary(r);
   printf("\n");
   

   printf("d\n");
   printGF5(d);
   printBinary(d);
   printf("\n");
}//wrong 6 > 6, 2 > 7, 1 > 6

