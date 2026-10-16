/* checkGF7.c is a addition and modulus of bit sliced
 * elements of GF7
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

void printGF7(long unsigned int r2, long unsigned int r1, long unsigned int r0){
   for (int i = 63; i >= 0; i--){
      long unsigned int tmp0 = (r0 >> i) & 1;
      long unsigned int tmp1 = (r1 >> i) & 1;
      long unsigned int tmp2 = (r2 >> i) & 1;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("%lu", tmp2*4 + tmp1*2 + tmp0);
   }
   printf("\n");
}
void printGF7Added(long unsigned int r3, long unsigned int r2, long unsigned int r1, long unsigned int r0){
   for (int i = 63; i >= 0; i--){
      long unsigned int tmp0 = (r0 >> i) & 1;
      long unsigned int tmp1 = (r1 >> i) & 1;
      long unsigned int tmp2 = (r2 >> i) & 1;
      long unsigned int tmp3 = (r3 >> i) & 1;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("%lu", tmp3*8 + tmp2*4 + tmp1*2 + tmp0);
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
      long unsigned int tmp = rand() % 7;
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
      long unsigned int tmp = rand() % 7;
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
   printf("\n");
   printBinary(b2);
   printBinary(b1);
   printBinary(b0);
   printf("\n");
   
   unsigned long long startTime = rdtsc();
   s0 = a0 ^ b0;//std bit packed addition
   s1 = a1 ^ b1;
   s2 = a2 ^ b2;
   s3 = (s2 & ((a1 & b1) | (s1 & (a0 & b0)))) ^ (a2 & b2);
   s2 ^= (s1 & (a0 & b0)) ^ a1 & b1;;
   s1 ^= a0 & b0;
//   s2 ^= a1 & b1;
   
//   r0 = s0 ^ s3;
//   r1 = s1;
//   r2 = s2;
//   r2 ^= r1 & (s0 & s3);
//   r1 ^= s0 & s3;
   r2 = (s2 ^ (s1 & (s0 & s3))) & ~(s0 & s1 & s2); //my meathod
   r1 = (s1 ^ (s0 & s3)) & ~(s0 & s1 & s2);
   r0 = (s0 ^ s3) & ~(s0 & s1 & s2);
   
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
   printf("s\n");
   printGF7Added(s3, s2, s1, s0);
   printBinary(s3);
   printBinary(s2);
   printBinary(s1);
   printBinary(s0);
   printf("\n");
   
   printf("r\n");
   printGF7(r2, r1, r0);
   printBinary(r2);
   printBinary(r1);
   printBinary(r0);
   printf("\n");
   
   
}//wrong 6 > 6, 2 > 7, 1 > 6

