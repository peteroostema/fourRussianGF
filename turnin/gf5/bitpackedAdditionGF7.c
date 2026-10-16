/* bitpackedAdditionGF7.c was my first attempt for addition and modulus
 * for bitpacked GF7 elements. It doesn't work and was never fixed
 * since a bitpacked four Russian matrix multiply was never made.
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
   long unsigned int a, b, r, n;
   
   printf("a\n");
   for (int i = 0; i < 16; i++){
      long unsigned int tmp = rand() % 7;
      printf("   %lu", tmp);
      a <<= 4;
      a |= tmp;
   }
   printf("\n");
   printBinary(a);
   printf("\n");
   
   printf("b\n");
   for (int i = 0; i < 16; i++){
      long unsigned int tmp = rand() % 7;
      printf("   %lu", tmp);
      b <<= 4;
      b |= tmp;
   }
   printf("\n");
   printBinary(b);
   printf("\n");
   //   a = rand() % 5;
   //   printf("a: %lu\n", a);
   //   b = rand() % 5;
   //   printf("b: %lu\n", b);
   
      int j = rand() % 100;
      printf("iteration: %d\n", j);
   int k = rand() % 10;
   
   unsigned long int s, s0, r0, r1, r2, rMod7Fix, rInt;
   
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < 1000; i++){
      a+=i;
   
      s = a + b;
      s0 = (s & 0x8888888888888888) >> 3;
      s &= 0x7777777777777777;
      rInt = s + s0;
   //   printf("r\n");
   //   printGF5(r);
   //   printBinary(r);
//      printf("\n");
//      printf(" ");
      r0 = rInt & 0x1111111111111111;
   //   printf("r0\n");
   //   printGF5(r0);
   //   printBinary(r0);
   //   printf("\n");
      r1 = (rInt & 0x2222222222222222) >> 1;
   //   printf("r1\n");
   //   printGF5(r1);
   //   printBinary(r1);
   //   printf("\n");
      r2 = (rInt & 0x4444444444444444) >> 2;
   //   printf("r2\n");
   //   printGF5(r2);
   //   printBinary(r2);
   //   printf("\n");
      rMod7Fix = r0 & r1 & r2;
      rMod7Fix |= rMod7Fix << 1; rMod7Fix |= rMod7Fix << 1;
   //   printf("r7Fix\n");
   //   printGF5(rMod7Fix);
   //   printBinary(rMod7Fix);
   //   printf("\n");
      r += (rInt & ~rMod7Fix) + i;
      r <<= 1;
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
   printf("r\n");
   printGF5(r);
   printBinary(r);
   printf("\n");
   startTime = rdtsc();
   for (int i = 0; i < j; i++){
      printf(" ");
   }
   endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
}//wrong 6 > 6, 2 > 7, 1 > 6

