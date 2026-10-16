/* bitpackedAdditionGF7.c was my first attempt for addition and modulus
 * for bitpacked GF5 elements. It doesn't work but was fixed in full
 * matrix multiply programs.
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
   printf("\n");
   printBinary(b);
   printf("\n");
   //   a = rand() % 5;
   //   printf("a: %lu\n", a);
   //   b = rand() % 5;
   //   printf("b: %lu\n", b);
//   int j = rand() % 1000;
//   printf("iteration: %d", j);
   
   unsigned long int s, s0, s1, s2, rMod5Fix;
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < 1000; i++){
      a++;
   
   s = a + b;
   
   
//   rFix8 = (rFix8 >> 1) | (rFix8 >> 3);
   
   s0 = s & 0x1111111111111111;
   s1 = (s & 0x2222222222222222) >> 1;
   s2 = (s & 0x4444444444444444) >> 2;
//   printf("s0\n");
//   printGF5(s0);
//   printBinary(s0);
//   printf("\n");
//   printf("s1\n");
//   printGF5(s1);
//   printBinary(s1);
//   printf("\n");
//   printf("s2\n");
//   printGF5(s2);
//   printBinary(s2);
//   printf("\n");
   
   rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
   rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//   
//   fFix7 = s0 & s1 & s2;
//   fFix6 = s2 & s1;
//   fFix5 = s2 & s0;

//      printf("r5Fix\n");
//      printGF5(rMod5Fix);
//      printBinary(rMod5Fix);
//      printf("\n");
   r ^= (s - rMod5Fix);
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
   printf("r\n");
   printGF5(r);
   printBinary(r);
   printf("\n");
   
   
}//wrong 6 > 6, 2 > 7, 1 > 6

