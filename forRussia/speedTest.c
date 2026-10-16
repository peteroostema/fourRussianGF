/* speedTest.c is an early immplementation I run speed tests on
 * Author: Peter Oostema
 * Date: Summer 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

int power(int base, int exp){
   int result = 1;
   for (int i = 0; i < exp; i++){
      result *= base;
   }
   return result;
}

static inline void createTable(unsigned int* T, unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber){
   for (int j = 0; j < n; j++){
      for (int l = 0; l < numOfTables; l++){
         T[j + TRows*l*n] = 0;
      }
   }
   for (int i = 1; i < TRows; i++){
      int rowsToXOR = t - 1;
      unsigned int position = i;
      while (rowsToXOR >= 0){
         //         printf("position: %d\t", position);
         if (position % 2){
            break;
         }
         position /= 2;
         rowsToXOR--;
      }
      //      printf("rowtoXor: %d\n", rowsToXOR);
      for (int j = 0; j < n; j++){
         for (int k = t - 1; k >= rowsToXOR; k--){
            for (int l = 0; l < numOfTables; l++){
               if (k == t - 1){
                  T[i*n + j + l*TRows*n] = T[(i - 1)*n + j + l*TRows*n] ^ Bt[(l*t + k)*n + j + tableNumber*32*n];
               } else{
                  T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
               }
            }
         }
      }
   }
}

static inline unsigned int readBits(unsigned int* A, int start, int length){
   unsigned int matrixStartID = start / 32;
   unsigned int startIntPos = start % 32;
   unsigned int result = A[matrixStartID] << startIntPos;
   return result >> (32 - length);
}

int main(){
   int m = 5000;  // rows in A
   int lb = 5000; // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = 4090;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
   int t = 8;  // four Russian divisor
   unsigned int* A = malloc(m*la*(sizeof(int)));
   
   int lilExtra = (32 - (lb % 32)) % 32;
   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
   unsigned int* C = calloc(m*n,(sizeof(int)));
   srand(time(NULL));
   
   for (int j = 0; j < la; j++){// init A column major order
      for (int i = 0; i < m; i++){
         if (j < la - 1){
            A[i + j*m] = rand();
         } else {
            A[i + j*m] = rand() >> ra;
            A[i + j*m] <<= ra;
         }
      }
   }
   
   for (int i = 0; i < lb; i++){// init B
      for (int j = 0; j < n; j++){
         if (j < n - 1){
            B[i*n + j] = rand();
         } else {
            B[i*n + j] = rand() >> rb;
            B[i*n + j] <<= rb;
         }
      }
      for (int i = 0; i < lilExtra; i++){
         for (int j = 0; j < n; j++){
            B[(i+lb)*n + j] = 0;
         }
      }
   }
   
   int TRows = power(2, t);
   int TSize = TRows * n;
   int numberOfTables = 32 / t;
   unsigned int* T = malloc(TSize*(numberOfTables)*sizeof(int));
   
   unsigned long long startTime = rdtsc();
   
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
         for (int j = 0; j < m; j++){
               unsigned int id1 = readBits(A, j*32 + i*32*m, t); //  + i*32*m//
               unsigned int id2 = readBits(A, t + j*32 + i*32*m, t);
               unsigned int id3 = readBits(A, 2*t + j*32 + i*32*m, t);
               unsigned int id4 = readBits(A, 3*t + j*32 + i*32*m, t);
               for (int l = 0; l < n; l++){
                  C[j*n + l] ^= T[id1*n + l];
                  C[j*n + l] ^= T[id2*n + TSize + l];
                  C[j*n + l] ^= T[id3*n + TSize*2 + l];
                  C[j*n + l] ^= T[id4*n + TSize*3 + l];
               }
         }
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   
   
   free(A);
   free(B);
   free(C);
   free(T);
}
