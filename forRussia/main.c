/* main.c immplements an early version of createTable, but doesn't use to to do a matrix multiplication
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printBinaryMatrix(unsigned int* matrix, int rows, int columns){
   int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         tmp = matrix[i*columns + j];
         for (int k = 0; k < 32; k++){
            if (tmp & 0x80000000){
               printf("1");
            } else {
               printf("0");
            }
            tmp <<= 1;
         }
         printf("\t");
      }
      printf(" %d", rowCounter);
      rowCounter++;
      printf("\n");
   }
   printf("\n");
}

int power(int base, int exp){
   int result = 1;
   for (int i = 0; i < exp; i++){
      result *= base;
   }
   return result;
}

void createTable(unsigned int* matrix, int size, int startCol, int t, unsigned int* T){
   unsigned int currentPower = 2;
   
   T[0] = 0;
   T[1] = 1;
   for (int i = 0; i < t; i++){
      for (int j = 0; j < currentPower; j++){
         T[currentPower + j] = T[currentPower - j];
         T[currentPower + j] |= currentPower;
      }
      currentPower *= 2;
   }
}

int main(){
   int m = 2;  // rows in A
   int lb = 6; // rows in B
   int la = lb / 32 + 1;  // columns in A
   int n = (64 - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (64 % 32)) % 32; // zeros to pad on B
   printf("m: %d\n", m);
   printf("lb: %d\n", lb);
   printf("la: %d\n", la);
   printf("n: %d\n", n);
   printf("ra: %d\n", ra);
   printf("rb: %d\n", rb);
   int t = 4;  // four Russian divisor
   unsigned int* A = malloc(m*la*(sizeof(int)));
   unsigned int* B = malloc(lb*n*(sizeof(int)));
   unsigned int* C = calloc(m*n,(sizeof(int)));
   int b = 0; // counter for printing
   srand(time(NULL));
   
   unsigned int* greyTable = malloc(power(2, t)*sizeof(int));
   unsigned int currentPower = 2;
   
   greyTable[0] = 0;
   greyTable[1] = 1;
   for (int i = 0; i < t; i++){
      for (int j = 0; j < currentPower; j++){
         greyTable[currentPower + j] = greyTable[currentPower - j - 1];
         greyTable[currentPower + j] |= currentPower;
      }
      currentPower *= 2;
   }
   printf("greyTable\n");
   printBinaryMatrix(greyTable, power(2, t), 1);
   
   int tmp =  rand();
   printf("tmp: %d\n", tmp);
   tmp >>= ra;
   printf("tmp: %d\n", tmp);
   tmp <<= ra;
   printf("tmp: %d\n", tmp);

   for (int i = 0; i < m; i++){// init A
      for (int j = 0; j < la; j++){
         if (j < la - 1){
            A[i*la + j] = rand();
         } else {
            A[i*la + j] = rand() >> ra;
            A[i*la + j] <<= ra;
         }
      }
   }
   printf("A\n");
   printBinaryMatrix(A, m, la);
   
   for (int i = 0; i < lb; i++){// init B
      for (int j = 0; j < n; j++){
         if (j < n - 1){
            B[i*n + j] = rand();
         } else {
            B[i*n + j] = rand() >> rb;
            B[i*n + j] <<= rb;
         }
      }
   }
   printf("B\n");
   printBinaryMatrix(B, lb, n);
   
//   for (int i = 0; i < m; i++){ // row
//      for (int j = 0; j < n; j++){ // column
//         for (int k = 0; k < la; k++){ // interate
//            C[i * n + j] = (A[i * la + k] & B[j + k * n]) ^ C[i * n + j];
//         }
//      }
//   }
   
   for (int i = 0; i < la; i++){
      int size = power(2, t);
      printf("2^t: %d\n", size);
      unsigned int* T = malloc(size*sizeof(int));
      createTable(B, size, i*t, t, T);
      free(T);
   }
   printf("C\n");
   printBinaryMatrix(C, m, n);
   
   free(A);
   free(B);
   free(C);
   free(greyTable);
}
