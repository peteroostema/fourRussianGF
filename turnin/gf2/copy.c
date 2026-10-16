/* copy.c has the first implementation of the method of four
 * Russians. It seems to be limited to multiples of 32 for sizes
 * of n, but accepts any value for t which was later taken out
 * for simplicity.
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

int findNextGreyRow (int position, int depth){
   printf("position: %d\t", position);
   if (position % 2){
      return depth;
   } else {
      return findNextGreyRow(position/2, depth + 1);
   }
   return -1; //error
}

void createTable(unsigned int* T, unsigned int* Bt, int TRows, int startRow, int n, int t){
   for (int j = 0; j < n; j++){
      T[j] = 0;
   }
   for (int i = 1; i < TRows; i++){
      unsigned int rowsToXOR = 0;
      unsigned int position = i;
      while (rowsToXOR <= t){
         //printf("position: %d\t", position);
         if (position % 2){
            break;
         }
         position /= 2;
         rowsToXOR++;
      }
      //printf("rowtoXor: %d\n", rowsToXOR);
      for (int j = 0; j < n; j++){
         for (int k = 0; k <= rowsToXOR; k++){
            if (k == 0){
               T[i*n + j] = T[(i - 1)*n + j] ^ Bt[(startRow + k)*n + j];
            } else{
               T[i*n + j] ^= Bt[(startRow + k)*n + j];
            }
            
         }
      }
      //T[i] = T[i - 1] ^ Bt[
   }
}

unsigned int readBits(unsigned int* A, int start, int length){
   unsigned int matrixStartID = start / 32;
   printf("matrixID: %d\n", matrixStartID);
   unsigned int startIntPos = start % 32;
   printf("intID: %d\n", startIntPos);
   printf("length: %d\n", length);
   if (startIntPos + length > 32){
      unsigned int leftSide = A[matrixStartID] << startIntPos;
      printf("A[0]: %d\n", A[0]);
      printf("leftSide: %d\n", leftSide);
      leftSide >>= (32 - length);
      printf("leftSide: %d\n", leftSide);
      unsigned int rightSide = A[matrixStartID + 1] >> (64 - startIntPos - length);
      printf("rightSide: %d\n", rightSide);
      return leftSide | rightSide;
   } else{
      unsigned result = A[matrixStartID] << startIntPos;
      return result >> (32 - length);
   }
}

int main(){
   int m = 3;  // rows in A
   int lb = 4; // rows in B
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
   for (int i = 1; i < t; i++){
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
   
   for (int i = 0; i < la/t + 1; i++){
      int TRows = power(2, t);
      printf("2^t: %d\n", TRows);
      unsigned int* T = malloc(TRows*n*sizeof(int));
      createTable(T, B, TRows, i*t, n, t);
      printBinaryMatrix(T, TRows, n);
      for (int j = 0; j < m; j++){
         unsigned int id = readBits(A, i*t + 32*j*la, t);
         printf("id: %d\n", id);
         for (int k = 0; k < n; k++){
            C[j*n + k] += T[id*n + k];
         }
      }
      free(T);
   }
   printf("C\n");
   printBinaryMatrix(C, m, n);
   
   free(A);
   free(B);
   free(C);
   free(greyTable);
}
