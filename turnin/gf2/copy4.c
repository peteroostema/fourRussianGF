/* copy4.c adds subTables to the four Russian lookup table.
 * Previously createTables only created tables for t rows at
 * a time. Now it creates tables for 32 rows at a time since
 * that many elemnts of A are read at one. subTables are 
 * created t rows at a time and there are 32/t subTables.
 *
 * The outer most loop is now la (number of integer columns)
 * Perviously it was tloops = lb/t (number of columns over t)
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// convert matrix for testing
void convertBinaryMatrixToFullSize(unsigned int* fullSizeMatrix, unsigned int* binaryMatrix, int rows, int columns, int binaryColumns){
   int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < binaryColumns; j++){
         tmp = binaryMatrix[i*binaryColumns + j];
         for (int k = 0; k < 32; k++){
            if (j*32 + k >= columns){
               break;
            }
            if (tmp & 0x80000000){
               fullSizeMatrix[i*columns + j * 32 + k] = 1;
            } else {
               fullSizeMatrix[i*columns + j * 32 + k] = 0;
            }
            tmp <<= 1;
         }
      }
      //rowCounter++;
   }
   //   printf("size: %d\n", rows*columns);
   //   for (int i = 0; i < rows; i++){
   //      for (int j = 0; j < columns; j++){
   //         printf("%d", fullSizeMatrix[i*columns + j]);
   //      }
   //      printf("\t");
   //      printf(" %d", rowCounter);
   //      printf("\n");
   //      rowCounter++;
   //   }
   //   printf("\n");
}
void convertBinaryMatrixToFullSizeColumnMajor(unsigned int* fullSizeMatrix, unsigned int* binaryMatrix, int rows, int columns, int binaryColumns){
   int tmp = 0;
   int rowCounter = 0;
   for (int j = 0; j < binaryColumns; j++){
      for (int i = 0; i < rows; i++){
         tmp = binaryMatrix[i + j*rows];
         for (int k = 0; k < 32; k++){
            if (j*32 + k >= columns){
               break;
            }
            if (tmp & 0x80000000){
               fullSizeMatrix[i*columns + j * 32 + k] = 1;
            } else {
               fullSizeMatrix[i*columns + j * 32 + k] = 0;
            }
            tmp <<= 1;
         }
      }
      //rowCounter++;
   }
   //   printf("size: %d\n", rows*columns);
   //   for (int i = 0; i < rows; i++){
   //      for (int j = 0; j < columns; j++){
   //         printf("%d", fullSizeMatrix[i*columns + j]);
   //      }
   //      printf("\t");
   //      printf(" %d", rowCounter);
   //      printf("\n");
   //      rowCounter++;
   //   }
   //   printf("\n");
}
int checkCorrectResult(unsigned int* C, unsigned int* Cf, int rows, int columns, int binaryColumns){
   int tmp = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < binaryColumns; j++){
         tmp = C[i*binaryColumns + j];
         for (int k = 0; k < 32; k++){
            if (j*32 + k >= columns){
               break;
            }
            if (tmp & 0x80000000){
               if (!(Cf[i*columns + j * 32 + k] == 1)){
                  return 0;
               }
            } else {
               if (!(Cf[i*columns + j * 32 + k] == 0)){
                  return 0;
               }
            }
            tmp <<= 1;
         }
      }
   }
   return 1;
}
void createTableOld(unsigned int* T, unsigned int* Bt, int TRows, int startRow, int n, int t){
   for (int j = 0; j < n; j++){
      T[j] = 0;
   }
   for (int i = 1; i < TRows; i++){
      int rowsToXOR = t - 1;
      unsigned int position = i;
      while (rowsToXOR >= 0){
         //printf("position: %d\t", position);
         if (position % 2){
            break;
         }
         position /= 2;
         rowsToXOR--;
      }
      //printf("rowtoXor: %d\n", rowsToXOR);
      for (int j = 0; j < n; j++){
         for (int k = t - 1; k >= rowsToXOR; k--){
            if (k == t - 1){
               T[i*n + j] = T[(i - 1)*n + j] ^ Bt[(startRow + k)*n + j];
            } else{
               T[i*n + j] ^= Bt[(startRow + k)*n + j];
            }
            
         }
      }
      //T[i] = T[i - 1] ^ Bt[
   }
}

void printBinary(int number){
   int tmp = number;
   for (int k = 0; k < 32; k++){
      if (tmp & 0x80000000){
         printf("1");
      } else {
         printf("0");
      }
      tmp <<= 1;
   }
   printf("\n");
}

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

void print3DBinaryMatrix(unsigned int* matrix, int rows, int columns, int depth){
   int tmp = 0;
   int rowCounter = 0;
   for (int l = 0; l < depth; l++){
      printf("subTable: %d\n", l);
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < columns; j++){
            tmp = matrix[i*columns + j + l*rows*columns];
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
      rowCounter = 0;
   }
   printf("\n");
}

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

int findNextGreyRow (int position, int depth){
   printf("position: %d\t", position);
   if (position % 2){
      return depth;
   } else {
      return findNextGreyRow(position/2, depth + 1);
   }
   return -1; //error
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
                  //T[i*n + j] = T[(i - 1)*n + j] ^ Bt[(startRow + k)*n + j];
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
   //   printf("matrixID: %d\n", matrixStartID);
   //   printf("intID: %d\n", startIntPos);
   //   printf("length: %d\n", length);
   unsigned int result = A[matrixStartID] << startIntPos;
   return result >> (32 - length);
}

int main(){
   int m = 512;  // rows in A
   int lb = 512; // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = 512;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
//      printf("m: %d\n", m);
//      printf("lb: %d\n", lb);
//      printf("la: %d\n", la);
//      printf("n: %d\n", n);
   //   printf("ra: %d\n", ra);
   //   printf("rb: %d\n", rb);
   int t = 8;  // four Russian divisor
   unsigned int* A = malloc(m*la*(sizeof(int)));
   
   int lilExtra = (32 - (lb % 32)) % 32;
   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
   //   printf("SizeB: %d\n", (lb + lilExtra)*n);
   unsigned int* C = calloc(m*n,(sizeof(int)));
   //int b = 0; // counter for printing
   srand(time(NULL));
   
   //   unsigned int* greyTable = malloc(power(2, t)*sizeof(int));
   //   unsigned int currentPower = 2;
   
   //   greyTable[0] = 0;
   //   greyTable[1] = 1;
   //   for (int i = 1; i < t; i++){
   //      for (int j = 0; j < currentPower; j++){
   //         greyTable[currentPower + j] = greyTable[currentPower - j - 1];
   //         greyTable[currentPower + j] |= currentPower;
   //      }
   //      currentPower *= 2;
   //   }
   //   printf("greyTable\n");
   //   printBinaryMatrix(greyTable, power(2, t), 1);
   
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
//   printf("A\n");
//   printBinaryMatrix(A, m, la);
   
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
//   printf("B\n");
//   printBinaryMatrix(B, lb + lilExtra, n);
   
   //   for (int i = 0; i < m; i++){ // row
   //      for (int j = 0; j < n; j++){ // column
   //         for (int k = 0; k < la; k++){ // interate
   //            C[i * n + j] = (A[i * la + k] & B[j + k * n]) ^ C[i * n + j];
   //         }
   //      }
   //   }
   
   int tLoops = lb / t;
   if (lb % t){
      tLoops++;
   }
   int TRows = power(2, t);
   int TSize = TRows * n;
   int numberOfTables = 32 / t;
   unsigned int* T = malloc(TSize*(numberOfTables)*sizeof(int));
   
   clock_t start = clock(), diff;
   unsigned long long startTime = rdtsc();
   
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
//      print3DBinaryMatrix(T, TRows, n, numberOfTables);
//      for (int l = 0; l < la; l++){
         for (int j = 0; j < m; j++){
            for (int k = 0; k < numberOfTables; k++){
               unsigned int id = readBits(A, k*t + j*32 + i*32*m, t); //  + i*32*m//
//               printf("APostition: %d\n", k*t + j*32 + i*32*m);
//               printf("id: %d, Table#: %d\n", id, k);
//               printf("CPosition: %d\n", j*n + i);
//               printf("TValue: ");
//               printBinary(T[id + TSize*k]);
//               printf("\n");
               for (int l = 0; l < n; l++){
                  C[j*n + l] ^= T[id*n + TSize*k + l];
               }
            }
//            printf("C\n");
//            printBinaryMatrix(C, m, n);
         }
//      }
   }
   
//      int columnCounter = -1;
//      for (int i = 0; i < n; i++){
//         createTableOld(T, B, TRows, TSize, n, t);
//   //      printf("T\n");
//   //      printBinaryMatrix(T, TRows, n);
//   //      if (!(i*t % 32)){
//   //         columnCounter++;
//   //      }
//         for (int j = 0; j < m; j++){
//            unsigned int id = readBits(A, (i*t % 32) + j*32 + columnCounter*32*m, t);
//   //         printf("id: %d\n", id);
//            for (int k = 0; k < n; k++){
//               C[j*n + k] ^= T[id*n + k];
//            }
//         }
//   //      printf("C\n");
//   //      printBinaryMatrix(C, m, n);
//      }
   unsigned long long endTime = rdtsc();
   diff = clock() - start;
   int milliSeconds = diff * 1000 / CLOCKS_PER_SEC;
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("Time %d ms\n", milliSeconds);
   
//      printf("C\n");
//      printBinaryMatrix(C, m, n);
   
   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   ///*
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   //   printf("Af\n");
   convertBinaryMatrixToFullSizeColumnMajor(Af, A, m, lb, la);
   //   printf("Bf\n");
   convertBinaryMatrixToFullSize(Bf, B, lb, nActual, n);
   
   start = clock();
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            //C[i * n + j] ^= (readBit(A, i * la * 32 + k) & readBit(B, j + k * n));
            Cf[i * nActual + j] ^= (Af[i * lb + k] & Bf[j + k * nActual]);
         }
      }
   }
   diff = clock() - start;
   milliSeconds = diff * 1000 / CLOCKS_PER_SEC;
   printf("Time %d ms\n", milliSeconds);
   
//      printf("Cf\n");
//      int rowCounter = 0;
//      for (int i = 0; i < m; i++){
//         for (int j = 0; j < nActual; j++){
//            printf("%d", Cf[i*nActual + j]);
//         }
//         printf("\t");
//         printf(" %d", rowCounter);
//         printf("\n");
//         rowCounter++;
//      }
//      printf("\n");
   
   printf("Correct = 1, Wrong = 0: %d\n", checkCorrectResult(C, Cf, m, nActual, n));
   //*/
   
   free(A);
   free(B);
   free(C);
   free(T);
   //   free(greyTable);
}
