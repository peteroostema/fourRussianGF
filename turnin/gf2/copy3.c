/* copy3.c changes A to be stored in column major order
 * A stores 32 elements from one row in an int then works
 * down the rows 32 columns at a time
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

static inline void createTable(unsigned int* T, unsigned int* Bt, int TRows, int startRow, int n, int t){
   for (int j = 0; j < n; j++){
      T[j] = 0;
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
            if (k == t - 1){
               T[i*n + j] = T[(i - 1)*n + j] ^ Bt[(startRow + k)*n + j];
            } else{
               T[i*n + j] ^= Bt[(startRow + k)*n + j];
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
   int m = 1024;  // rows in A
   int lb = 1024; // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = 1024;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
   //   printf("m: %d\n", m);
   //   printf("lb: %d\n", lb);
   //   printf("la: %d\n", la);
   //   printf("n: %d\n", n);
   //   printf("ra: %d\n", ra);
   //   printf("rb: %d\n", rb);
   int t = 8;  // four Russian divisor
   unsigned int* A = malloc(m*la*(sizeof(int)));
   
   int lilExtra = (t - ((lb*n) % t)) % t;
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
   //      printf("i: %d\n", i);
   //      printf("2^t: %d\n", TRows);
   unsigned int* T = malloc(TRows*n*sizeof(int));
   int columnCounter = -1;
   clock_t start = clock(), diff;
   for (int i = 0; i < tLoops; i++){
      createTable(T, B, TRows, i*t, n, t);
//      printf("T\n");
//      printBinaryMatrix(T, TRows, n);
      if (!(i*t % 32)){
         columnCounter++;
      }
      for (int j = 0; j < m; j++){
//         printf("bitToRead: %d\n", (i*t % 32) + j*32 + columnCounter*32*m);
//         printf("i: %d, t: %d, column: %d\n", i, t, columnCounter);
         unsigned int id = readBits(A, (i*t % 32) + j*32 + columnCounter*32*m, t);
         
//         printf("id: %d\n", id);
         ////id = (TRows - id) % TRows;
         ////printf("id: %d\n", id);
         for (int k = 0; k < n; k++){
            C[j*n + k] ^= T[id*n + k];
         }
      }
//      printf("C\n");
//      printBinaryMatrix(C, m, n);
   }
   diff = clock() - start;
   int milliSeconds = diff * 1000 / CLOCKS_PER_SEC;
   printf("Time %d ms\n", milliSeconds);
   
//   printf("C\n");
//   printBinaryMatrix(C, m, n);
   
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
   
//   printf("Cf\n");
//   int rowCounter = 0;
//   for (int i = 0; i < m; i++){
//      for (int j = 0; j < nActual; j++){
//         printf("%d", Cf[i*nActual + j]);
//      }
//      printf("\t");
//      printf(" %d", rowCounter);
//      printf("\n");
//      rowCounter++;
//   }
//   printf("\n");
   
   printf("Correct = 1, Wrong = 0: %d\n", checkCorrectResult(C, Cf, m, nActual, n));
   //*/
   
   free(A);
   free(B);
   free(C);
   free(T);
   //   free(greyTable);
}
