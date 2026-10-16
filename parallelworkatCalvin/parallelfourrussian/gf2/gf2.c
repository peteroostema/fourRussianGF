/* gf2.c implements a finite field matrix multiply for gf2 in parallel
 *   It uses the meathod of four russians for a O(lg(t)) speed up.
 *   Computation is parallelized over the outer most loop that
 *   complely rewrites the matrix C.
 *   The pragma can be found on line 402
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 * Revised: December 2017 for Calvin College, CS 374 honors
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>
#include <malloc.h>

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

static inline void createTable(unsigned int* T, unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){
      for (int j = 0; j < n; j++){
         T[i*TRows*n + j] = 0;
      }
   }
   //   printf("nq256: %d\n", nq256);
   //   printf("TFirstRow0:\n");
   //   printBinaryMatrix(T, TRows, n);
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
      int leftoverStart = nq256*8;
      //            printf("rowtoXor: %d\n", rowsToXOR);
      for (int k = t - 1; k >= rowsToXOR; k--){
         for (int j = 0; j < nq256; j++){
            for (int l = 0; l < numOfTables; l++){
               if (k == t - 1){
                  _mm256_storeu_si256(T + i*n + j*8 + l*TRows*n, _mm256_xor_si256 (_mm256_load_si256 (T + (i - 1)*n + j*8 + l*TRows*n), _mm256_load_si256 (Bt + (l*t + k)*n + j*8 + tableNumber*32*n)));
                  //                  T[i*n + j + l*TRows*n] = T[(i - 1)*n + j + l*TRows*n] ^ Bt[(l*t + k)*n + j + tableNumber*32*n];
               } else{
                  //                  T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
                  _mm256_storeu_si256(T + i*n + j*8 + l*TRows*n, _mm256_xor_si256 (_mm256_load_si256 (T + i*n + j*8 + l*TRows*n), _mm256_load_si256 (Bt + (l*t + k)*n + j*8 + tableNumber*32*n)));
                  //                  __m256i idk = _mm256_load_si256 (Bt + (l*t + k)*n + j + tableNumber*32*n);
                  //                  unsigned int* idk2 = (unsigned int*) &(idk);
                  //                  printf("BRead: %d\n", (l*t + k)*n + j + tableNumber*32*n);
                  //                  printBinaryMatrix(idk2, 8, 1);
                  //                  __m256i tstuff = _mm256_load_si256 (T + i*n + j*8 + l*TRows*n);
                  //                  unsigned int* tint = (unsigned int*) &(tstuff);
                  //                  printf("TRead: %d\n", i*n + j*8 + l*TRows*n);
                  //                  printBinaryMatrix(tint, 8, 1);
               }
            }
            //            T[i*n + j + l*TRows*n] = T[(i - 1)*n + j + l*TRows*n] ^ Bt[(l*t + k)*n + j + tableNumber*32*n];
            //            T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
            //            T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
            //            T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
         }
         for (int j = 0; j < nr256; j++){
            for (int l = 0; l < numOfTables; l++){
               if (k == t - 1){
                  T[i*n + j + l*TRows*n + leftoverStart] = T[(i - 1)*n + j + l*TRows*n + leftoverStart] ^ Bt[(l*t + k)*n + j + tableNumber*32*n + leftoverStart];
               } else{
                  T[i*n + j + l*TRows*n + leftoverStart] ^= Bt[(l*t + k)*n + j + tableNumber*32*n + leftoverStart];// consider j in B since it goes one by one row
               }
            }
         }
      }
   }
   //   printf("TFirstRow?:\n");
   //   printBinaryMatrix(T, TRows, n);
}

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

void convertBinaryMatrixToFullSize(unsigned int* fullSizeMatrix, unsigned int* binaryMatrix, int rows, int columns, int binaryColumns){
   int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < binaryColumns; j++){
         tmp = binaryMatrix[i*binaryColumns + j];
         if (j < binaryColumns - 1){
            for (int k = 0; k < 32; k++){
               if (tmp & 0x80000000){
                  fullSizeMatrix[i*columns + j * 32 + k] = 1;
               } else {
                  fullSizeMatrix[i*columns + j * 32 + k] = 0;
               }
               tmp <<= 1;
            }
         }
         else {
            for (int k = 0; k < (columns - 1) % 32 + 1; k++){
               if (tmp & 0x80000000){
                  fullSizeMatrix[i*columns + j * 32 + k] = 1;
               } else {
                  fullSizeMatrix[i*columns + j * 32 + k] = 0;
               }
               tmp <<= 1;
            }
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




void printBinaryColumnMatrix(unsigned int* matrix, int rows, int columns){
   int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         tmp = matrix[i + j*rows];
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

int power(int base, int exp){
   int result = 1;
   for (int i = 0; i < exp; i++){
      result *= base;
   }
   return result;
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

int findNextGreyRow (int position, int t){
   for (int i = 1; i <= t; i++) {
      printf("position: %d\t", position);
      if (position % 2){
         return i;
      }
      position /= 2;
   }
   
   return -1; //error
}

unsigned int readBits(unsigned int* A, int start, int length){
   unsigned int matrixStartID = start / 32;
   unsigned int startIntPos = start % 32;
   unsigned result = A[matrixStartID] << startIntPos;
   return result >> (32 - length);
}


unsigned int readBit(unsigned int* A, int start){
   
   unsigned int matrixStartID = start / 32;
   //      printf("matrixID: %d\n", matrixStartID);
   //matrixEndID = end / 32;
   unsigned int startIntPos = start % 32;
   //      printf("intID: %d\n", startIntPos);
   //   printf("length: %d\n", length);
   //endPosInt = end % 32;
   
   ////printf("there");
   unsigned int result = A[matrixStartID] << startIntPos;
   //   printf("result1: %d\n", result);
   result >>= 31;
   //   printf("result2: %d\n", result);
   //   printf("result3: %d\n", result << (32 - startIntPos));
   return result << (32 - startIntPos);
}

int main(int argc, char *argv[]){
   srand(time(NULL));
   
   int m = atoi(argv[1]);  // rows in A
   int lb = atoi(argv[2]); // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = atoi(argv[3]);
   int nThreads = atoi(argv[4]);
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
   printf("m: %d\n", m);
//   printf("lb: %d\n", lb);
//   printf("la: %d\n", la);
//   printf("n: %d\n", n);
//   printf("ra: %d\n", ra);
//   printf("rb: %d\n", rb);
   int t = 8;  // four Russian divisor
   printf("t: %d\n", t);
   unsigned int* A = memalign(64, m*la*(sizeof(int)));
   
   int lilExtra = (32 - (lb % 32)) % 32;
   unsigned int* B = memalign(64, (lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
//   printf("SizeB: %d\n", (lb + lilExtra)*n);
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
//   printf("A\n");
//   printBinaryColumnMatrix(A, m, la);
   
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

   
   int TRows = 1 << t;
   int TSize = TRows * n;
   int numberOfTables = 32 / t;
   unsigned int* T = memalign(64, TSize*(numberOfTables)*sizeof(int));
   int k = 256;
   int readsPerColumn = m*32 / k;
   int rowsLeftover = (m*32 % k) / 32;
   int nq256 = nActual / 256;
   int nr256 = ((n - 1) % 8) + 1;
   if (!(nActual % 256)){
      nr256 = 0;
   }

   printf("256 bit reads: %d\n", readsPerColumn);
   printf("leftOver: %d\n", rowsLeftover);
   printf("B is %d * 256bits + %d ints\n", nq256, nr256);
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i, nq256, nr256);
      //print3DBinaryMatrix(T, TRows, n, numberOfTables);
      //printf("loop in la: %d\n", i);
#pragma omp parallel for num_threads(nThreads)
      for (int j = 0; j < readsPerColumn; j++){
         __m256i id1, id2, id3, id4;
         id4 = _mm256_load_si256(A + j*8 + i*m);
         id1 = _mm256_srli_epi32 (id4, 24);
         register unsigned int* id1Int = (unsigned int*) &(id1);
         id4 = _mm256_slli_epi32(id4, 8);
         id2 = _mm256_srli_epi32 (id4, 24);
         register unsigned int* id2Int = (unsigned int*) &(id2);
         id4 = _mm256_slli_epi32(id4, 8);
         id3 = _mm256_srli_epi32 (id4, 24);
         register unsigned int* id3Int = (unsigned int*) &(id3);
         id4 = _mm256_slli_epi32(id4, 8);
         id4 = _mm256_srli_epi32 (id4, 24);
         register unsigned int* id4Int = (unsigned int*) &(id4);

         //register __m256i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5,
         //                 xmm6, xmm7;
         register __m256i xmm0 = _mm256_setzero_si256();
         register __m256i xmm1 = _mm256_setzero_si256();
         register __m256i xmm2 = _mm256_setzero_si256();
         register __m256i xmm3 = _mm256_setzero_si256();
         register __m256i xmm4 = _mm256_setzero_si256();
         register __m256i xmm5 = _mm256_setzero_si256();
         register __m256i xmm6 = _mm256_setzero_si256();
         register __m256i xmm7 = _mm256_setzero_si256();
         
         register int* subDiv1 = T;
         register int* subDiv2 = T + TSize;
         register int* subDiv3 = T + 2*TSize;
         register int* subDiv4 = T + 3*TSize;
         
         for (int k = 0; k < nq256; k++){ //nActual / 256
            //printf("loop in nq256: %d\n", k);
            xmm0 = _mm256_load_si256 (subDiv1 + (*(id1Int + 0))*n);
            xmm1 = _mm256_load_si256 (subDiv1 + (*(id1Int + 1))*n);
            xmm2 = _mm256_load_si256 (subDiv1 + (*(id1Int + 2))*n);
            xmm3 = _mm256_load_si256 (subDiv1 + (*(id1Int + 3))*n);
            xmm4 = _mm256_load_si256 (subDiv1 + (*(id1Int + 4))*n);
            xmm5 = _mm256_load_si256 (subDiv1 + (*(id1Int + 5))*n);
            xmm6 = _mm256_load_si256 (subDiv1 + (*(id1Int + 6))*n);
            xmm7 = _mm256_load_si256 (subDiv1 + (*(id1Int + 7))*n);
            //            printf("C\n");
            //            printBinaryMatrix(C, m, n);
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv2 + (*(id2Int + 0))*n));
            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv2 + (*(id2Int + 1))*n));
            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv2 + (*(id2Int + 2))*n));
            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv2 + (*(id2Int + 3))*n));
            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv2 + (*(id2Int + 4))*n));
            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv2 + (*(id2Int + 5))*n));
            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv2 + (*(id2Int + 6))*n));
            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv2 + (*(id2Int + 7))*n));
            
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv3 + (*(id3Int + 0))*n));
            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv3 + (*(id3Int + 1))*n));
            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv3 + (*(id3Int + 2))*n));
            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv3 + (*(id3Int + 3))*n));
            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv3 + (*(id3Int + 4))*n));
            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv3 + (*(id3Int + 5))*n));
            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv3 + (*(id3Int + 6))*n));
            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv3 + (*(id3Int + 7))*n));
            
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv4 + (*(id4Int + 0))*n));
            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv4 + (*(id4Int + 1))*n));
            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv4 + (*(id4Int + 2))*n));
            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv4 + (*(id4Int + 3))*n));
            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv4 + (*(id4Int + 4))*n));
            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv4 + (*(id4Int + 5))*n));
            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv4 + (*(id4Int + 6))*n));
            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv4 + (*(id4Int + 7))*n));
            
            //try column major C and or subDiv like reference variables
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + 8*n*j + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (xmm1, _mm256_load_si256 (C + 8*n*j + n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (xmm2, _mm256_load_si256 (C + 8*n*j + 2*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (xmm3, _mm256_load_si256 (C + 8*n*j + 3*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (xmm4, _mm256_load_si256 (C + 8*n*j + 4*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (xmm5, _mm256_load_si256 (C + 8*n*j + 5*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (xmm6, _mm256_load_si256 (C + 8*n*j + 6*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (xmm7, _mm256_load_si256 (C + 8*n*j + 7*n + 8*k)));
            
            
            subDiv1 += 8;
            subDiv2 += 8;
            subDiv3 += 8;
            subDiv4 += 8;
         }

         for (int l = 0; l < nr256; l++){
            C[8*n*j + 8*nq256 + l] ^= T[id1Int[0]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + n + l] ^= T[id1Int[1]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 2*n + l] ^= T[id1Int[2]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 3*n + l] ^= T[id1Int[3]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 4*n + l] ^= T[id1Int[4]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 5*n + l] ^= T[id1Int[5]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 6*n + l] ^= T[id1Int[6]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + 7*n + l] ^= T[id1Int[7]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + l] ^= T[id2Int[0]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + n + l] ^= T[id2Int[1]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 2*n + l] ^= T[id2Int[2]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 3*n + l] ^= T[id2Int[3]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 4*n + l] ^= T[id2Int[4]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 5*n + l] ^= T[id2Int[5]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 6*n + l] ^= T[id2Int[6]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + 7*n + l] ^= T[id2Int[7]*n + l + 8*nq256 + TSize];
            C[8*n*j + 8*nq256 + l] ^= T[id3Int[0]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + n + l] ^= T[id3Int[1]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 2*n + l] ^= T[id3Int[2]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 3*n + l] ^= T[id3Int[3]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 4*n + l] ^= T[id3Int[4]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 5*n + l] ^= T[id3Int[5]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 6*n + l] ^= T[id3Int[6]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + 7*n + l] ^= T[id3Int[7]*n + l + 8*nq256 + 2*TSize];
            C[8*n*j + 8*nq256 + l] ^= T[id4Int[0]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + n + l] ^= T[id4Int[1]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 2*n + l] ^= T[id4Int[2]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 3*n + l] ^= T[id4Int[3]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 4*n + l] ^= T[id4Int[4]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 5*n + l] ^= T[id4Int[5]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 6*n + l] ^= T[id4Int[6]*n + l + 8*nq256 + 3*TSize];
            C[8*n*j + 8*nq256 + 7*n + l] ^= T[id4Int[7]*n + l + 8*nq256 + 3*TSize];
            
//            T[id2*n + TSize + l];
//            C[(j+readsPerColumn)*8*n + l] ^= T[id3*n + TSize*2 + l];
//            C[(j+readsPerColumn)*8*n + l] ^= T[id4*n + TSize*3 + l];
         }
         //printf("C\n");
         //printBinaryMatrix(C, m, n);
      }
      for (int j = 0; j < rowsLeftover; j++){
         unsigned int id1 = readBits(A, j*32 + i*32*m + 256*readsPerColumn, t); //  + i*32*m//
         unsigned int id2 = readBits(A, t + j*32 + i*32*m + 256*readsPerColumn, t);
         unsigned int id3 = readBits(A, 2*t + j*32 + i*32*m + 256*readsPerColumn, t);
         unsigned int id4 = readBits(A, 3*t + j*32 + i*32*m + 256*readsPerColumn, t);
         //               printf("APostition: %d\n", k*t + j*32 + i*32*m);
         //               printf("id: %d, Table#: %d\n", id, k);
         //               printf("CPosition: %d\n", j*n + i);
         //               printf("TValue: ");
         //               printBinary(T[id + TSize*k]);
         //               printf("\n");
         for (int l = 0; l < n; l++){
            C[(j+8*readsPerColumn)*n + l] ^= T[id1*n + l];
            C[(j+8*readsPerColumn)*n + l] ^= T[id2*n + TSize + l];
            C[(j+8*readsPerColumn)*n + l] ^= T[id3*n + TSize*2 + l];
            C[(j+8*readsPerColumn)*n + l] ^= T[id4*n + TSize*3 + l];
         }
      }
      //      printf("\n");
      
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu s\n", (endTime - startTime));
   
   //   printf("C\n");
   //   printBinaryMatrix(C, m, n);

   
   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   /*
//   startTime = rdtsc();
   unsigned int* Af = memalign(64, m*lb*sizeof(int));
   unsigned int* Bf = memalign(64, lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   //   printf("Af\n");
   convertBinaryMatrixToFullSizeColumnMajor(Af, A, m, lb, la);
   //   printf("Bf\n");
   convertBinaryMatrixToFullSize(Bf, B, lb, nActual, n);
   
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            //C[i * n + j] ^= (readBit(A, i * la * 32 + k) & readBit(B, j + k * n));
            Cf[i * nActual + j] ^= (Af[i * lb + k] & Bf[j + k * nActual]);
         }
      }
   }
//   endTime = rdtsc();
//   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   
//            printf("Cf\n");
//            int rowCounter = 0;
//            for (int i = 0; i < m; i++){
//               for (int j = 0; j < nActual; j++){
//                  printf("%d", Cf[i*nActual + j]);
//               }
//               printf("\t");
//               printf(" %d", rowCounter);
//               printf("\n");
//               rowCounter++;
//            }
//            printf("\n");
   
   printf("Correct = 1, Wrong = 0: %d\n", checkCorrectResult(C, Cf, m, nActual, n));
   //*/
   //   printf("B\n");
   //   printBinaryMatrix(B, lb + lilExtra, n);
   //   printf("BfirstRow: \n");
   //   printBinaryMatrix(B, 8, 1);
   //   printf("A\n");
   //   printBinaryColumnMatrix(A, m, la);
   //
   //   __m256i xorTest = _mm256_xor_si256 (_mm256_load_si256 (A), _mm256_load_si256 (B));
   //   unsigned int* xortest = (unsigned int*) &xorTest;
   //   printBinaryMatrix(xortest, 8, 1);
}
