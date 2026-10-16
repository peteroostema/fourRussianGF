/* hopefullyTheLast.c was my first attempt at a kernel 
 * version of my program which I was hoping was the last,
 * but wasn't.
 *
 * This implementation hits about 5% of peak (16 cycles).
 * Later attempts use different unrolling and ways of 
 * parsing indexs from A
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>

#define VXOR(srca,srcb,dest)             \
asm volatile(                               \
"vpxor %[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
: [vsrca] "x"(srca),                    \
[vsrcb] "x"(srcb)                     \
);

#define VLOAD(srcaddr,dest)           \
asm volatile(                            \
"vmovdqu %[addr], %[rvdest]\n"       \
: [rvdest] "=x"(dest)                \
: [addr] "m"(*srcaddr)               \
);

void print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, innerColumns, int outerColumns);
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
      _mm256_storeu_si256(T + i*TRows*8, _mm256_setzero_si256());
   }
   for (int i = 1; i < TRows; i++){
      int rowsToXOR = t - 1;
      unsigned int position = i;
      while (rowsToXOR >= 0){
         if (position % 2){
            break;
         }
         position /= 2;
         rowsToXOR--;
      }
      for(int e = 0; e < nq256; e++){
      for (int k = t - 1; k >= rowsToXOR; k--){
         for (int l = 0; l < numOfTables; l++){
            if (k == t - 1){
               _mm256_storeu_si256(T + l*TRows*8 + i*8 + e*numOfTables*TRows*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*8 + (i - 1)*8 + e*numOfTables*TRows*8), _mm256_load_si256 (Bt + (l*t + k)*n + tableNumber*64*n + 8*e)));
            } else{
               _mm256_storeu_si256(T + l*TRows*8 + i*8 + e*numOfTables*TRows*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*8 + i*8 + e*numOfTables*TRows*8), _mm256_load_si256 (Bt + (l*t + k)*n + tableNumber*64*n + e*8)));
            }
         }
      }
      }
   }
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

void print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns){
   int tmp = 0;
   int rowCounter = 0;
   int subtableSize = rows*(majorColumns*innerColumns + outerColumns);
   for (int l = 0; l < depth; l++){
      printf("subTable: %d\n", l);
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < majorColumns; j++){
            for (int m = 0; m < innerColumns; m++){
               tmp = matrix[l*subtableSize + j*rows*innerColumns + i*innerColumns + m];
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
         }
         for (int m = 0; m < outerColumns; m++){
            tmp = matrix[l*subtableSize + majorColumns*rows*innerColumns + m*rows + i];
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
   unsigned int startIntPos = start % 32;
   unsigned int result = A[matrixStartID] << startIntPos;
   result >>= 31;
   return result << (32 - startIntPos);
}

int main(){
   srand(time(NULL));
   
   int m = 4;  // rows in A
   int lb = 64; // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = 256;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
   //   printf("m: %d\n", m);
   //   printf("lb: %d\n", lb);
   //   printf("la: %d\n", la);
   //   printf("n: %d\n", n);
   //   printf("ra: %d\n", ra);
   //   printf("rb: %d\n", rb);
   
   //   printf("t: %d\n", t);
   unsigned int* A = malloc(m*la*(sizeof(int)));
   
   int lilExtra = (32 - (lb % 32)) % 32;
   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
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
   
   int t = 4;  // four Russian divisor
   unsigned int mask = 0xF;
   int TRows = 1 << t;
   int TSize = TRows * 8;
   int numberOfTables = 64 / t;
   
   int k = 256;
   int readsPerColumn = m*32 / k;
   int rowsLeftover = (m*32 % k) / 32;
   int nq256 = nActual / 256;
   int nr256 = ((n - 1) % 8) + 1;
   if (!(nActual % 256)){
      nr256 = 0;
   }
   unsigned int* T = malloc(nq256*TSize*(numberOfTables)*sizeof(int));
   printf("256 bit reads: %d\n", readsPerColumn);
   printf("leftOver: %d\n", rowsLeftover);
   printf("B is %d * 256bits + %d ints\n", nq256, nr256);
   
   createTable(T, B, TRows, n, t, numberOfTables, 0, nq256, nr256);
   //       printBinaryMatrix(T, 8, 1);
   //       printBinaryMatrix(T, TRows, n);
//   print4DBinaryMatrix(T, TRows, 1, numberOfTables*nq256, 8, 0);
   //print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, innerColumns, int outerColumns
   
   __m256i xmm0, xmm1, xmm2, xmm3, ymm0, ymm1, ymm2, ymm3;
   register unsigned int* subDiv = T;
   unsigned int* cPosition = C;
   
//   unsigned int* cIntermediateTest = calloc(m*n,(sizeof(int)));
   
   unsigned long long startTime = rdtsc();
   for (int k = 0; k < nq256; k++){ //nActual / 256
//      createTable(T, B, TRows, n, t, numberOfTables, 0, nq256, nr256);
      for (int j = 0; j < m/4; j++){
         for (int i = 0; i < la; i++){
            xmm0 = _mm256_setzero_si256();
            xmm1 = _mm256_setzero_si256();
            xmm2 = _mm256_setzero_si256();
            xmm3 = _mm256_setzero_si256();
//            unsigned int a0 = *(A + i*m + 0);
//            unsigned int a1 = *(A + i*m + 1);
//            unsigned int a2 = *(A + i*m + 2);
//            unsigned int a3 = *(A + i*m + 3);
            for (int ii = 0; ii < numberOfTables/2; ii++){
               unsigned int index0 = ((*(A + i*m + 0) >> (28 - ((ii) * t))) & mask) << 3; //variable t (32 - ((ii+1) * t)))
               VLOAD((subDiv + index0), ymm0);
               VXOR(xmm0, ymm0, xmm0);
               
               unsigned int index1 = ((*(A + i*m + 1) >> (28 - ((ii) * t))) & mask) << 3;
               VLOAD((subDiv + index1), ymm1);
               VXOR(xmm1, ymm1, xmm1);
               
               unsigned int index2 = ((*(A + i*m + 2) >> (28 - ((ii) * t))) & mask) << 3;
               VLOAD((subDiv + index2), ymm2);
               VXOR(xmm2, ymm2, xmm2);
               
               unsigned int index3 = ((*(A + i*m + 3) >> (28 - ((ii) * t))) & mask) << 3;
               VLOAD((subDiv + index3), ymm3);
               VXOR(xmm3, ymm3, xmm3);
               
               subDiv += TSize;
//                           printf("id0: %u\n", index0/8);
//                           printf("id1: %u\n", index1/8);
//                           printf("id2: %u\n", index2/8);
//                           printf("id3: %u\n", index3/8);
//               
//               _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, xmm0);
//               _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, xmm1);
//               _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, xmm2);
//               _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, xmm3);
//               printf("CIntermediate\n");
//               printBinaryMatrix(cIntermediateTest, m, n);
            }
//            for (int ii = numberOfTables/2; ii < numberOfTables; ii++){
//               unsigned int index0 = ((*(A + m + 0) >> (28 - (ii*4))) & mask) << 2;
//               VLOAD((B + index0), ymm0);
//               VXOR(xmm0, ymm0, xmm0);
//               
//               unsigned int index1 = ((*(A + m + 1) >> (28 - (ii*4))) & mask) << 2;
//               VLOAD((B + index1), ymm1);
//               VXOR(xmm1, ymm1, xmm1);
//               
//               unsigned int index2 = ((*(A + m + 2) >> (28 - (ii*4))) & mask) << 2;
//               VLOAD((B + index2), ymm2);
//               VXOR(xmm2, ymm2, xmm2);
//               
//               unsigned int index3 = ((*(A + m + 3) >> (28 - (ii*4))) & mask) << 2;
//               VLOAD((B + index3), ymm3);
//               VXOR(xmm3, ymm3, xmm3);
//                  printf("id0: %u\n", index0/4);
//                  printf("id1: %u\n", index1/4);
//                  printf("id2: %u\n", index2/4);
//                  printf("id3: %u\n", index3/4);
//            }
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + 8*n*j + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (xmm1, _mm256_load_si256 (C + 8*n*j + n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (xmm2, _mm256_load_si256 (C + 8*n*j + 2*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (xmm3, _mm256_load_si256 (C + 8*n*j + 3*n + 8*k)));
         } // i
         //cPosition += 8*n;
      } // j
      //cPostition += 8;
   } // k
//            id4 = _mm256_load_si256(A + j*8 + i*m);// address A by integer// run on CMU cluster
//            id1 = _mm256_srli_epi32 (id4, 24);
//            id1 = _mm256_slli_epi32(id1, 3);
//            register unsigned int* id1Int = (unsigned int*) &(id1);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id2 = _mm256_srli_epi32 (id4, 24);
//            id2 = _mm256_slli_epi32(id2, 3);
//            register unsigned int* id2Int = (unsigned int*) &(id2);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id3 = _mm256_srli_epi32 (id4, 24);
//            id3 = _mm256_slli_epi32(id3, 3);
//            register unsigned int* id3Int = (unsigned int*) &(id3);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id4 = _mm256_srli_epi32 (id4, 24);
//            id4 = _mm256_slli_epi32(id4, 3);
//            register unsigned int* id4Int = (unsigned int*) &(id4);
//            //         printf("id1: \n");
//            //         printBinary(id1Int[0]); printBinary(id1Int[1]); printBinary(id1Int[2]); printBinary(id1Int[3]); printBinary(id1Int[4]); printBinary(id1Int[5]); printBinary(id1Int[6]); printBinary(id1Int[7]);
//            //         printf("id2: \n");
//            //         printBinary(id2Int[0]); printBinary(id2Int[1]); printBinary(id2Int[2]); printBinary(id2Int[3]); printBinary(id2Int[4]); printBinary(id2Int[5]); printBinary(id2Int[6]); printBinary(id2Int[7]);
//            //         printf("id3: \n");
//            //         printBinary(id3Int[0]); printBinary(id3Int[1]); printBinary(id3Int[2]); printBinary(id3Int[3]); printBinary(id3Int[4]); printBinary(id3Int[5]); printBinary(id3Int[6]); printBinary(id3Int[7]);
//            //         printf("id4: \n");
//            //         printBinary(id4Int[0]); printBinary(id4Int[1]); printBinary(id4Int[2]); printBinary(id4Int[3]); printBinary(id4Int[4]); printBinary(id4Int[5]); printBinary(id4Int[6]); printBinary(id4Int[7]);
//            printf("id0: %u\n", id1Int[0]/8);
//            printf("id1: %u\n", id1Int[1]/8);
//            printf("id2: %u\n", id1Int[2]/8);
//            printf("id3: %u\n", id1Int[3]/8);
//            printf("id4: %u\n", id1Int[4]/8);
//            printf("id5: %u\n", id1Int[5]/8);
//            printf("id6: %u\n", id1Int[6]/8);
//            printf("id7: %u\n", id1Int[7]/8);
//            printf("id0: %u\n", id2Int[0]/8);
//            printf("id1: %u\n", id2Int[1]/8);
//            printf("id2: %u\n", id2Int[2]/8);
//            printf("id3: %u\n", id2Int[3]/8);
//            printf("id4: %u\n", id2Int[4]/8);
//            printf("id5: %u\n", id2Int[5]/8);
//            printf("id6: %u\n", id2Int[6]/8);
//            printf("id7: %u\n", id2Int[7]/8);
//            register __m256i xmm0 = _mm256_setzero_si256();
//            register __m256i xmm1 = _mm256_setzero_si256();
//            register __m256i xmm2 = _mm256_setzero_si256();
//            register __m256i xmm3 = _mm256_setzero_si256();
//            register __m256i xmm4 = _mm256_setzero_si256();
//            register __m256i xmm5 = _mm256_setzero_si256();
//            register __m256i xmm6 = _mm256_setzero_si256();
//            register __m256i xmm7 = _mm256_setzero_si256();
//            
//            subDiv1 = T;
//            subDiv2 = T + TSize;
//            subDiv3 = T + 2*TSize;
//            subDiv4 = T + 3*TSize;
//            
//            //            id = _mm256_slli_epi32 (id, int imm8);
//            
//            //prefetch C
//            //            xmm0 = _mm256_setzero_si256();
//            //            xmm1 = _mm256_setzero_si256();
//            //            xmm2 = _mm256_setzero_si256();
//            //            xmm3 = _mm256_setzero_si256();
//            //            xmm4 = _mm256_setzero_si256();
//            //            xmm5 = _mm256_setzero_si256();
//            //            xmm6 = _mm256_setzero_si256();
//            //            xmm7 = _mm256_setzero_si256();
//            
//            //            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv1 + (*(id1Int + 0))*n));
//            //            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv1 + (*(id1Int + 1))*n));
//            //            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv1 + (*(id1Int + 2))*n));
//            //            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv1 + (*(id1Int + 3))*n));
//            //            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv1 + (*(id1Int + 4))*n));
//            //            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv1 + (*(id1Int + 5))*n));
//            //            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv1 + (*(id1Int + 6))*n));
//            //            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv1 + (*(id1Int + 7))*n));
//            //            xmm0 = _mm256_load_si256 (subDiv1 + (_mm256_extract_epi32 (id1, 0))*n);
//            xmm0 = _mm256_load_si256 (subDiv1 + (*(id1Int + 0)));
//            xmm1 = _mm256_load_si256 (subDiv1 + (*(id1Int + 1)));
//            xmm2 = _mm256_load_si256 (subDiv1 + (*(id1Int + 2)));
//            xmm3 = _mm256_load_si256 (subDiv1 + (*(id1Int + 3)));
//            xmm4 = _mm256_load_si256 (subDiv1 + (*(id1Int + 4)));
//            xmm5 = _mm256_load_si256 (subDiv1 + (*(id1Int + 5)));
//            xmm6 = _mm256_load_si256 (subDiv1 + (*(id1Int + 6)));
//            xmm7 = _mm256_load_si256 (subDiv1 + (*(id1Int + 7)));//id1Int++
//            //            id1Int -= 8;
//            //            printf("C\n");
//            //            printBinaryMatrix(C, m, n);
//            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv2 + (*(id2Int + 0))));
//            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv2 + (*(id2Int + 1))));
//            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv2 + (*(id2Int + 2))));
//            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv2 + (*(id2Int + 3))));
//            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv2 + (*(id2Int + 4))));
//            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv2 + (*(id2Int + 5))));
//            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv2 + (*(id2Int + 6))));
//            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv2 + (*(id2Int + 7))));
//            
//            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv3 + (*(id3Int + 0))));
//            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv3 + (*(id3Int + 1))));
//            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv3 + (*(id3Int + 2))));
//            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv3 + (*(id3Int + 3))));
//            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv3 + (*(id3Int + 4))));
//            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv3 + (*(id3Int + 5))));
//            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv3 + (*(id3Int + 6))));
//            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv3 + (*(id3Int + 7))));
//            
//            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv4 + (*(id4Int + 0))));
//            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv4 + (*(id4Int + 1))));
//            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv4 + (*(id4Int + 2))));
//            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv4 + (*(id4Int + 3))));
//            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv4 + (*(id4Int + 4))));
//            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv4 + (*(id4Int + 5))));
//            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv4 + (*(id4Int + 6))));
//            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv4 + (*(id4Int + 7))));
//            //            C[8*n*j + n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[1] ));//wrong idInt
//            //            C[8*n*j + 2*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[2] ));//wrong idInt
//            //            C[8*n*j + 3*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[3] ));//wrong idInt
//            //            C[8*n*j + 4*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[4] ));//wrong idInt
//            //            C[8*n*j + 5*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[5] ));//wrong idInt
//            //            C[8*n*j + 6*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[6] ));//wrong idInt
//            //            C[8*n*j + 7*n + k] ^=  (_mm256_load_si256(T + idInt[7] ));//wrong idInt
//            
//            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + 8*n*j + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (xmm1, _mm256_load_si256 (C + 8*n*j + n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (xmm2, _mm256_load_si256 (C + 8*n*j + 2*n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (xmm3, _mm256_load_si256 (C + 8*n*j + 3*n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (xmm4, _mm256_load_si256 (C + 8*n*j + 4*n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (xmm5, _mm256_load_si256 (C + 8*n*j + 5*n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (xmm6, _mm256_load_si256 (C + 8*n*j + 6*n + 8*k)));
//            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (xmm7, _mm256_load_si256 (C + 8*n*j + 7*n + 8*k)));
//            
//            //            _mm256_storeu_si256(C + 8*n*j + 8*k, xmm0);
//            //            _mm256_storeu_si256(C + 8*n*j + n + 8*k, xmm1);
//            //            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, xmm2);
//            //            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, xmm3);
//            //            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, xmm4);
//            //            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, xmm5);
//            //            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, xmm6);
//            //            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, xmm7);
//            
//            //            subDiv1 += 8*TRows;
//            //            subDiv2 += 8*TRows;
//            //            subDiv3 += 8*TRows;
//            //            subDiv4 += 8*TRows;
//         }
////         for (int j = 0; j < rowsLeftover; j++){
////            unsigned int id1 = readBits(A, j*32 + i*32*m + 256*readsPerColumn, t); //  + i*32*m//
////            unsigned int id2 = readBits(A, t + j*32 + i*32*m + 256*readsPerColumn, t);
////            unsigned int id3 = readBits(A, 2*t + j*32 + i*32*m + 256*readsPerColumn, t);
////            unsigned int id4 = readBits(A, 3*t + j*32 + i*32*m + 256*readsPerColumn, t);
////            printf("id1: %d\n", id1);
////            printf("id2: %d\n", id2);
////            printf("id3: %d\n", id3);
////            printf("id4: %d\n", id4);
////            //                        printf("APostition: %d\n", k*t + j*32 + i*32*m);
////            //                        printf("id: %d, Table#: %d\n", id, k);
////            //                        printf("CPosition: %d\n", j*n + i);
////            //                        printf("TValue: ");
////            //                        printBinary(T[id + TSize*k]);
////            //                        printf("\n");
////            register __m256i xmm0 = _mm256_setzero_si256();
////            
////            xmm0 = _mm256_load_si256 (subDiv1 + id1*8);
////            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv2 + id2*8));
////            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv3 + id3*8));
////            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv4 + id4*8));
////            _mm256_storeu_si256(C + (j+8*readsPerColumn)*n + k*8, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + (j+8*readsPerColumn)*n + k*8)));
////         }
//      }
//   }
//   //   for (int k = 0; k < nr256; k += 8){//just go one for all of nr 256
//   if (nr256){
//      for (int i = 0; i < la; i++){
////         createTable(T, B, TRows, n, t, numberOfTables, i, nq256, nr256, nq256);
//         //       printBinaryMatrix(T, 8, 1);
//         //       printBinaryMatrix(T, TRows, n);
//         print4DBinaryMatrix(T, TRows, 0, numberOfTables, 0, nr256);
//         //         print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, innerColumns, int outerColumns);
//         
//         for (int j = 0; j < readsPerColumn; j++){
//            id4 = _mm256_load_si256(A + j*8 + i*m);
//            id1 = _mm256_srli_epi32 (id4, 24);
//            //            id1 = _mm256_slli_epi32(id1, 3);
//            register unsigned int* id1Int = (unsigned int*) &(id1);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id2 = _mm256_srli_epi32 (id4, 24);
//            //            id2 = _mm256_slli_epi32(id2, 3);
//            register unsigned int* id2Int = (unsigned int*) &(id2);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id3 = _mm256_srli_epi32 (id4, 24);
//            //            id3 = _mm256_slli_epi32(id3, 3);
//            register unsigned int* id3Int = (unsigned int*) &(id3);
//            id4 = _mm256_slli_epi32(id4, 8);
//            id4 = _mm256_srli_epi32 (id4, 24);
//            //            id4 = _mm256_slli_epi32(id4, 3);
//            register unsigned int* id4Int = (unsigned int*) &(id4);
//            printf("id0: %u\n", id1Int[0]);
//            printf("id1: %u\n", id1Int[1]);
//            printf("id2: %u\n", id1Int[2]);
//            printf("id3: %u\n", id1Int[3]);
//            printf("id4: %u\n", id1Int[4]);
//            printf("id5: %u\n", id1Int[5]);
//            printf("id6: %u\n", id1Int[6]);
//            printf("id7: %u\n", id1Int[7]);
//            printf("id0: %u\n", id2Int[0]);
//            printf("id1: %u\n", id2Int[1]);
//            printf("id2: %u\n", id2Int[2]);
//            printf("id3: %u\n", id2Int[3]);
//            printf("id4: %u\n", id2Int[4]);
//            printf("id5: %u\n", id2Int[5]);
//            printf("id6: %u\n", id2Int[6]);
//            printf("id7: %u\n", id2Int[7]);
//            for (int l = 0; l < nr256; l++){
//               C[8*n*j + 8*nq256 + l] ^= *(subDiv1 + *(id1Int + 0)); //T[id1Int[0]*n + l + 8*nq256];
//               C[8*n*j + 8*nq256 + n + l] ^= *(subDiv1 + *(id1Int + 1));
//               C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv1 + *(id1Int + 2));
//               C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv1 + *(id1Int + 3));
//               C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv1 + *(id1Int + 4));
//               C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv1 + *(id1Int + 5));
//               C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv1 + *(id1Int + 6));
//               C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv1 + *(id1Int + 7));
//               C[8*n*j + 8*nq256 + l] ^= *(subDiv2 + *(id2Int + 0));
//               C[8*n*j + 8*nq256 + n + l] ^= *(subDiv2 + *(id2Int + 1));
//               C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv2 + *(id2Int + 2));
//               C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv2 + *(id2Int + 3));
//               C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv2 + *(id2Int + 4));
//               C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv2 + *(id2Int + 5));
//               C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv2 + *(id2Int + 6));
//               C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv2 + *(id2Int + 7));
//               C[8*n*j + 8*nq256 + l] ^= *(subDiv3 + *(id3Int + 0));
//               C[8*n*j + 8*nq256 + n + l] ^= *(subDiv3 + *(id3Int + 1));
//               C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv3 + *(id3Int + 2));
//               C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv3 + *(id3Int + 3));
//               C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv3 + *(id3Int + 4));
//               C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv3 + *(id3Int + 5));
//               C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv3 + *(id3Int + 6));
//               C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv3 + *(id3Int + 7));
//               C[8*n*j + 8*nq256 + l] ^= *(subDiv4 + *(id4Int + 0));
//               C[8*n*j + 8*nq256 + n + l] ^= *(subDiv4 + *(id4Int + 1));
//               C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv4 + *(id4Int + 2));
//               C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv4 + *(id4Int + 3));
//               C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv4 + *(id4Int + 4));
//               C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv4 + *(id4Int + 5));
//               C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv4 + *(id4Int + 6));
//               C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv4 + *(id4Int + 7)); //T[id4Int[7]*n + l + 8*nq256 + 3*TSize]
//               
//               //            T[id2*n + TSize + l];
//               //            C[(j+readsPerColumn)*8*n + l] ^= T[id3*n + TSize*2 + l];
//               //            C[(j+readsPerColumn)*8*n + l] ^= T[id4*n + TSize*3 + l];
//               
//               subDiv1 += TRows;
//               subDiv2 += TRows;
//               subDiv3 += TRows;
//               subDiv4 += TRows;
//            }
//            subDiv1 = T;
//            subDiv2 = T + TSize;
//            subDiv3 = T + 2*TSize;
//            subDiv4 = T + 3*TSize;
//         }
//         for (int j = 0; j < rowsLeftover; j++){
//            unsigned int id1 = readBits(A, j*32 + i*32*m + 256*readsPerColumn, t); //  + i*32*m//
//            unsigned int id2 = readBits(A, t + j*32 + i*32*m + 256*readsPerColumn, t);
//            unsigned int id3 = readBits(A, 2*t + j*32 + i*32*m + 256*readsPerColumn, t);
//            unsigned int id4 = readBits(A, 3*t + j*32 + i*32*m + 256*readsPerColumn, t);
//            printf("id1: %d\n", id1);
//            printf("id2: %d\n", id2);
//            printf("id3: %d\n", id3);
//            printf("id4: %d\n", id4);
//            //                        printf("APostition: %d\n", k*t + j*32 + i*32*m);
//            //                        printf("id: %d, Table#: %d\n", id, k);
//            //                        printf("CPosition: %d\n", j*n + i);
//            //                        printf("TValue: ");
//            //                        printBinary(T[id + TSize*k]);
//            //                        printf("\n");
//            
//            subDiv1 = T;
//            subDiv2 = T + TSize;
//            subDiv3 = T + 2*TSize;
//            subDiv4 = T + 3*TSize;
//            for (int k = 0; k < nr256; k++){
//               C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv1 + k*TRows + id1);//T[id1 + l];
//               C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv2 + k*TRows + id2);//T[id2 + TSize + l];
//               C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv3 + k*TRows + id3);//T[id3 + TSize*2 + l];
//               C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv4 + k*TRows + id4);//T[id4 + TSize*3 + l];
//            }
//            subDiv1 += TRows;
//            subDiv2 += TRows;
//            subDiv3 += TRows;
//            subDiv4 += TRows;
//            
//         }
//         subDiv1 = T;
//         subDiv2 = T + TSize;
//         subDiv3 = T + 2*TSize;
//         subDiv4 = T + 3*TSize;
//         //         printf("C\n");
//         //         printBinaryMatrix(C, m, n);
//         
//      }
//   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", (endTime - startTime));
   
//   printf("C\n");
//   printBinaryMatrix(C, m, n);
   
   
   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   ///*
   //   startTime = rdtsc();
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
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