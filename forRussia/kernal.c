/* kernal.c resembles my final immplementation, but
 * runs with a bus error, so its probably not useful
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>
//#include <malloc.h>

#define VXOR(srca,srcb,dest)             \
__asm__ volatile(                               \
"vpxor %[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
: [vsrca] "x"(srca),                    \
[vsrcb] "x"(srcb)                     \
);

#define VLOAD(srcaddr,dest)           \
__asm__ volatile(                            \
"vmovdqu %[addr], %[rvdest]\n"       \
: [rvdest] "=x"(dest)                \
: [addr] "m"(*srcaddr)               \
);

#define and(srca, dest)        \
__asm__ volatile(                 \
"and  %[vsrca], %[rvdest]\n"       \
: [rvdest] "=x"(dest)                \
:[vsrca] "m"(srca)               \
);


static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

void print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns);
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
     _mm256_storeu_si256((T + i*TRows*8), _mm256_setzero_si256());
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
		 _mm256_storeu_si256((T + l*TRows*8 + i*8 + e*numOfTables*TRows*8), _mm256_xor_si256 (_mm256_loadu_si256((T + l*TRows*8 + (i - 1)*8 + e*numOfTables*TRows*8)), _mm256_loadu_si256((Bt + (l*t + k)*n + tableNumber*64*n + 8*e))));
               } else{
                  _mm256_storeu_si256(T + l*TRows*8 + i*8 + e*numOfTables*TRows*8, _mm256_xor_si256 (_mm256_loadu_si256 (T + l*TRows*8 + i*8 + e*numOfTables*TRows*8), _mm256_loadu_si256 (Bt + (l*t + k)*n + tableNumber*64*n + e*8)));
               }
            }
         }
      }
   }
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
void convertBinaryMatrixToFullSizeColumnMajor(unsigned int* fullSizeMatrix, long unsigned int* binaryMatrix, int rows, int columns, int binaryColumns){
   long unsigned int tmp = 0;
   int rowCounter = 0;
   for (int j = 0; j < binaryColumns; j++){
      for (int i = 0; i < rows; i++){
         tmp = binaryMatrix[i + j*rows];
         for (int k = 0; k < 64; k++){
            if (j*64 + k >= columns){
               break;
            }
            if (tmp & 0x8000000000000000){
               fullSizeMatrix[i*columns + j * 64 + k] = 1;
            } else {
               fullSizeMatrix[i*columns + j * 64 + k] = 0;
            }
            tmp <<= 1;
         }
      }
      //rowCounter++;
   }
//      printf("size: %d\n", rows*columns);
//      for (int i = 0; i < rows; i++){
//         for (int j = 0; j < columns; j++){
//            printf("%d", fullSizeMatrix[i*columns + j]);
//         }
//         printf("\t");
//         printf(" %d", rowCounter);
//         printf("\n");
//         rowCounter++;
//      }
//      printf("\sn");
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



void printBinaryColumnMatrix(long unsigned int* matrix, int rows, int columns){
   long unsigned tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         tmp = matrix[i + j*rows];
         for (int k = 0; k < 64; k++){
            if (tmp & 0x8000000000000000){
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
   int la = (lb - 1) / 64 + 1;  // columns in A
   int nActual = 512;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = (64 - (lb % 64)) % 64; // zeros to pad on A
   int rb = (64 - (nActual % 64)) % 64; // zeros to pad on B
      printf("m: %d\n", m);
      printf("lb: %d\n", lb);
      printf("la: %d\n", la);
      printf("n: %d\n", n);
      printf("ra: %d\n", ra);
      printf("rb: %d\n", rb);
   
   //   printf("t: %d\n", t);
      //   long unsigned int* A = malloc(m*la*(sizeof(long unsigned int)));
      long unsigned int* A = malloc(m*la*(sizeof(long unsigned int)));
   
   int lilExtra = (64 - (lb % 64)) % 64;
   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(long unsigned int))); // pad rows for createTable
   //   printf("SizeB: %d\n", (lb + lilExtra)*n);
   unsigned int* C = calloc(m*n,(sizeof(long unsigned int)));
   srand(time(NULL));
   
   for (int j = 0; j < la; j++){// init A column major order
      for (int i = 0; i < m; i++){
//         if (j < la - 1){
            A[i + j*m] = ((long unsigned) rand()) << 32;
            A[i + j*m] += rand();
//         } else {
//            if (ra > 32){
//               A[i + j*m] = ((long unsigned) rand()) << 32;
//               A[i + j*m] += rand() >> (ra - 32);
//               A[i + j*m] <<= (ra - 32);
//            } else {
//               A[i + j*m] = rand() >> ra;
//               A[i + j*m] <<= ra;
//            }
//         }
      }
   }
//      printf("A\n");
//      printBinaryColumnMatrix(A, m, la);
   
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
//      printf("B\n");
//      printBinaryMatrix(B, lb + lilExtra, n);
   
   int t = 8;  // four Russian divisor
   unsigned int mask = 0x7f8; //0x7f8
   unsigned int startingMask = 0xff;
   int TRows = 1 << t;
   int TSize = TRows * 8;
   int numberOfTables = 64 / t;
   
   int k = 256;
   int readsPerColumn = m*32 / k;
   int rowsLeftover = (m*32 % k) / 32;
   int nq512 = nActual / 512;
   int nq256 = nActual / 256;
   int nr256 = ((n - 1) % 8) + 1;
   if (!(nActual % 256)){
      nr256 = 0;
   }
   unsigned int* T = malloc(nq256*TSize*(numberOfTables)*sizeof(int));
   printf("number of tables: %d\n", numberOfTables);
   printf("TSize: %d\n", TSize);
   printf("FullTablesSize: %d\n", nq256*TSize*(numberOfTables)*sizeof(int));
   printf("B is %d * 256bits + %d ints\n", nq256, nr256);
   
   createTable(T, B, TRows, n, t, numberOfTables, 0, nq256, nr256);
   //       printBinaryMatrix(T, 8, 1);
   //       printBinaryMatrix(T, TRows, n);
//         print4DBinaryMatrix(T, TRows, 1, numberOfTables*nq256, 8, 0);
//         print4DBinaryMatrix(T+numberOfTables*TSize, TRows, 1, numberOfTables*nq256, 8, 0);
   //print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, innerColumns, int outerColumns
   
   __m256i xmm0, xmm1, xmm2, xmm3, ymm0, ymm1, ymm2, ymm3;
   __m256i xmm4, xmm5, xmm6, xmm7, ymm4, ymm5, ymm6, ymm7;
   register volatile  unsigned int* subDiv = T + (numberOfTables - 1)*TSize;
   register volatile unsigned int* subDiv2 = T + (numberOfTables - 1)*TSize + numberOfTables*TSize;
   //unsigned int* cPosition = C;
   
//         unsigned int* cIntermediateTest = calloc(m*n,(sizeof(int)));
   
   unsigned long long startTime = rdtsc();
   for (int k = 0; k < 1; k++){ //nActual / 256
     //      createTable(T, B, TRows, n, t, numberOfTables, 0, nq256, nr256);
     for (int j = 0; j < 1; j++){
        
       for (int i = 0; i < 1000; i++){
          register long unsigned int a0 = *(A + 0);
          register long unsigned int a1 = *(A + 1);
          register long unsigned int a2 = *(A + 2);
          register long unsigned int a3 = *(A + 3);
            
          register unsigned int index0 = (a0 & startingMask) << 3; //variable t (32 - ((ii+1) * t)))
          VLOAD((subDiv + index0), xmm0);
	  //          VLOAD((subDiv2 + index0), xmm4);
          register unsigned int index1 = (a1 & startingMask) << 3;// consider hardcoding more of these values
          VLOAD((subDiv + index1), xmm1);
	  //          VLOAD((subDiv2 + index1), xmm5);
          register unsigned int index2 = (a2 & startingMask) << 3;
          VLOAD((subDiv + index2), xmm2);
	  //          VLOAD((subDiv2 + index2), xmm6);
          register unsigned int index3 = (a3 & startingMask) << 3;
          VLOAD((subDiv + index3), xmm3);
	  //          VLOAD((subDiv2 + index3), xmm7);          

          /*
          xmm0 = _mm256_setzero_si256();
          xmm1 = _mm256_setzero_si256();
          xmm2 = _mm256_setzero_si256();
          xmm3 = _mm256_setzero_si256();
          xmm4 = _mm256_setzero_si256();
          xmm5 = _mm256_setzero_si256();
          xmm6 = _mm256_setzero_si256();
          xmm7 = _mm256_setzero_si256();
          */
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
//                      printf("id0: %u, full %u\n", index0/8, index0);
//                      printf("id1: %u, full %u\n", index1/8, index1);
//                      printf("id2: %u, full %u\n", index2/8, index2);
//                      printf("id3: %u, full %u\n", index3/8, index3);
//          
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, xmm0); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k + 8, xmm4);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, xmm1); _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k + 8, xmm5);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, xmm2); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k + 8, xmm6);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, xmm3); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k + 8, xmm7);
//                      printf("CIntermediate\n");
//                      printBinaryMatrix(cIntermediateTest, m, n);
          
          a0 = a0 >> 5; a1 = a1 >> 5; a2 = a2 >> 5; a3 = a3 >> 5;
          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);

          index1 = (a1 & mask);   a1 = a1 >> 8; 
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);

          index2 = (a2 & mask);   a2 = a2 >> 8; 
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);

          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);

          subDiv -= TSize;
          subDiv2 -= TSize;
                  
//                      printf("id0: %u, full %u\n", index0/8, index0);
//                      printf("id1: %u, full %u\n", index1/8, index1);
//                      printf("id2: %u, full %u\n", index2/8, index2);
//                      printf("id3: %u, full %u\n", index3/8, index3);
//                      
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k + 8, ymm4);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1); _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k + 8, ymm5);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k + 8, ymm6);
//                      _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3); _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k + 8, ymm7);
//                      printf("CIntermediate\n");
//                      printBinaryMatrix(cIntermediateTest, m, n);
          

          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);   a1 = a1 >> 8;
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);   a2 = a2 >> 8;
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
          //            printf("id0: %u, full %u\n", index0/8, index0);
          //            printf("id1: %u, full %u\n", index1/8, index1);
          //            printf("id2: %u, full %u\n", index2/8, index2);
          //            printf("id3: %u, full %u\n", index3/8, index3);
          //
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3);
          //            printf("CIntermediate\n");
          //            printBinaryMatrix(cIntermediateTest, m, n);
                  
          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);   a1 = a1 >> 8;
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);   a2 = a2 >> 8;
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
          //            printf("id0: %u, full %u\n", index0/8, index0);
          //            printf("id1: %u, full %u\n", index1/8, index1);
          //            printf("id2: %u, full %u\n", index2/8, index2);
          //            printf("id3: %u, full %u\n", index3/8, index3);
          //
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3);
          //            printf("CIntermediate\n");
          //            printBinaryMatrix(cIntermediateTest, m, n);
                  
          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);   a1 = a1 >> 8;
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);   a2 = a2 >> 8;
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
          //            printf("id0: %u, full %u\n", index0/8, index0);
          //            printf("id1: %u, full %u\n", index1/8, index1);
          //            printf("id2: %u, full %u\n", index2/8, index2);
          //            printf("id3: %u, full %u\n", index3/8, index3);
          //
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3);
          //            printf("CIntermediate\n");
          //            printBinaryMatrix(cIntermediateTest, m, n);
                  
          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);   a1 = a1 >> 8;
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);   a2 = a2 >> 8;
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
          //            printf("id0: %u, full %u\n", index0/8, index0);
          //            printf("id1: %u, full %u\n", index1/8, index1);
          //            printf("id2: %u, full %u\n", index2/8, index2);
          //            printf("id3: %u, full %u\n", index3/8, index3);
          //
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3);
          //            printf("CIntermediate\n");
          //            printBinaryMatrix(cIntermediateTest, m, n);
                  
          index0 = (a0 & mask); 	 a0 = a0 >> 8;
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);   a1 = a1 >> 8;
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);   a2 = a2 >> 8;
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);   a3 = a3 >> 8;
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          subDiv -= TSize;
          subDiv2 -= TSize;
                  
          //            printf("id0: %u, full %u\n", index0/8, index0);
          //            printf("id1: %u, full %u\n", index1/8, index1);
          //            printf("id2: %u, full %u\n", index2/8, index2);
          //            printf("id3: %u, full %u\n", index3/8, index3);
          //
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 8*k, ymm0);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + n + 8*k, ymm1);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 2*n + 8*k, ymm2);
          //            _mm256_storeu_si256(cIntermediateTest + 4*n*j + 3*n + 8*k, ymm3);
          //            printf("CIntermediate\n");
          //            printBinaryMatrix(cIntermediateTest, m, n);
                  
          index0 = (a0 & mask);
          VLOAD((subDiv + index0), ymm0);
          VLOAD((subDiv2 + index0), ymm4);
          VXOR(xmm0, ymm0, xmm0);
          VXOR(xmm4, ymm4, xmm4);
          
          index1 = (a1 & mask);
          VLOAD((subDiv + index1), ymm1);
          VLOAD((subDiv2 + index1), ymm5);
          VXOR(xmm1, ymm1, xmm1);
          VXOR(xmm5, ymm5, xmm5);
          
          index2 = (a2 & mask);
          VLOAD((subDiv + index2), ymm2);
          VLOAD((subDiv2 + index2), ymm6);
          VXOR(xmm2, ymm2, xmm2);
          VXOR(xmm6, ymm6, xmm6);
          
          index3 = (a3 & mask);
          VLOAD((subDiv + index3), ymm3);
          VLOAD((subDiv2 + index3), ymm7);
          VXOR(xmm3, ymm3, xmm3);
          VXOR(xmm7, ymm7, xmm7);
          
          //unsigned long long endTime = rdtsc();
          //printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
          //printf("cycles: %llu\n", (endTime - startTime));             

	  VLOAD(C, ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C, xmm0);

	  VLOAD((C+n), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + n, xmm1);

	  VLOAD((C+2*n), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + 2*n, xmm2);

	  VLOAD((C+3*n), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + 3*n, xmm3);

	  VLOAD((C+8), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + 8, xmm4);

	  VLOAD((C+n+8), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + n + 8, xmm5);

	  VLOAD((C+2*n+8), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + 2*n + 8, xmm6);

	  VLOAD((C+3*n+8), ymm0);
	  VXOR(ymm0, xmm0, xmm0);
          _mm256_storeu_si256(C + 3*n + 8, xmm7);

         } // i
         //cPosition += 8*n;
      } // j
      //cPostition += 8;
   } // k
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", (endTime - startTime));
   
//      printf("C\n");
//      printBinaryMatrix(C, m, n);
   
   
   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   ///*
   //   startTime = rdtsc();
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
//      printf("Af\n");
   convertBinaryMatrixToFullSizeColumnMajor(Af, A, m, lb, la);
//   for (int i = 0; i < m*lb; i++){
//      printf("%d", Af[i]);
//   } printf("\n");
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
