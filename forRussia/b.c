/* b.c was used for random testing to put later
 * into other files, each test typically stands
 * on its own and is commented out as a block.
 * Some of these tests include: updating the creatTable
 * function,
 * properly padding matrixes, A, B, C
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>

#define PI 3.1415
#define circleArea(r) (PI*r*r)

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

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

static inline void createTableOld(unsigned int* T, unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
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

static inline void createTable(unsigned int* T, unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){// get rid of n store column majorish
   int leftoverStart = nq256*8*TRows;
   for (int i = 0; i < numOfTables; i++){
      for (int j = 0; j < nq256; j++){// rewrite to init all 256 bits things to zero // check print 3d for 256x2^txnx32/t
         _mm256_storeu_si256(T + i*TRows*n + j*TRows*8, _mm256_setzero_si256());
      }
      for (int j = 0; j < nr256; j++){
         T[leftoverStart + i*TRows*n + j*TRows] = 0;
      }
   }
   printf("init zero\n");
   //   print4DBinaryMatrix(T, TRows, nq256, numOfTables, 8, nr256);
   //   printBinaryMatrix(T, 8, 1);
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
      //            printf("rowtoXor: %d\n", rowsToXOR);
      for (int k = t - 1; k >= rowsToXOR; k--){
         for (int j = 0; j < nq256; j++){
            for (int l = 0; l < numOfTables; l++){
               if (k == t - 1){
                  _mm256_storeu_si256(T + l*TRows*n + i*8 + j*TRows*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*n + (i - 1)*8 + j*TRows*8), _mm256_load_si256 (Bt + (l*t + k)*n + j*8 + tableNumber*32*n)));
                  //                  T[i*n + j + l*TRows*n] = T[(i - 1)*n + j + l*TRows*n] ^ Bt[(l*t + k)*n + j + tableNumber*32*n];
               } else{
                  //                  T[i*n + j + l*TRows*n] ^= Bt[(l*t + k)*n + j + tableNumber*32*n];
                  _mm256_storeu_si256(T + l*TRows*n + i*8 + j*TRows*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*n + i*8 + j*TRows*8), _mm256_load_si256 (Bt + (l*t + k)*n + j*8 + tableNumber*32*n)));
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
                  T[l*TRows*n + i + j*TRows + leftoverStart] = T[l*TRows*n + (i-1) + j*TRows + leftoverStart] ^ Bt[(l*t + k)*n + j + tableNumber*32*n + nq256*8];
               } else{
                  T[l*TRows*n + i + j*TRows + leftoverStart] ^= Bt[(l*t + k)*n + j + tableNumber*32*n + nq256*8];// consider j in B since it goes one by one row
               }
            }
         }
      }
   }
   //   printf("TFirstRow?:\n");
   //   printBinaryMatrix(T, TRows, n);
}

int main(){
//   for (int i = 1; i < 258; i++){
//      printf("i: %d  r: ", i);
//      printf("%d%d\n", (i + 31) / 32 % 9, (((i % 256 - 1) / 32 + 1 - 1) % 8) + 1);
//   }
   
//   srand(time(NULL));
//   unsigned int* A = malloc(8*sizeof(int));
//   unsigned int* B1 = malloc(8*sizeof(int));
//   unsigned int* B2 = malloc(8*sizeof(int));
//   unsigned int* B3 = malloc(8*sizeof(int));
//   unsigned int* B4 = malloc(8*sizeof(int));
//   for (int i = 0; i < 8; i++){
//      A[i] = rand();
//   }
//   for (int i = 0; i < 8; i++){
//      B1[i] = 0x80808000;
//      B1[i] +=  i*4;
//      B2[i] = 0x80808000;
//      B2[i] +=  i*4 + 1;
//      B3[i] = 0x80808000;
//      B3[i] +=  i*4 + 2;
//      B4[i] = 0x80808000;
//      B4[i] +=  i*4 + 3;
//   }
////   __m256i a = _mm256_setzero_si256();
////   __m256i b = _mm256_setzero_si256();
////   __m256i c = _mm256_setzero_si256();
//   __m256i a = _mm256_load_si256(A);
//   __m256i b1 = _mm256_load_si256(B1);
//   __m256i b2 = _mm256_load_si256(B2);
//   __m256i b3 = _mm256_load_si256(B3);
//   __m256i b4 = _mm256_load_si256(B4);
//   __m256i c1 = _mm256_shuffle_epi8(a, b1);
//   __m256i c2 = _mm256_shuffle_epi8(a, b2);
//   __m256i c3 = _mm256_shuffle_epi8(a, b3);
//   __m256i c4 = _mm256_shuffle_epi8(a, b4);
//   unsigned int* aa = (unsigned int*) &(a);
////   unsigned int* bb = (unsigned int*) &(b);
//   unsigned int* C1 = (unsigned int*) &(c1);
//   unsigned int* C2 = (unsigned int*) &(c2);
//   unsigned int* C3 = (unsigned int*) &(c3);
//   unsigned int* C4 = (unsigned int*) &(c4);
//   printf("A\n");
//   printBinaryMatrix(aa, 8, 1);
////   printf("B\n");
////   printBinaryMatrix(bb, 8, 1);
//   printf("C1\n");
//   printBinaryMatrix(C1, 8, 1);
//   printf("C2\n");
//   printBinaryMatrix(C2, 8, 1);
//   printf("C3\n");
//   printBinaryMatrix(C3, 8, 1);
//   printf("C4\n");
//   printBinaryMatrix(C4, 8, 1);
//   printf("extract: %d\n", _mm256_extract_epi32 (a, 3));
//   
//   _mm256_storeu_si256(aa, _mm256_setzero_si256());
//   printBinaryMatrix(aa, 8, 1);
   
//   int lb = 16384; // rows in B
//   int nActual = 16384;
//   int n = (nActual - 1) / 32 + 1;  // columns in B
//   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
//   int t = 8;  // four Russian divisor
//   
//   int lilExtra = (32 - (lb % 32)) % 32;
//   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
//   srand(time(NULL));
//   
//   for (int i = 0; i < lb; i++){// init B
//      for (int j = 0; j < n; j++){
//         if (j < n - 1){
//            B[i*n + j] = rand();
//         } else {
//            B[i*n + j] = rand() >> rb;
//            B[i*n + j] <<= rb;
//         }
//      }
//      for (int i = 0; i < lilExtra; i++){
//         for (int j = 0; j < n; j++){
//            B[(i+lb)*n + j] = 0;
//         }
//      }
//   }
//   int TRows = 1 << t;
//   int TSize = TRows * n;
//   int numberOfTables = 32 / t;
//   int nq256 = nActual / 256;
//   int nr256 = ((n - 1) % 8) + 1;
//   if (!(nActual % 256)){
//      nr256 = 0;
//   }
//   unsigned int* T = malloc(TSize*(numberOfTables)*sizeof(int));
//   
//   unsigned long long startTime = rdtsc();
//   createTable(T, B, TRows, n, t, numberOfTables, 1, nq256, nr256);
//   unsigned long long endTime = rdtsc();
//   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
//   
//   startTime = rdtsc();
//   createTable(T, B, TRows, n, t, numberOfTables, 1, nq256, nr256);
//   endTime = rdtsc();
//   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   
   printf("%f\n", circleArea(2));
}