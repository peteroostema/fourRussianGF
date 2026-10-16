/* scriptTest , 2, 3.c are versions that compute after taking in
 * each dimension as a command line argument.
 * They work with scriptTest.sh to test what sizes they correclty
 * compute
 *
 * run by: ./a.out 512 512 512
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>

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
   int a = 7683478;
   //   a >>= 10;
   //   a <<= 10;
   //   printf("%d", a);
   //   printf("\n");
   //   printf("%d", ~a);
   //   printf("\n");
   
   //   a = 0x80000000;
   //   printBinary(a);
   
   //   a = rand(); printf("rand: %d\n", a); printBinary(a);
   //   a = rand(); printf("rand: %d\n", a); printBinary(a);
   //   a = rand(); printf("rand: %d\n", a); printBinary(a);
   //   a = rand(); printf("rand: %d\n", a); printBinary(a);
   //   a = rand(); printf("rand: %d\n", a); printBinary(a);
   
   
   //   for (int i = 1; i < 32; i++){
   //      printf("%d\n", findNextGreyRow(i, 5));
   //   }
   
   //   int t = 5;
   //   int size = 32;
   //   for (int i = 1; i < size; i++){
   //      unsigned int rowToXOR = 1;
   //      unsigned int position = i;
   //      while (rowToXOR <= t) {
   //         printf("position: %d\t", position);
   //         if (position % 2){
   //            break;
   //         }
   //         position /= 2;
   //         rowToXOR++;
   //      }
   //      printf("rowtoXor: %d\n", rowToXOR);
   //   }
   
   //   unsigned int b = 0;
   //   printf("%d\n", b);
   //   b = 2000000000;
   //   printf("%d\n", b);
   //   b = -2;
   //   printf("%d\n", b);
   
   //   unsigned int* A = malloc(8*sizeof(int));
   //   for (int i = 0; i < 8; i++) {
   //      A[i] = i;
   //   }
   //   printBinaryMatrix(A, 8, 1);
   //   printf("%d\n", readBits(A, 126, 4));
   //   printf("%d\n", readBits(A, 126, 2));
   //   free(A);
   
   //   int tmp =  rand();
   //   printf("tmp: %d\n", tmp);
   //   tmp >>= ra;
   //   printf("tmp: %d\n", tmp);
   //   tmp <<= ra;
   //   printf("tmp: %d\n", tmp);
   
   int m = atoi(argv[1]);  // rows in A
   int lb = atoi(argv[2]); // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = atoi(argv[3]);
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
//      printf("m: %d\n", m);
//      printf("lb: %d\n", lb);
//      printf("la: %d\n", la);
//      printf("n: %d\n", n);
//      printf("ra: %d\n", ra);
//      printf("rb: %d\n", rb);
   int t = 8;  // four Russian divisor
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
   //
   //   int TRows = power(2, t);
   //   int TSize = TRows * n;
   //   int numberOfTables = 32 / t;
   //   unsigned int* T = malloc(TSize*(numberOfTables)*sizeof(int));
   //   printf("subTableSize: %d\n", TSize);
   //   printf("numberOfTables: %d\n", numberOfTables);
   //   printf("Size of Table: %lu\n", TSize*(numberOfTables)*sizeof(int));
   //   createTable(T, B, TRows, TSize, n, t, numberOfTables);
   //   print3DBinaryMatrix(T, TRows, n, numberOfTables);
   //   //   sleep(10);
   //   free(T);
   
   //   for (int i = 0; i < m; i++){// init C // test
   //      for (int j = 0; j < n; j++){
   //         if (j < n - 1){
   //            C[i*n + j] = rand();
   //         } else {
   //            C[i*n + j] = rand();
   //         }
   //      }
   //   }
   //   printf("C\n");
   //   printBinaryMatrix(C, m, n);
   
   //   unsigned int* Af = malloc(m*lb*sizeof(int));
   //   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   //   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   //   printf("Af\n");
   //   convertBinaryMatrixToFullSize(Af, A, m, lb, la);
   //   printf("Bf\n");
   //   convertBinaryMatrixToFullSize(Bf, B, lb, nActual, n);
   //
   //   for (int i = 0; i < m; i++){ // row
   //      for (int j = 0; j < nActual; j++){ // column
   //         for (int k = 0; k < lb; k++){ // interate
   //            //C[i * n + j] ^= (readBit(A, i * la * 32 + k) & readBit(B, j + k * n));
   //            Cf[i * nActual + j] ^= (Af[i * la + k] & Bf[j + k * nActual]);
   //         }
   //      }
   //   }
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
   
   //   printf("%d", 0*4 % 32);
   
   //   unsigned long long time = rdtsc();
   //   for (int i = 0; i < 200000000; i++){
   //      for (int j = 0; j < 20; j++){
   //               a++;
   //      }
   //   }
   
   //   unsigned long long laterTime = rdtsc();
   //   printf("time: %llu\n", time);
   //   printf("time: %llu\n", laterTime);
   //   printf("time: %f\n", (float)(laterTime - time) / 2800000000.0);
   //   printf("time: %llu\n", (laterTime - time) % 2800000000);
   //   printf("CLOCKS_PER_SEC: %d\n", CLOCKS_PER_SEC);
   
   
   
   //   /* Initialize the two argument vectors */
   //   __m256 evens = _mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
   //   __m256 odds = _mm256_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0);
   //
   //   /* Compute the difference between the two vectors */
   //   __m256 result = _mm256_sub_ps(evens, odds);
   //
   //   /* Display the elements of the result vector */
   //   float* f = (float*)&result;
   //   printf("%f %f %f %f %f %f %f %f\n",
   //            f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
   //
   //   return 0;
   
   int TRows = 1 << t;
   int TSize = TRows * n;
   int numberOfTables = 32 / t;
   unsigned int* T = malloc(TSize*(numberOfTables)*sizeof(int));
   int k = 256;
   int readsPerColumn = m*32 / k;
   int rowsLeftover = (m*32 % k) / 32;
   int nq256 = nActual / 256;
   int nr256 = ((n - 1) % 8) + 1;
   if (!(nActual % 256)){
      nr256 = 0;
   }
   __m256i id1;
   __m256i id2;
   __m256i id3;
   __m256i id4;
   //   printf("256 bit reads: %d\n", readsPerColumn);
   //   printf("leftOver: %d\n", rowsLeftover);
   //   printf("B is %d * 256bits + %d ints\n", nq256, nr256);
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i, nq256, nr256);
      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
      for (int j = 0; j < readsPerColumn; j++){
         id4 = _mm256_load_si256(A + j*8 + i*m);
         //         unsigned int* idInt = (unsigned int*) &id4;
         //         printf("Integers: \n");
         //         printBinary(idInt[0]); printBinary(idInt[1]); printBinary(idInt[2]); printBinary(idInt[3]); printBinary(idInt[4]); printBinary(idInt[5]); printBinary(idInt[6]); printBinary(idInt[7]);
         //         printf("IntegersAgain killme: %u %u %u %u %u %u %u %u\n", idInt[0], idInt[1], idInt[2], idInt[3], idInt[4], idInt[5], idInt[6], idInt[7]);
         //__m256i id = _mm256_load_ps(A + j/8 + i*m);
         id1 = _mm256_srli_epi32 (id4, 24);
         unsigned int* id1Int = (unsigned int*) &(id1);
         id4 = _mm256_slli_epi32(id4, 8);
         id2 = _mm256_srli_epi32 (id4, 24);
         unsigned int* id2Int = (unsigned int*) &(id2);
         id4 = _mm256_slli_epi32(id4, 8);
         id3 = _mm256_srli_epi32 (id4, 24);
         unsigned int* id3Int = (unsigned int*) &(id3);
         id4 = _mm256_slli_epi32(id4, 8);
         id4 = _mm256_srli_epi32 (id4, 24);
         unsigned int* id4Int = (unsigned int*) &(id4);
         //         printf("id1: \n");
         //         printBinary(id1Int[0]); printBinary(id1Int[1]); printBinary(id1Int[2]); printBinary(id1Int[3]); printBinary(id1Int[4]); printBinary(id1Int[5]); printBinary(id1Int[6]); printBinary(id1Int[7]);
         //         printf("id2: \n");
         //         printBinary(id2Int[0]); printBinary(id2Int[1]); printBinary(id2Int[2]); printBinary(id2Int[3]); printBinary(id2Int[4]); printBinary(id2Int[5]); printBinary(id2Int[6]); printBinary(id2Int[7]);
         //         printf("id3: \n");
         //         printBinary(id3Int[0]); printBinary(id3Int[1]); printBinary(id3Int[2]); printBinary(id3Int[3]); printBinary(id3Int[4]); printBinary(id3Int[5]); printBinary(id3Int[6]); printBinary(id3Int[7]);
         //         printf("id4: \n");
         //         printBinary(id4Int[0]); printBinary(id4Int[1]); printBinary(id4Int[2]); printBinary(id4Int[3]); printBinary(id4Int[4]); printBinary(id4Int[5]); printBinary(id4Int[6]); printBinary(id4Int[7]);
         //                  printf("id0: %u\n", id1Int[0]);
         //                  printf("id1: %u\n", id1Int[1]);
         //                  printf("id2: %u\n", id1Int[2]);
         //                  printf("id3: %u\n", id1Int[3]);
         //                  printf("id4: %u\n", id1Int[4]);
         //                  printf("id5: %u\n", id1Int[5]);
         //                  printf("id6: %u\n", id1Int[6]);
         //                  printf("id7: %u\n", id1Int[7]);
         //                  printf("id0: %u\n", id2Int[0]);
         //                  printf("id1: %u\n", id2Int[1]);
         //                  printf("id2: %u\n", id2Int[2]);
         //                  printf("id3: %u\n", id2Int[3]);
         //                  printf("id4: %u\n", id2Int[4]);
         //                  printf("id5: %u\n", id2Int[5]);
         //                  printf("id6: %u\n", id2Int[6]);
         //                  printf("id7: %u\n", id2Int[7]);
         for (int k = 0; k < nq256; k++){ //nActual / 256
            //            id = _mm256_slli_epi32 (id, int imm8);
            
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 8*k), _mm256_load_si256 (T + id1Int[0]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + n + 8*k), _mm256_load_si256 (T + id1Int[1]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 2*n + 8*k), _mm256_load_si256 (T + id1Int[2]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 3*n + 8*k), _mm256_load_si256 (T + id1Int[3]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 4*n + 8*k), _mm256_load_si256 (T + id1Int[4]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 5*n + 8*k), _mm256_load_si256 (T + id1Int[5]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 6*n + 8*k), _mm256_load_si256 (T + id1Int[6]*n + k*8)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 7*n + 8*k), _mm256_load_si256 (T + id1Int[7]*n + k*8)));
            //            printf("C\n");
            //            printBinaryMatrix(C, m, n);
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 8*k), _mm256_load_si256 (T + id2Int[0]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + n + 8*k), _mm256_load_si256 (T + id2Int[1]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 2*n + 8*k), _mm256_load_si256 (T + id2Int[2]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 3*n + 8*k), _mm256_load_si256 (T + id2Int[3]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 4*n + 8*k), _mm256_load_si256 (T + id2Int[4]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 5*n + 8*k), _mm256_load_si256 (T + id2Int[5]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 6*n + 8*k), _mm256_load_si256 (T + id2Int[6]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 7*n + 8*k), _mm256_load_si256 (T + id2Int[7]*n + k*8 + TSize)));
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 8*k), _mm256_load_si256 (T + id3Int[0]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + n + 8*k), _mm256_load_si256 (T + id3Int[1]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 2*n + 8*k), _mm256_load_si256 (T + id3Int[2]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 3*n + 8*k), _mm256_load_si256 (T + id3Int[3]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 4*n + 8*k), _mm256_load_si256 (T + id3Int[4]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 5*n + 8*k), _mm256_load_si256 (T + id3Int[5]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 6*n + 8*k), _mm256_load_si256 (T + id3Int[6]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 7*n + 8*k), _mm256_load_si256 (T + id3Int[7]*n + k*8 + 2*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 8*k), _mm256_load_si256 (T + id4Int[0]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + n + 8*k), _mm256_load_si256 (T + id4Int[1]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 2*n + 8*k), _mm256_load_si256 (T + id4Int[2]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 3*n + 8*k), _mm256_load_si256 (T + id4Int[3]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 4*n + 8*k), _mm256_load_si256 (T + id4Int[4]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 5*n + 8*k), _mm256_load_si256 (T + id4Int[5]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 6*n + 8*k), _mm256_load_si256 (T + id4Int[6]*n + k*8 + 3*TSize)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (_mm256_load_si256 (C + 8*n*j + 7*n + 8*k), _mm256_load_si256 (T + id4Int[7]*n + k*8 + 3*TSize)));
            //            C[8*n*j + n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[1] ));//wrong idInt
            //            C[8*n*j + 2*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[2] ));//wrong idInt
            //            C[8*n*j + 3*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[3] ));//wrong idInt
            //            C[8*n*j + 4*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[4] ));//wrong idInt
            //            C[8*n*j + 5*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[5] ));//wrong idInt
            //            C[8*n*j + 6*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[6] ));//wrong idInt
            //            C[8*n*j + 7*n + k] ^=  (_mm256_load_si256(T + idInt[7] ));//wrong idInt
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
         //         printf("C\n");
         //         printBinaryMatrix(C, m, n);
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
//   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   
   //      printf("C\n");
   //      printBinaryMatrix(C, m, n);
   //
   //   __m256i idk = _mm256_load_si256 (B + (0*t + 0)*n + 0*8 + 0*32*n);
   //   unsigned int* idk2 = (unsigned int*) &(idk);
   //   printf("BfirstRow: \n");
   //   printBinary(idk2[0]); printBinary(idk2[1]); printBinary(idk2[2]); printBinary(idk2[3]); printBinary(idk2[4]); printBinary(idk2[5]); printBinary(idk2[6]); printBinary(idk2[7]);
   //   printf("B is %d * 256bits + %d long\n", nActual / 256, nActual % 256);
   //   __m256i tstuff = _mm256_load_si256 (T + (0 - 1)*0 + 0*8 + 0*TRows*n);
   //   unsigned int* tInt = (unsigned int*) &(tstuff);
   //   printf("TFirstRow: \n");
   //   printBinary(tInt[0]); printBinary(tInt[1]); printBinary(tInt[2]); printBinary(tInt[3]); printBinary(tInt[4]); printBinary(tInt[5]); printBinary(tInt[6]); printBinary(tInt[7]);
   
   //   unsigned int* D = malloc(8*sizeof(int));
   //   __m256i d = _mm256_load_si256(A);
   //   D = (unsigned int*)&d;
   //   printBinaryMatrix(D, 4, 1);
   
   //   err_t bli_init( void );
   //   err_t bli_finalize( void );
   
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
   
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            //C[i * n + j] ^= (readBit(A, i * la * 32 + k) & readBit(B, j + k * n));
            Cf[i * nActual + j] ^= (Af[i * lb + k] & Bf[j + k * nActual]);
         }
      }
   }
   
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