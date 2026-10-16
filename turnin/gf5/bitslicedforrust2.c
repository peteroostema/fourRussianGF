/* bitslicedforrust2.c is a bitsliced Four Russian matrix multiply with t = 2
 * Since it isn't unrolled and doesn't use intrinsics I've included full 
 * documentation in this file
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printGF5(long unsigned int r0, long unsigned int r1, long unsigned int r2);

/* Make a table of 64 rows of B sectioned off by t
 * @param: T, lookup table
 * @param: Bt, B (2nd factor in matrix multiply)
 * @param: TRows, 5^t
 * @param: n, matrix dimension n
 * @param: t, four Russian divisor
 * @param: numOfTables, 64/t # of sections of every combination of t rows of B
 * @param: tableNumber, which section of 64 row of B are being looked at
 * @param: nq256, # of 256 element sections n has (used in intrinsic solutions)
 * @param: nr256, remainder elements after 256 sections (used in intrinsic solutions)
 * result: each subTable in T contains every possible GF5 addition of t rows of B,
 *         T contains all needed info from 64 rows of B
 */
static inline void createTable(long unsigned int* T, long unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){// init the first row of each subTable to 0
      for (int j = 0; j < n; j++){
         *(T + (i*TRows + j)*3 + 0) = 0;
         *(T + (i*TRows + j)*3 + 1) = 0;
         *(T + (i*TRows + j)*3 + 2) = 0;
      }
   }
   int index = 0;
   
   for (int i = 0; i < numOfTables; i++){ // work through each subTable primarily
      long unsigned int b00, b01, b02, b10, b11, b12, t0, t1, t2, s0, s1, s2, s3;
//      printf("subtable: %d,   index: %d\n", i, index);
      
      for (int o = 0; o < n; o++){
         index = o*3 + i*TRows*n*3; // index of T to start with
         // this implementation has t = 2 so read in two rows of B
         b00 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 0]; // b0 is the row after b1 because A is read backwards,
         b01 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 1]; // that way you don't rearange the two elements read
         b02 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 2]; // from A at once.
         b10 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 0];
         b11 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 1];
         b12 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 2];
         for (int l = 0; l < 5; l++){// add seconed row after every 5 of the first
            for (int m = 0; m < 5; m++){ // add first row of B 5 times
               t0 = T[(index + 0)];//read in a set of 64 elements
               t1 = T[(index + 1)];                                       // Note: It's possible to create this table with fewer
               t2 = T[(index + 2)];                                       // operations by doing one addition for the first row
//               printGF5(t0, t1, t2);                                    // a double for the second and negating the two created
               index += 3*n; // results goes one row down                 // for the remaining needed two, but this works only
                                                                          // when 0 of the second row has been added. So this makes
               s0 = t0 ^ b00;//bit sliced addition                        // the creation take just as long using this fact from
               s1 = t1 ^ b01;                                             // what I've seen
               s2 = t2 ^ b02;
               s2 ^= s1 & (t0 & b00);
               s1 ^= t0 & b00;
               s2 ^= t1 & b01;
               s3 = t2 & b02;
               T[(index + 2)] = (s2 & ~s1 & ~s0); //my meathod
               T[(index + 1)] = (s3) | (s1 & s0) | (~s2 & s1);
               T[(index + 0)] = (~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2);
            } // adding 5 times should return the row to 0
            //   so row 5 will be exactly the second row read from B
            t0 = T[(index + 0)];
            t1 = T[(index + 1)];
            t2 = T[(index + 2)];
            s0 = t0 ^ b10;// bit sliced addition
            s1 = t1 ^ b11;
            s2 = t2 ^ b12;
            s2 ^= s1 & (t0 & b10);
            s1 ^= t0 & b10;
            s2 ^= t1 & b11;
            s3 = t2 & b12;
            T[(index + 2)] = (s2 & ~s1 & ~s0); //my meathod
            T[(index + 1)] = (s3) | (s1 & s0) | (~s2 & s1);
            T[(index + 0)] = (~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2);
         }
      }
   }
//   printf("\nfinal index: %d\n", index);
}

/* print out a lookup table (most of the time only 3 dimensions of the possilbe 4 are used)
 * @param: matrix, lookup table
 * @param: rows, # of rows in each subTable
 * @param: majorColumn, not used in 3D (once used to specify how many 256 columns the matrix has)
 * @param: innerColumns, # of register length columns for a 3D matrix (# of integers in 256 bit block for 4D)
 * @param: outerColumns, not used in 3D (remainder integer columns for 4D)
 * result: table printed out one subTable after another
 */
void print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   int subtableSize = rows*(majorColumns*innerColumns + outerColumns);
   for (int l = 0; l < depth; l++){ // subTables are printed out one after another
      printf("subTable: %d\nsubtableSize: %d\n", l, subtableSize);
      
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < innerColumns; j++){ // print an entire row before another
            tmp0 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 0];
            tmp1 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 1];
            tmp2 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 2];
            for (int i = 63; i >= 0; i--){
               long unsigned print0 = ((tmp0 >> i) & 1);
               long unsigned print1 = ((tmp1 >> i) & 1);
               long unsigned print2 = ((tmp2 >> i) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               printf("%lu", printInt); //"   %lu"
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

/* print out a matrix of GF5 elements
 *  only different from printBinaryMatrix as it needs to read from 3
 *  integers at a time and add up the matching bits
 */
void printGF5Matrix(long unsigned int* matrix, int rows, int columns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         tmp0 = matrix[(i*columns + j)*3 + 0];
         tmp1 = matrix[(i*columns + j)*3 + 1];
         tmp2 = matrix[(i*columns + j)*3 + 2];
         for (int i = 63; i >= 0; i--){
            long unsigned print0 = ((tmp0 >> i) & 1);
            long unsigned print1 = ((tmp1 >> i) & 1);
            long unsigned print2 = ((tmp2 >> i) & 1);
            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
            printf("%lu", printInt);
         }
         printf("\t");
      }
      printf(" %d", rowCounter);
      rowCounter++;
      printf("\n");
   }
   printf("\n");
}

/* print out the binary representation of a 64 bit integer
 */
void printBinary(long unsigned int number){
   long unsigned int tmp = number;
   for (int k = 0; k < 64; k++){
      if (tmp & 0x8000000000000000){
         printf("1");
      } else {
         printf("0");
      }
      tmp <<= 1;
   }
   printf("\n");
}

/* print out 64 bitsliced elements in GF5
 * @param: r0, r1, r2, the first second and third bit of 64 elements
 * result: 64 GF5 elements printed out right next to each other
 */
void printGF5(long unsigned int r0, long unsigned int r1, long unsigned int r2){
   long unsigned tmp0 = r0;
   long unsigned tmp1 = r1;
   long unsigned tmp2 = r2;
   for (int i = 63; i >= 0; i--){
      long unsigned print0 = ((tmp0 >> i) & 1);
      long unsigned print1 = ((tmp1 >> i) & 1);
      long unsigned print2 = ((tmp2 >> i) & 1);
      long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("%lu", printInt); //"   %lu"
   }
   printf("\n");
}

/* converts a GF5 matrix stored in bitsliced row major order in 64 bit integers
 * @param: fullSizeMatrix, matrix where each elements stored in a seperate unsignedint
 * @param: binaryMatrix, matrix with sliced elements
 * @param: rows, # of rows in the matrixes
 * @param: columns, # of columns in the matrixes
 * @param: binaryColumns, # of 32 bit unsignedints to store a row of binaryMatrix
 * result: fullSizeMatrix contains every element from binaryMatrix stored in a seperate unsignedint
 */
void convertGF5MatrixToFullSize(unsigned int* fullSizeMatrix, long unsigned int* gf5Matrix, int rows, int columns, int integerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp0 = gf5Matrix[(i*integerColumns + j)*3 + 0];// read 64
         tmp1 = gf5Matrix[(i*integerColumns + j)*3 + 1];// elements from
         tmp2 = gf5Matrix[(i*integerColumns + j)*3 + 2];// the matrix
         for (int k = 63; k >= 0; k--){
            long unsigned print0 = ((tmp0 >> k) & 1);
            long unsigned print1 = ((tmp1 >> k) & 1);
            long unsigned print2 = ((tmp2 >> k) & 1);
            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0; // add them up
            fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
         }
      }
   }
}

/* check the results from a definition and a four Russian matrix multiplication
 * @param: rows, # of rows in the matrixes
 * @param: columns, # of columns in the matrixes
 * @param: integerColumns, # of columns divided up into 64 bit integers
 */
int checkCorrectResult(long unsigned int* C, unsigned int* Cf, int rows, int columns, int integerColumns){
   long unsigned int tmp = 0;
   long unsigned int tmp0, tmp1, tmp2;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp0 = C[(i*integerColumns + j)*3 + 0];
         tmp1 = C[(i*integerColumns + j)*3 + 1];
         tmp2 = C[(i*integerColumns + j)*3 + 2];
         for (int k = 63; k >= 0; k--){
            long unsigned print0 = ((tmp0 >> k) & 1);
            long unsigned print1 = ((tmp1 >> k) & 1);
            long unsigned print2 = ((tmp2 >> k) & 1);
            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
            
            if (j*64 + (64 - k) >= columns){
               break;
            }
            if (!(printInt == Cf[i*columns + j * 64 + (63 - k)])){
               printf("row: %d\ncolumn: %d\n", i, j);
               printf("C: %lu\n", printInt);
               return 0;
            }
         }
      }
   }
   return 1;
}

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

int main(){
   srand(time(NULL));
   int m = 32;  // rows in A, and C
   int lb = 32; // rows in B
   int la = (lb - 1) / 64 + 1;  // integer columns in A
   int nActual = 32; // columns in B, and C
   int n = (nActual - 1) / 64 + 1;  // integer columns in B
   int t = 2;  // four Russian divisor
   
   printf("m: %d\n", m);
   printf("lb: %d\n", lb);
   printf("la: %d\n", la);
   printf("n: %d\n", n);
   printf("t: %d\n", t);
   int lPadding = (64 - (lb % 64)) % 64; // pad to the end of a 64 bit integer
   int nPadding = (64 - (nActual % 64)) % 64; // ''
   printf("Lpadding: %d\n", lPadding);
   printf("Npadding: %d\n", nPadding);
   
   long unsigned int* A = malloc(m*la*3*(sizeof(long unsigned int)));
   long unsigned int* B = malloc((lb + lPadding)*n*3*(sizeof(long unsigned int)));
   long unsigned int* C = calloc(m*n*3, (sizeof(long unsigned int)));
   
   for (int i = 0; i < m; i++){// init A
      for (int j = 0; j < la; j++){
         long unsigned int a2 = 0, a1 = 0, a0 = 0;
         if (j < la - 1){
            for (int i = 0; i < 64; i++){
               long unsigned int tmp = rand() % 5;
//               printf("%lu", tmp);
               a2 <<= 1;
               a1 <<= 1;
               a0 <<= 1;
               if (tmp & 0b100){
                  a2 |= 1;
               }
               if (tmp & 0b010){
                  a1 |= 1;
               }
               if (tmp & 0b001){
                  a0 |= 1;
               }
            }
         } else{
            for (int i = 0; i < 64 - lPadding; i++){
               long unsigned int tmp = rand() % 5;
//               printf("%lu", tmp);
               a2 <<= 1;
               a1 <<= 1;
               a0 <<= 1;
               if (tmp & 0b100){
                  a2 |= 1;
               }
               if (tmp & 0b010){
                  a1 |= 1;
               }
               if (tmp & 0b001){
                  a0 |= 1;
               }
            }
            a2 <<= lPadding;
            a1 <<= lPadding;
            a0 <<= lPadding;
//            printf("\n");
//            printBinary(a2);
//            printBinary(a1);
//            printBinary(a0);
//            printf("\n");
         }
         A[(i*la + j)*3 + 2] = a2;
         A[(i*la + j)*3 + 1] = a1;
         A[(i*la + j)*3 + 0] = a0;
      }
   }
//   printf("A\n");
//   printGF5Matrix(A, m, la);
   
   for (int i = 0; i < lb; i++){// init B
      for (int j = 0; j < n; j++){
         long unsigned int b2 = 0, b1 = 0, b0 = 0;
         if (j < n - 1){
            for (int i = 0; i < 64; i++){
               long unsigned int tmp = rand() % 5;
//               printf("%lu", tmp);
               b2 <<= 1;
               b1 <<= 1;
               b0 <<= 1;
               if (tmp & 0b100){
                  b2 |= 1;
               }
               if (tmp & 0b010){
                  b1 |= 1;
               }
               if (tmp & 0b001){
                  b0 |= 1;
               }
            }
         } else{
            for (int i = 0; i < 64 - nPadding; i++){
               long unsigned int tmp = rand() % 5;
//               printf("%lu", tmp);
               b2 <<= 1;
               b1 <<= 1;
               b0 <<= 1;
               if (tmp & 0b100){
                  b2 |= 1;
               }
               if (tmp & 0b010){
                  b1 |= 1;
               }
               if (tmp & 0b001){
                  b0 |= 1;
               }
            }
            b2 <<= nPadding;
            b1 <<= nPadding;
            b0 <<= nPadding;
//                        printf("\n");
//                        printBinary(b2);
//                        printBinary(b1);
//                        printBinary(b0);
//                        printf("\n");
         }
         B[(i*n + j)*3 + 2] = b2;
         B[(i*n + j)*3 + 1] = b1;
         B[(i*n + j)*3 + 0] = b0;
      }
   }
//   printf("B\n");
//   printGF5Matrix(B, lb, n);
   
   int TRows = 1; // 5^t rows to store every combination of t rows of B
   for (int i = 0; i < t; i++){
      TRows *= 5;
   }
   int TSize = TRows*n;
   int numberOfTables = 64 / t;
   printf("TRows: %d\n# of tables: %d\n", TRows, numberOfTables);
   long unsigned int* T = malloc((TSize*(numberOfTables)+n)*3*sizeof(long int));
   
   long unsigned int bits0, bits1, bits2, index, index0, index1, index2;
   long unsigned int c0, c1, c2, t0, t1, t2, s0, s1, s2, s3, rMod5Fix;
   int counter = 0;
   
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
      unsigned long long compensate = rdtsc();
      createTable(T, B, TRows, n, t, numberOfTables, i, 0, 0);
//      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
      startTime += (rdtsc() - compensate);
      for (int j = 0; j < m; j++){
         bits0 = A[(j*la + i)*3 + 0];// read 64 sliced elements out of A
         bits1 = A[(j*la + i)*3 + 1];
         bits2 = A[(j*la + i)*3 + 2];
//         printGF5(bits0, bits1, bits2);
         for (int k = numberOfTables - 1; k >= 0; k--){
            
            index0 = bits0 >> k*2;// read two elements out of A
            index1 = bits1 >> k*2;
            index2 = bits2 >> k*2;
//            printGF5(index0, index1, index2);
            index = ((index2 & 2) << 3) + ((index2 & 2) << 1) + ((index1 & 2) << 2) + (index1 & 2) + ((index0 & 2) << 1) + ((index0 & 2) >> 1) + ((index2 & 1) << 2) + ((index1 & 1) << 1) + (index0 & 1); // convert the bits into indexes
//                        printf("id: %lu\t%d\tk: %d\n", index, counter++, k);
            for (int l = 0; l < n; l++){
               t0 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 0];
               t1 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 1];
               t2 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 2];
//               printf("t: ");printGF5(t0, t1, t2);
               c0 = C[(j*n + l)*3 + 0];
               c1 = C[(j*n + l)*3 + 1];
               c2 = C[(j*n + l)*3 + 2];
//               printf("c: ");printGF5(c0, c1, c2);
               
               s0 = t0 ^ c0;// bit sliced addition
               s1 = t1 ^ c1;
               s2 = t2 ^ c2;
               s2 ^= s1 & (t0 & c0);
               //   s3 ^= s2 & (a1 & b1);
               s1 ^= t0 & c0;
               s2 ^= t1 & c1;
               s3 = t2 & c2;
               
               C[(j*n + l)*3 + 2] = (s2 & ~s1 & ~s0); //my meathod for modulus
               C[(j*n + l)*3 + 1] = (s3) | (s1 & s0) | (~s2 & s1);
               C[(j*n + l)*3 + 0] = (~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2);
//               printf("r: ");printGF5((~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2), (s3) | (s1 & s0) | (~s2 & s1), (s2 & ~s1 & ~s0));
            }
         }
      }
   }
   unsigned long long endTime = rdtsc();
      printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
//   printf("C\n");
//   printGF5Matrix(C, m, n);
   
   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   
//   printf("Af\n");
   convertGF5MatrixToFullSize(Af, A, m, lb, la);
//   printf("Bf\n");
   convertGF5MatrixToFullSize(Bf, B, lb, nActual, n);
   
   startTime = rdtsc();
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            Cf[i * nActual + j] += (Af[i * lb + k] * Bf[j + k * nActual]);
         }
         Cf[i * nActual + j] %= 5;
      }
   }
   endTime = rdtsc();
   //   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
//   printf("Cf\n"); // print Cf
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
   
   free(A);
   free(B);
   free(C);
   free(T);
}



