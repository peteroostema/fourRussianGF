/* gf5para 1, 2, 4.c implement a parallel FF MM for gf5 with
 * t = 1, 2, 4
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 * Revised: December 2017, for Calvin College CS 374 honors
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>

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

void printGF5Matrix256(long unsigned int* matrix, int rows, int columns);

static inline void createTable(long unsigned int* T, long unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){
      for (int j = 0; j < n; j++){
         *(T + (i*TRows*n + j)*3 + 0) = 0;
         *(T + (i*TRows*n + j)*3 + 1) = 0;
         *(T + (i*TRows*n + j)*3 + 2) = 0;
      }
   }
   int index = 0;
   
//   long unsigned int* test = malloc(6*4*3*sizeof(long int));
   
   for (int i = 0; i < numOfTables; i++){
      __m256i b00, b01, b02, b10, b11, b12, b20, b21, b22, b30, b31, b32, t0, t1, t2, s0, s1, s2, s3, r0, r1, r2;
      //      printf("subtable: %d,   index: %d\n", i, index);
      
      for (int o = 0; o < nq256; o++){
         index = o*3*4 + i*TRows*n*3;
//                  printf("o: %d\tindex: %d\n", o, index);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + o*4)*3 + 0), b30);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + o*4)*3 + 4), b31);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + o*4)*3 + 8), b32);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 1*n + o*4)*3 + 0), b20);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 1*n + o*4)*3 + 4), b21);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 1*n + o*4)*3 + 8), b22);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 2*n + o*4)*3 + 0), b10);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 2*n + o*4)*3 + 4), b11);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 2*n + o*4)*3 + 8), b12);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 3*n + o*4)*3 + 0), b00);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 3*n + o*4)*3 + 4), b01);
         VLOAD((Bt + (tableNumber*t*n*numOfTables + i*t*n + 3*n + o*4)*3 + 8), b02);
         //         b00 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 0];
         //         b01 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 1];
         //         b02 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 2];
         //         b10 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 0];
         //         b11 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 1];
         //         b12 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 2];
         
//                  _mm256_storeu_si256(test, b00);
//                  _mm256_storeu_si256(test + 4, b01);
//                  _mm256_storeu_si256(test + 8, b02);
//                  _mm256_storeu_si256(test + 12, b10);
//                  _mm256_storeu_si256(test + 16, b11);
//                  _mm256_storeu_si256(test + 20, b12);
//                  printf("BLoad\n");
//                  printGF5Matrix256(test, 2, 1);
         for (int j = 0; j < 5; j++){
            for (int k = 0; k < 5; k++){
               for (int l = 0; l < 5; l++){
                  for (int m = 0; m < 5; m++){
                     VLOAD((T + index + 0), t0);
                     VLOAD((T + index + 4), t1);
                     VLOAD((T + index + 8), t2);
                     //                  _mm256_storeu_si256(test, t0);
                     //                  _mm256_storeu_si256(test + 4, t1);
                     //                  _mm256_storeu_si256(test + 8, t2);
                     //                  printf("TableLoad\n");
                     //                  printGF5Matrix256(test, 1, 4);
                     index += 3*n;
                     
                     s0 = _mm256_xor_si256 (b00, t0);
                     s1 = _mm256_xor_si256 (b01, t1);
                     s2 = _mm256_xor_si256 (b02, t2);
                     s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(b00, t0), s1), s2);
                     s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(b00, t0));
                     s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(b01, t1));
                     s3 = _mm256_and_si256(b02, t2);
                     r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
                     r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
                     r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
                     
                     //                  _mm256_storeu_si256(test, s0);
                     //                  _mm256_storeu_si256(test + 4, s1);
                     //                  _mm256_storeu_si256(test + 8, s2);
                     //                  printf("AddTest\n");
                     //                  printGF5Matrix256(test, 1, 4);
                     
                     _mm256_storeu_si256((T + index + 0), r0);
                     _mm256_storeu_si256((T + index + 4), r1);
                     _mm256_storeu_si256((T + index + 8), r2);
                  }
//                  VLOAD((T + index + 0), t0);
//                  VLOAD((T + index + 4), t1);
//                  VLOAD((T + index + 8), t2);
                  //                  _mm256_storeu_si256(test, t0);
                  //                  _mm256_storeu_si256(test + 4, t1);
                  //                  _mm256_storeu_si256(test + 8, t2);
                  //                  printf("TableLoad\n");
                  //                  printGF5Matrix256(test, 1, 4);
                  
                  
                  s0 = _mm256_xor_si256 (b10, r0);
                  s1 = _mm256_xor_si256 (b11, r1);
                  s2 = _mm256_xor_si256 (b12, r2);
                  s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(b10, r0), s1), s2);
                  s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(b10, r0));
                  s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(b11, r1));
                  s3 = _mm256_and_si256(b12, r2);
                  r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
                  r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
                  r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
                  
                  //                  _mm256_storeu_si256(test, s0);
                  //                  _mm256_storeu_si256(test + 4, s1);
                  //                  _mm256_storeu_si256(test + 8, s2);
                  //                  printf("AddTest\n");
                  //                  printGF5Matrix256(test, 1, 4);
                  
                  _mm256_storeu_si256((T + index + 0), r0);
                  _mm256_storeu_si256((T + index + 4), r1);
                  _mm256_storeu_si256((T + index + 8), r2);
               }
//               VLOAD((T + index + 0), t0);
//               VLOAD((T + index + 4), t1);
//               VLOAD((T + index + 8), t2);
               //                  _mm256_storeu_si256(test, t0);
               //                  _mm256_storeu_si256(test + 4, t1);
               //                  _mm256_storeu_si256(test + 8, t2);
               //                  printf("TableLoad\n");
               //                  printGF5Matrix256(test, 1, 4);
               
               
               s0 = _mm256_xor_si256 (b20, r0);
               s1 = _mm256_xor_si256 (b21, r1);
               s2 = _mm256_xor_si256 (b22, r2);
               s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(b20, r0), s1), s2);
               s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(b20, r0));
               s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(b21, r1));
               s3 = _mm256_and_si256(b22, r2);
               r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
               r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
               r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
               
               //                  _mm256_storeu_si256(test, s0);
               //                  _mm256_storeu_si256(test + 4, s1);
               //                  _mm256_storeu_si256(test + 8, s2);
               //                  printf("AddTest\n");
               //                  printGF5Matrix256(test, 1, 4);
               
               _mm256_storeu_si256((T + index + 0), r0);
               _mm256_storeu_si256((T + index + 4), r1);
               _mm256_storeu_si256((T + index + 8), r2);
            }
//            VLOAD((T + index + 0), t0);
//            VLOAD((T + index + 4), t1);
//            VLOAD((T + index + 8), t2);
            //                  _mm256_storeu_si256(test, t0);
            //                  _mm256_storeu_si256(test + 4, t1);
            //                  _mm256_storeu_si256(test + 8, t2);
            //                  printf("TableLoad\n");
            //                  printGF5Matrix256(test, 1, 4);
            
            
            s0 = _mm256_xor_si256 (b30, r0);
            s1 = _mm256_xor_si256 (b31, r1);
            s2 = _mm256_xor_si256 (b32, r2);
            s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(b30, r0), s1), s2);
            s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(b30, r0));
            s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(b31, r1));
            s3 = _mm256_and_si256(b32, r2);
            r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
            r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
            r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
            
            //                  _mm256_storeu_si256(test, s0);
            //                  _mm256_storeu_si256(test + 4, s1);
            //                  _mm256_storeu_si256(test + 8, s2);
            //                  printf("AddTest\n");
            //                  printGF5Matrix256(test, 1, 4);
            
            _mm256_storeu_si256((T + index + 0), r0);
            _mm256_storeu_si256((T + index + 4), r1);
            _mm256_storeu_si256((T + index + 8), r2);
         }
      }
   }
   //   printf("\nfinal index: %d\n", index);
}

void print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   int subtableSize = rows*(majorColumns*innerColumns + outerColumns)*4;
   for (int l = 0; l < depth; l++){
      printf("subTable: %d\nsubtableSize: %d\n", l, subtableSize);
      
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < innerColumns; j++){
            for (int k = 0; k < 4; k++){
               //            tmp = matrix[l*subtableSize + j + i*innerColumns];
               //            for (int k = 15; k >= 0; k--){
               //               long unsigned tmp2 = (tmp >> k*4) & 0xf;
               //               //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
               //               printf("%lu", tmp2); //"   %lu"
               //            }
               //            tmp0 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 0];
               //            tmp1 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 1];
               //            tmp2 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 2];
               tmp0 = matrix[(l*subtableSize + i*innerColumns*4 + j*4)*3 + k + 0];
               tmp1 = matrix[(l*subtableSize + i*innerColumns*4 + j*4)*3 + k + 4];
               tmp2 = matrix[(l*subtableSize + i*innerColumns*4 + j*4)*3 + k + 8];
               for (int i = 63; i >= 0; i--){
                  long unsigned print0 = ((tmp0 >> i) & 1);
                  long unsigned print1 = ((tmp1 >> i) & 1);
                  long unsigned print2 = ((tmp2 >> i) & 1);
                  long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
                  //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
                  printf("%lu", printInt); //"   %lu"
               }
               
               //         printf("\n");
               //         for (int k = 0; k < 64; k++){
               //            if (tmp & 0x8000000000000000){
               //               printf("1");
               //            } else {
               //               printf("0");
               //            }
               //            tmp <<= 1;
               //         }
               printf("\t");
            }
         }
         printf(" %d", rowCounter);
         rowCounter++;
         printf("\n");
      }
      
      rowCounter = 0;
   }
   printf("\n");
}

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
            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
            printf("%lu", printInt); //"   %lu"
         }
         //         printf("\n");
         //         for (int k = 0; k < 64; k++){
         //            if (tmp & 0x8000000000000000){
         //               printf("1");
         //            } else {
         //               printf("0");
         //            }
         //            tmp <<= 1;
         //         }
         printf("\t");
      }
      printf(" %d", rowCounter);
      rowCounter++;
      printf("\n");
   }
   printf("\n");
}

//void printGF5Matrix256(long unsigned int* matrix, int rows, int columns){
//   long unsigned int tmp0, tmp1, tmp2;
//   int rowCounter = 0;
//   for (int i = 0; i < rows; i++){
//      for (int j = 0; j < columns; j++){
//         tmp0 = matrix[(i*columns)*3 + j + 0];
//         tmp1 = matrix[(i*columns)*3 + j + 4];
//         tmp2 = matrix[(i*columns)*3 + j + 8];
//         for (int i = 63; i >= 0; i--){
//            long unsigned print0 = ((tmp0 >> i) & 1);
//            long unsigned print1 = ((tmp1 >> i) & 1);
//            long unsigned print2 = ((tmp2 >> i) & 1);
//            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
//            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
//            printf("%lu", printInt); //"   %lu"
//         }
//         //         printf("\n");
//         //         for (int k = 0; k < 64; k++){
//         //            if (tmp & 0x8000000000000000){
//         //               printf("1");
//         //            } else {
//         //               printf("0");
//         //            }
//         //            tmp <<= 1;
//         //         }
//         printf("\t");
//      }
//      printf(" %d", rowCounter);
//      rowCounter++;
//      printf("\n");
//   }
//   printf("\n");
//}

void printGF5Matrix256(long unsigned int* matrix, int rows, int columns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
//   printf("columns: %d\n", columns);
//   printf("column: %d\n", columns%4);
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         for (int k = 0; k < 4; k++){
            tmp0 = matrix[(i*columns*4 + j*4)*3 + k + 0];
            tmp1 = matrix[(i*columns*4 + j*4)*3 + k + 4];
            tmp2 = matrix[(i*columns*4 + j*4)*3 + k + 8];
            for (int i = 63; i >= 0; i--){
               long unsigned print0 = ((tmp0 >> i) & 1);
               long unsigned print1 = ((tmp1 >> i) & 1);
               long unsigned print2 = ((tmp2 >> i) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
               printf("%lu", printInt); //"   %lu"
            }
            //         printf("\n");
            //         for (int k = 0; k < 64; k++){
            //            if (tmp & 0x8000000000000000){
            //               printf("1");
            //            } else {
            //               printf("0");
            //            }
            //            tmp <<= 1;
            //         }
            printf("\t");
         }
      }
//      for (int j = 0; j < columns%4; j++){
//         tmp0 = matrix[(columns + i*columns)*3 + j + 0];
//         tmp1 = matrix[(columns + i*columns)*3 + j + 4];
//         tmp2 = matrix[(columns + i*columns)*3 + j + 8];
//         for (int i = 63; i >= 0; i--){
//            long unsigned print0 = ((tmp0 >> i) & 1);
//            long unsigned print1 = ((tmp1 >> i) & 1);
//            long unsigned print2 = ((tmp2 >> i) & 1);
//            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
//            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
//            printf("%lu", printInt); //"   %lu"
//         }
//         //         printf("\n");
//         //         for (int k = 0; k < 64; k++){
//         //            if (tmp & 0x8000000000000000){
//         //               printf("1");
//         //            } else {
//         //               printf("0");
//         //            }
//         //            tmp <<= 1;
//         //         }
//         printf("\t");
//      }
      printf(" %d", rowCounter);
      rowCounter++;
      printf("\n");
   }
   printf("\n");
}

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
//if (j < la - 1){
//   for (int i = 0; i < 64; i++){
//      long unsigned int tmp = rand() % 5;
//      //               printf("%lu", tmp);
//      a2 <<= 1;
//      a1 <<= 1;
//      a0 <<= 1;
//      if (tmp & 0b100){
//         a2 |= 1;
//      }
//      if (tmp & 0b010){
//         a1 |= 1;
//      }
//      if (tmp & 0b001){
//         a0 |= 1;
//      }
//   }
//} else{
//   for (int i = 0; i < 64 - lPadding; i++){
//      long unsigned int tmp = rand() % 5;
void convertGF5MatrixToFullSize(unsigned int* fullSizeMatrix, long unsigned int* gf5Matrix, int rows, int columns, int integerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp0 = gf5Matrix[(i*integerColumns + j)*3 + 0];
         tmp1 = gf5Matrix[(i*integerColumns + j)*3 + 1];
         tmp2 = gf5Matrix[(i*integerColumns + j)*3 + 2];
         if (j < integerColumns - 1){
            for (int k = 63; k >= 0; k--){
               long unsigned print0 = ((tmp0 >> k) & 1);
               long unsigned print1 = ((tmp1 >> k) & 1);
               long unsigned print2 = ((tmp2 >> k) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
            }
         } else {
            for (int k = 63; k >= (64 - (columns % 64)) % 64; k--){
               long unsigned print0 = ((tmp0 >> k) & 1);
               long unsigned print1 = ((tmp1 >> k) & 1);
               long unsigned print2 = ((tmp2 >> k) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
            }
         }
      }
      //rowCounter++;
   }
//         printf("size: %d\n", rows*columns);
//         for (int i = 0; i < rows; i++){
//            for (int j = 0; j < columns; j++){
//               printf("%d", fullSizeMatrix[i*columns + j]);
//            }
//            printf("\t");
//            printf(" %d", rowCounter);
//            printf("\n");
//            rowCounter++;
//         }
//         printf("\n");
}
void convertGF5MatrixToFullSize256(unsigned int* fullSizeMatrix, long unsigned int* gf5Matrix, int rows, int columns, int integerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         for (int l = 0; l < 4; l++){
            tmp0 = gf5Matrix[(i*integerColumns*4 + j*4)*3 + l + 0];
            tmp1 = gf5Matrix[(i*integerColumns*4 + j*4)*3 + l + 4];
            tmp2 = gf5Matrix[(i*integerColumns*4 + j*4)*3 + l + 8];
            for (int k = 63; k >= 0; k--){
               long unsigned print0 = ((tmp0 >> k) & 1);
               long unsigned print1 = ((tmp1 >> k) & 1);
               long unsigned print2 = ((tmp2 >> k) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               fullSizeMatrix[i*columns + j*256 + l*64 + (63 - k)] = printInt;
            }
         }
      }
      //rowCounter++;
   }
//         printf("size: %d\n", rows*columns);
//         for (int i = 0; i < rows; i++){
//            for (int j = 0; j < columns; j++){
//               printf("%d", fullSizeMatrix[i*columns + j]);
//            }
//            printf("\t");
//            printf(" %d", rowCounter);
//            printf("\n");
//            rowCounter++;
//         }
//         printf("\n");
}
//for (int i = 0; i < rows; i++){
//   for (int j = 0; j < columns/4; j++){
//      for (int k = 0; k < 4; k++){
//         tmp0 = matrix[(i*columns + j*4)*3 + k + 0];
//         tmp1 = matrix[(i*columns + j*4)*3 + k + 4];
//         tmp2 = matrix[(i*columns + j*4)*3 + k + 8];
int checkCorrectResult(long unsigned int* C, unsigned int* Cf, int rows, int columns, int integerColumns){
   long unsigned int tmp = 0;
   long unsigned int tmp0, tmp1, tmp2;
   //   for (int i = 0; i < rows; i++){
   //      for (int j = 0; j < integerColumns; j++){
   //         tmp = C[i*integerColumns + j];
   //         for (int k = 15; k >= 0; k--){
   //            if (j*16 + (16 - k) >= columns){
   //               break;
   //            }
   //            tmp2 = (tmp >> k*4) & 0xf;
   //            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
   //            if (!(tmp2 == Cf[i*columns + j * 16 + (15 - k)])){
   //               return 0;
   //            }
   //         }
   //      }
   //   }
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         for (int l = 0; l < 4; l++){
            tmp0 = C[(i*integerColumns*4 + j*4)*3 + l + 0];
            tmp1 = C[(i*integerColumns*4 + j*4)*3 + l + 4];
            tmp2 = C[(i*integerColumns*4 + j*4)*3 + l + 8];
            for (int k = 63; k >= 0; k--){
               long unsigned print0 = ((tmp0 >> k) & 1);
               long unsigned print1 = ((tmp1 >> k) & 1);
               long unsigned print2 = ((tmp2 >> k) & 1);
               long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
               //            fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
               
               if (j*256 + 64*l + (64 - k) >= columns){
                  break;
               }
               //            tmp2 = (tmp >> k*4) & 0xf;
               //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
               if (!(printInt == Cf[i*columns + j*256 + l*64 + (63 - k)])){
                  printf("row: %d\ncolumn: %d\nk: %d\nl: %d\n", i, j, k, l);
                  printf("C: %lu\n", printInt);
                  printf("Cf: %u\n", Cf[i*columns + j*256 + l*64 + (63 - k)]);
                  return 0;
               }
            }
         }
      }
      //rowCounter++;
   }
   return 1;
}

static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

int main(int argc, char *argv[]){
   srand(time(NULL));
   int m = atoi(argv[1]);  // rows in A
   int lb = atoi(argv[2]); // rows in B
   int la = (lb - 1) / 64 + 1;  // integer columns in A
   int nActual = atoi(argv[3]);
   int nThreads = atoi(argv[4]);
   int n = (nActual - 1) / 64 + 1;  // integer columns in B
   int nLength = n + (4 - n%4);
   if (!(n%4)){
      nLength = n;
   }
   //   int ra = 64 - (lb % 64); // zeros to pad on A
   //   int rb = (64 - (nActual % 64)) % 64; // zeros to pad on B
   int t = 4;  // four Russian divisor
   
   printf("m: %d\n", m);
   printf("lb: %d\n", lb);
   printf("la: %d\n", la);
   printf("n: %d\n", n);
   printf("t: %d\n", t);
   //   printf("ra: %d\n", ra);
   //   printf("rb: %d\n", rb);
   int lPadding = (64 - (lb % 64)) % 64;
   int nPadding = (64 - (nActual % 64)) % 64;
   int mLength = m + m%4;
   printf("Lpadding: %d\n", lPadding);
   printf("Npadding: %d\n", nPadding);
   printf("nLength: %d\n", nLength);
   
   int TRows = 1;
   for (int i = 0; i < t; i++){
      TRows *= 5;
   }
   int TSize = TRows*nLength;
   int numberOfTables = 64 / t;
   int nq256 = (nActual - 1) / 256 + 1;
   int nr256 = n%4;//((n - 1) % 4) + 1;
   if (nr256 == 0){//   if (!(nActual % 256)){
      nr256 = 4;//      nr256 = 0;
   }//   }
   
   long unsigned int* A = malloc(mLength*la*3*(sizeof(long unsigned int)));
   long unsigned int* B = malloc((lb + lPadding)*nq256*4*3*(sizeof(long unsigned int)));
   long unsigned int* C = calloc(m*nq256*4*3, (sizeof(long unsigned int)));
   
   for (int i = 0; i < m; i++){
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
//      printf("A\n");
//      printGF5Matrix(A, m, la);
   
   for (int i = 0; i < lb; i++){
      for (int j = 0; j < nq256; j++){
         for (int k = 0; k < 4; k++){
            long unsigned int b2 = 0, b1 = 0, b0 = 0;
            if ( (j == nq256 - 1) && (k == nr256)){
               break;
            }
            if (j*4 + k < n - 1){
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
            B[(i*nq256*4 + j*4)*3 + k + 8] = b2;
            B[(i*nq256*4 + j*4)*3 + k + 4] = b1;
            B[(i*nq256*4 + j*4)*3 + k + 0] = b0;
         }
      }
   }
//      printf("B\n");
//      printGF5Matrix256(B, lb, nq256);
   //   printGF5(B[0], B[4], B[8]);
   //   printGF5(B[1], B[5], B[9]);
   
   
   printf("TRows: %d\n# of tables: %d\n", TRows, numberOfTables);
   printf("nq256: %d\tnr256: %d\n", nq256, nr256);
   long unsigned int* T = malloc((TSize*(numberOfTables)+n)*3*sizeof(long int));
   //   createTable(T, B, TRows, n, t, numberOfTables, 0, 0, 0);
//      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
   //print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns)
   

   int counter = 0;
    long unsigned int* test = malloc(6*4*3*sizeof(long int));
   
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
//      unsigned long long compensate = rdtsc();
      createTable(T, B, TRows, nLength, t, numberOfTables, i, nq256, nr256);
//            print4Dgf5Matrix(T, TRows, 1, numberOfTables, nq256, 0);
//      startTime += (rdtsc() - compensate);
      //      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
#pragma omp parallel for num_threads(nThreads)
      for (int j = 0; j < m; j+=4){
   long unsigned int bits0, bits1, bits2, bits00, bits01, bits02, bits10, bits11, bits12, bits20, bits21, bits22, bits30, bits31, bits32, index0, index00, index01, index02, index1, index10, index11, index12, index2, index20, index21, index22, index3, index30, index31, index32;
   __m256i c00, c01, c02, c10, c11, c12, c20, c21, c22, c30, c31, c32, t00, t01, t02, t10, t11, t12, t20, t21, t22, t30, t31, t32, s0, s1, s2, s3, r0, r1, r2;
//         compensate = rdtsc();
//         bits0 = A[(j*la + i + 0)*3 + 0];
//         bits1 = A[(j*la + i + 0)*3 + 1];
//         bits2 = A[(j*la + i + 0)*3 + 2];
         bits00 = A[((j + 0)*la + i)*3 + 0];
         bits01 = A[((j + 0)*la + i)*3 + 1];
         bits02 = A[((j + 0)*la + i)*3 + 2];
         bits10 = A[((j + 1)*la + i)*3 + 0];
         bits11 = A[((j + 1)*la + i)*3 + 1];
         bits12 = A[((j + 1)*la + i)*3 + 2];
         bits20 = A[((j + 2)*la + i)*3 + 0];
         bits21 = A[((j + 2)*la + i)*3 + 1];
         bits22 = A[((j + 2)*la + i)*3 + 2];
         bits30 = A[((j + 3)*la + i)*3 + 0];
         bits31 = A[((j + 3)*la + i)*3 + 1];
         bits32 = A[((j + 3)*la + i)*3 + 2];
//         startTime += (rdtsc() - compensate);
         //         bits00 = A[(j*la*4 + i + 0)*3 + 0];
         //         bits01 = A[(j*la*4 + i + 0)*3 + 1];
         //         bits02 = A[(j*la*4 + i + 0)*3 + 2];
         //         bits10 = A[(j*la*4 + i + 1)*3 + 0];
         //         bits11 = A[(j*la*4 + i + 1)*3 + 1];
         //         bits12 = A[(j*la*4 + i + 1)*3 + 2];
         //         bits20 = A[(j*la*4 + i + 2)*3 + 0];
         //         bits21 = A[(j*la*4 + i + 2)*3 + 1];
         //         bits22 = A[(j*la*4 + i + 2)*3 + 2];
         //         bits30 = A[(j*la*4 + i + 3)*3 + 0];
         //         bits31 = A[(j*la*4 + i + 3)*3 + 1];
         //         bits32 = A[(j*la*4 + i + 3)*3 + 2];
         //         printGF5(bits0, bits1, bits2);
         for (int k = numberOfTables - 1; k >= 0; k--){
//            compensate = rdtsc();
            
            index00 = bits00 >> k*4;
            index01 = bits01 >> k*4;
            index02 = bits02 >> k*4;
            index10 = bits10 >> k*4;
            index11 = bits11 >> k*4;
            index12 = bits12 >> k*4;
            index20 = bits20 >> k*4;
            index21 = bits21 >> k*4;
            index22 = bits22 >> k*4;
            index30 = bits30 >> k*4;
            index31 = bits31 >> k*4;
            index32 = bits32 >> k*4;
            
            
            //            printGF5(index0, index1, index2);
            index0 = (index02 & 8) * 62.5 + (index01 & 8) * 31.25 + (index00 & 8) * 15.625
            + (index02 & 4) * 25 + (index01 & 4) * 12.5 + (index00 & 4) * 6.25
            + ((index02 & 2) << 3) + ((index02 & 2) << 1) + ((index01 & 2) << 2) + (index01 & 2) + ((index00 & 2) << 1) + ((index00 & 2) >> 1)
            + ((index02 & 1) << 2) + ((index01 & 1) << 1) + (index00 & 1);
            index1 =(index12 & 8) * 62.5 + (index11 & 8) * 31.25 + (index10 & 8) * 15.625
            + (index12 & 4) * 25 + (index11 & 4) * 12.5 + (index10 & 4) * 6.25
            + ((index12 & 2) << 3) + ((index12 & 2) << 1) + ((index11 & 2) << 2) + (index11 & 2) + ((index10 & 2) << 1) + ((index10 & 2) >> 1)
            + ((index12 & 1) << 2) + ((index11 & 1) << 1) + (index10 & 1);
            index2 =(index22 & 8) * 62.5 + (index21 & 8) * 31.25 + (index20 & 8) * 15.625
            + (index22 & 4) * 25 + (index21 & 4) * 12.5 + (index20 & 4) * 6.25
            + ((index22 & 2) << 3) + ((index22 & 2) << 1) + ((index21 & 2) << 2) + (index21 & 2) + ((index20 & 2) << 1) + ((index20 & 2) >> 1)
            + ((index22 & 1) << 2) + ((index21 & 1) << 1) + (index20 & 1);
            index3 =(index32 & 8) * 62.5 + (index31 & 8) * 31.25 + (index30 & 8) * 15.625
            + (index32 & 4) * 25 + (index31 & 4) * 12.5 + (index30 & 4) * 6.25
            + ((index32 & 2) << 3) + ((index32 & 2) << 1) + ((index31 & 2) << 2) + (index31 & 2) + ((index30 & 2) << 1) + ((index30 & 2) >> 1)
            + ((index32 & 1) << 2) + ((index31 & 1) << 1) + (index30 & 1);
//               printf("id: %lu, %lu, %lu, %lu\t%d\n", index0, index1, index2, index3, counter++);
//            startTime += (rdtsc() - compensate);
            for (int l = 0; l < nq256; l++){
               
               VLOAD((T + (index0*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 0), t00);
               VLOAD((T + (index0*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 4), t01);
               VLOAD((T + (index0*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 8), t02);
               VLOAD((T + (index1*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 0), t10);
               VLOAD((T + (index1*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 4), t11);
               VLOAD((T + (index1*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 8), t12);
               VLOAD((T + (index2*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 0), t20);
               VLOAD((T + (index2*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 4), t21);
               VLOAD((T + (index2*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 8), t22);
               VLOAD((T + (index3*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 0), t30);
               VLOAD((T + (index3*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 4), t31);
               VLOAD((T + (index3*nLength + TSize*(numberOfTables - k - 1) + l*4)*3 + 8), t32);
               VLOAD((C + ((j + 0)*nLength + l*4)*3 + 0), c00);
               VLOAD((C + ((j + 0)*nLength + l*4)*3 + 4), c01);
               VLOAD((C + ((j + 0)*nLength + l*4)*3 + 8), c02);
               VLOAD((C + ((j + 1)*nLength + l*4)*3 + 0), c10);
               VLOAD((C + ((j + 1)*nLength + l*4)*3 + 4), c11);
               VLOAD((C + ((j + 1)*nLength + l*4)*3 + 8), c12);
               VLOAD((C + ((j + 2)*nLength + l*4)*3 + 0), c20);
               VLOAD((C + ((j + 2)*nLength + l*4)*3 + 4), c21);
               VLOAD((C + ((j + 2)*nLength + l*4)*3 + 8), c22);
               VLOAD((C + ((j + 3)*nLength + l*4)*3 + 0), c30);
               VLOAD((C + ((j + 3)*nLength + l*4)*3 + 4), c31);
               VLOAD((C + ((j + 3)*nLength + l*4)*3 + 8), c32);
               
               s0 = _mm256_xor_si256 (c00, t00);
               s1 = _mm256_xor_si256 (c01, t01);
               s2 = _mm256_xor_si256 (c02, t02);
               s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(c00, t00), s1), s2);
               s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(c00, t00));
               s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(c01, t01));
               s3 = _mm256_and_si256(c02, t02);
               r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
               r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
               r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
               
               _mm256_storeu_si256((C + (j*nLength + l*4)*3 + 0), r0);
               _mm256_storeu_si256((C + (j*nLength + l*4)*3 + 4), r1);
               _mm256_storeu_si256((C + (j*nLength + l*4)*3 + 8), r2);
               
//               _mm256_storeu_si256(test, t00);
//               _mm256_storeu_si256(test + 4, t01);
//               _mm256_storeu_si256(test + 8, t02);
//               _mm256_storeu_si256(test + 12, r0);
//               _mm256_storeu_si256(test + 16, r1);
//               _mm256_storeu_si256(test + 20, r2);
//               printf("TableLoad\n");
//               printGF5Matrix256(test, 2, 1);
               
               s0 = _mm256_xor_si256 (c10, t10);
               s1 = _mm256_xor_si256 (c11, t11);
               s2 = _mm256_xor_si256 (c12, t12);
               s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(c10, t10), s1), s2);
               s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(c10, t10));
               s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(c11, t11));
               s3 = _mm256_and_si256(c12, t12);
               r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
               r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
               r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
               
               _mm256_storeu_si256((C + ((j + 1)*nLength + l*4)*3 + 0), r0);
               _mm256_storeu_si256((C + ((j + 1)*nLength + l*4)*3 + 4), r1);
               _mm256_storeu_si256((C + ((j + 1)*nLength + l*4)*3 + 8), r2);
               
//               _mm256_storeu_si256(test, t10);
//               _mm256_storeu_si256(test + 4, t11);
//               _mm256_storeu_si256(test + 8, t12);
//               _mm256_storeu_si256(test + 12, r0);
//               _mm256_storeu_si256(test + 16, r1);
//               _mm256_storeu_si256(test + 20, r2);
//               printf("TableLoad\n");
//               printGF5Matrix256(test, 2, 1);
               
               s0 = _mm256_xor_si256 (c20, t20);
               s1 = _mm256_xor_si256 (c21, t21);
               s2 = _mm256_xor_si256 (c22, t22);
               s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(c20, t20), s1), s2);
               s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(c20, t20));
               s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(c21, t21));
               s3 = _mm256_and_si256(c22, t22);
               r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
               r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
               r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
               
               _mm256_storeu_si256((C + ((j + 2)*nLength + l*4)*3 + 0), r0);
               _mm256_storeu_si256((C + ((j + 2)*nLength + l*4)*3 + 4), r1);
               _mm256_storeu_si256((C + ((j + 2)*nLength + l*4)*3 + 8), r2);
               
               s0 = _mm256_xor_si256 (c30, t30);
               s1 = _mm256_xor_si256 (c31, t31);
               s2 = _mm256_xor_si256 (c32, t32);
               s2 = _mm256_xor_si256 ( _mm256_and_si256( _mm256_and_si256(c30, t30), s1), s2);
               s1 = _mm256_xor_si256 ( s1, _mm256_and_si256(c30, t30));
               s2 = _mm256_xor_si256 ( s2, _mm256_and_si256(c31, t31));
               s3 = _mm256_and_si256(c32, t32);
               r2 = _mm256_andnot_si256 (s0, _mm256_andnot_si256(s1, s2));
               r1 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_and_si256(s0, s1), _mm256_andnot_si256(s2, s1)), s3);
               r0 = _mm256_or_si256 ( _mm256_or_si256 ( _mm256_andnot_si256(s2, s0), _mm256_andnot_si256(s2, s3)), _mm256_and_si256(_mm256_andnot_si256(s0, s1), s2));
               
               _mm256_storeu_si256((C + ((j + 3)*nLength + l*4)*3 + 0), r0);
               _mm256_storeu_si256((C + ((j + 3)*nLength + l*4)*3 + 4), r1);
               _mm256_storeu_si256((C + ((j + 3)*nLength + l*4)*3 + 8), r2);
               
               //               t0 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 0];
               //               t1 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 1];
               //               t2 = T[(index*n + TSize*(numberOfTables - k - 1) + l)*3 + 2];
               ////               printf("t: ");printGF5(t0, t1, t2);
               //               c0 = C[(j*n + l)*3 + 0];
               //               c1 = C[(j*n + l)*3 + 1];
               //               c2 = C[(j*n + l)*3 + 2];
               ////               printf("c: ");printGF5(c0, c1, c2);
               //               s0 = t0 ^ c0;// bit packed addition
               //               s1 = t1 ^ c1;
               //               s2 = t2 ^ c2;
               //               s2 ^= s1 & (t0 & c0);
               //               //   s3 ^= s2 & (a1 & b1);
               //               s1 ^= t0 & c0;
               //               s2 ^= t1 & c1;
               //               s3 = t2 & c2;
               //               C[(j*n + l)*3 + 2] = (s2 & ~s1 & ~s0); //my meathod
               //               C[(j*n + l)*3 + 1] = (s3) | (s1 & s0) | (~s2 & s1);
               //               C[(j*n + l)*3 + 0] = (~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2);
               //               printf("r: ");printGF5((~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2), (s3) | (s1 & s0) | (~s2 & s1), (s2 & ~s1 & ~s0));
               
               //               VLOAD(((index*n + TSize*(numberOfTables - k - 1) + l)*3 + 0), t0);
               //               VLOAD((subDiv + index01), xmm1);
               //               VLOAD((subDiv + index02), xmm2);
            }
         }
         //                     printf("C\n");
         //                     printGF5Matrix(C, m, n);
      }
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
//      printf("C\n");
//      printGF5Matrix256(C, m, nq256);
   
/*   ///////////////////////////////////////////////////
   // testing                                        /
   ///////////////////////////////////////////////////
   
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   
//      printf("Af\n");
   convertGF5MatrixToFullSize(Af, A, m, lb, la);
//      printf("Bf\n");
   convertGF5MatrixToFullSize256(Bf, B, lb, nActual, nq256);
   
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
   
   printf("Correct = 1, Wrong = 0: %d\n", checkCorrectResult(C, Cf, m, nActual, nq256));*/
   
   free(A);
   free(B);
   free(C);
   free(T);
}



