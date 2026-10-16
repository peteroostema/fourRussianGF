/* slicedt2gf7.c is a modification of t2fourrusgf5.c to work
 * with GF7 addition
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void printBinary(long unsigned int number);
void printGF5(long unsigned int r0, long unsigned int r1, long unsigned int r2);
void printGF5Added(long unsigned int r0, long unsigned int r1, long unsigned int r2, long unsigned int r3);
static inline void createTable(long unsigned int* T, long unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){
      for (int j = 0; j < n; j++){
         *(T + (i*TRows + j)*3 + 0) = 0;
         *(T + (i*TRows + j)*3 + 1) = 0;
         *(T + (i*TRows + j)*3 + 2) = 0;
      }
   }
   int index = 0;
   
   for (int i = 0; i < numOfTables; i++){
      long unsigned int b00, b01, b02, b10, b11, b12, t0, t1, t2, s0, s1, s2, s3;
//      printf("subtable: %d,   index: %d\n", i, index);
      
      for (int o = 0; o < n; o++){
         index = o*3 + i*TRows*n*3;
         b00 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 0];
         b01 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 1];
         b02 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 2];
         b10 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 0];
         b11 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 1];
         b12 = Bt[(tableNumber*t*n*numOfTables + i*t*n + 0 + o)*3 + 2];
         
//         printf("index: %d\n", (tableNumber*t*n*numOfTables + i*t*n + 1*n + o)*3 + 0);
//         printf("b0\n");
//         printGF5(b00, b01, b02);
//         printBinary(b02);
//         printBinary(b01);
//         printBinary(b00);
         for (int l = 0; l < 7; l++){
            for (int m = 0; m < 7; m++){
//               printf("index: %d\n", index);
               
//               s = T[index] + b0;
//               index += n;
//               s0 = s & 0x1111111111111111;
//               s1 = (s & 0x2222222222222222) >> 1;
//               s2 = (s & 0x4444444444444444) >> 2;
//               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
//               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//               T[index] = (s - rMod5Fix);
               
               t0 = T[(index + 0)];
               t1 = T[(index + 1)];
               t2 = T[(index + 2)];
//               printGF5(t0, t1, t2);
               index += 3*n;
               s0 = t0 ^ b00;// bit sliced addition
               s1 = t1 ^ b01;
               s2 = t2 ^ b02;
               s3 = (s2 & ((t1 & b01) | (s1 & (t0 & b00)))) ^ (t2 & b02);
               s2 ^= (s1 & (t0 & b00)) ^ t1 & b01;
               s1 ^= t0 & b00;
               T[(index + 2)] = (s2 ^ (s1 & (s0 & s3))) & ~(s0 & s1 & s2); //my meathod
               T[(index + 1)] = (s1 ^ (s0 & s3)) & ~(s0 & s1 & s2);
               T[(index + 0)] = (s0 ^ s3) & ~(s0 & s1 & s2);
               
//               printf("t\n");
//               printGF5(t0, t1, t2);
//               printf("b0\n");
//               printGF5(b00, b01, b02);
//               printf("r\n");
//               printGF5(s0 ^ s3, s1 ^ (s0 & s3),  s2 ^ (s1 & (s0 & s3)));
               
            }
            t0 = T[(index + 0)];
            t1 = T[(index + 1)];
            t2 = T[(index + 2)];
            s0 = t0 ^ b10;// bit sliced addition
            s1 = t1 ^ b11;
            s2 = t2 ^ b12;
            s3 = (s2 & ((t1 & b11) | (s1 & (t0 & b10)))) ^ (t2 & b12);
            s2 ^= (s1 & (t0 & b10)) ^ t1 & b11;
            s1 ^= t0 & b10;
            T[(index + 2)] = (s2 ^ (s1 & (s0 & s3))) & ~(s0 & s1 & s2); //my meathod
            T[(index + 1)] = (s1 ^ (s0 & s3)) & ~(s0 & s1 & s2);
            T[(index + 0)] = (s0 ^ s3) & ~(s0 & s1 & s2);
//            printf("index: %d\n", index);
//                           printf("t\n");
//                           printGF5(t0, t1, t2);
//                           printf("b1\n");
//                           printGF5(b10, b11, b12);
//                           printf("s\n");
//                           printGF5Added(s0, s1, s2, s3);
//                           printf("r\n");
//                           printGF5((s0 ^ s3) & ~(s0 & s1 & s2), (s1 ^ (s0 & s3)) & ~(s0 & s1 & s2),  (s2 ^ (s1 & (s0 & s3))) & ~(s0 & s1 & s2));
         }
      }
   }
//   printf("\nfinal index: %d\n", index);
}

void print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   int subtableSize = rows*(majorColumns*innerColumns + outerColumns);
   for (int l = 0; l < depth; l++){
      printf("subTable: %d\nsubtableSize: %d\n", l, subtableSize);
      
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < innerColumns; j++){
//            tmp = matrix[l*subtableSize + j + i*innerColumns];
//            for (int k = 15; k >= 0; k--){
//               long unsigned tmp2 = (tmp >> k*4) & 0xf;
//               //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
//               printf("%lu", tmp2); //"   %lu"
//            }
            tmp0 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 0];
            tmp1 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 1];
            tmp2 = matrix[(l*subtableSize + i*innerColumns + j)*3 + 2];
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
void printGF5Added(long unsigned int r0, long unsigned int r1, long unsigned int r2, long unsigned int r3){
   long unsigned tmp0 = r0;
   long unsigned tmp1 = r1;
   long unsigned tmp2 = r2;
   long unsigned tmp3 = r3;
   for (int i = 63; i >= 0; i--){
      long unsigned print0 = ((tmp0 >> i) & 1);
      long unsigned print1 = ((tmp1 >> i) & 1);
      long unsigned print2 = ((tmp2 >> i) & 1);
      long unsigned print3 = ((tmp3 >> i) & 1);
      long unsigned printInt = (print3 << 3) | (print2 << 2) | (print1 << 1) | print0;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("%lu", printInt); //"   %lu"
   }
   printf("\n");
}

void convertGF5MatrixToFullSize(unsigned int* fullSizeMatrix, long unsigned int* gf5Matrix, int rows, int columns, int integerColumns){
   long unsigned int tmp0, tmp1, tmp2;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp0 = gf5Matrix[(i*integerColumns + j)*3 + 0];
         tmp1 = gf5Matrix[(i*integerColumns + j)*3 + 1];
         tmp2 = gf5Matrix[(i*integerColumns + j)*3 + 2];
         for (int k = 63; k >= 0; k--){
            long unsigned print0 = ((tmp0 >> k) & 1);
            long unsigned print1 = ((tmp1 >> k) & 1);
            long unsigned print2 = ((tmp2 >> k) & 1);
            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
            fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
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
         tmp0 = C[(i*integerColumns + j)*3 + 0];
         tmp1 = C[(i*integerColumns + j)*3 + 1];
         tmp2 = C[(i*integerColumns + j)*3 + 2];
         for (int k = 63; k >= 0; k--){
            long unsigned print0 = ((tmp0 >> k) & 1);
            long unsigned print1 = ((tmp1 >> k) & 1);
            long unsigned print2 = ((tmp2 >> k) & 1);
            long unsigned printInt = (print2 << 2) | (print1 << 1) | print0;
//            fullSizeMatrix[i*columns + j*64 + (63 - k)] = printInt;
            
            if (j*64 + (64 - k) >= columns){
               break;
            }
//            tmp2 = (tmp >> k*4) & 0xf;
            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
            if (!(printInt == Cf[i*columns + j * 64 + (63 - k)])){
               printf("row: %d\ncolumn: %d\n", i, j);
               printf("C: %lu\n", printInt);
               return 0;
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

int main(){
   srand(time(NULL));
   int m = 512;  // rows in A
   int lb = 512; // rows in B
   int la = (lb - 1) / 64 + 1;  // integer columns in A
   int nActual = 512;
   int n = (nActual - 1) / 64 + 1;  // integer columns in B
   //   int ra = 64 - (lb % 64); // zeros to pad on A
   //   int rb = (64 - (nActual % 64)) % 64; // zeros to pad on B
   int t = 2;  // four Russian divisor
   
   printf("m: %d\n", m);
   printf("lb: %d\n", lb);
   printf("la: %d\n", la);
   printf("n: %d\n", n);
   printf("t: %d\n", t);
   //   printf("ra: %d\n", ra);
   //   printf("rb: %d\n", rb);
   int lPadding = (64 - (lb % 64)) % 64;
   int nPadding = (64 - (nActual % 64)) % 64;
   printf("Lpadding: %d\n", lPadding);
   printf("Npadding: %d\n", nPadding);
   
   long unsigned int* A = malloc(m*la*3*(sizeof(long unsigned int)));
   long unsigned int* B = malloc((lb + lPadding)*n*3*(sizeof(long unsigned int)));
   long unsigned int* C = calloc(m*n*3, (sizeof(long unsigned int)));
   
   for (int i = 0; i < m; i++){
      for (int j = 0; j < la; j++){
         long unsigned int a2 = 0, a1 = 0, a0 = 0;
         if (j < la - 1){
            for (int i = 0; i < 64; i++){
               long unsigned int tmp = rand() % 7;
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
               long unsigned int tmp = rand() % 7;
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
   
   for (int i = 0; i < lb; i++){
      for (int j = 0; j < n; j++){
         long unsigned int b2 = 0, b1 = 0, b0 = 0;
         if (j < n - 1){
            for (int i = 0; i < 64; i++){
               long unsigned int tmp = rand() % 7;
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
               long unsigned int tmp = rand() % 7;
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
   
   int TRows = 1;
   for (int i = 0; i < t; i++){
      TRows *= 7;
   }
   int TSize = TRows*n;
   int numberOfTables = 64 / t;
   printf("TRows: %d\n# of tables: %d\n", TRows, numberOfTables);
   long unsigned int* T = malloc((TSize*(numberOfTables)+n)*3*sizeof(long int));
   //   createTable(T, B, TRows, n, t, numberOfTables, 0, 0, 0);
//   print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
   //print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns)
   
   long unsigned int bits0, bits1, bits2, index, index0, index1, index2;
   long unsigned int c0, c1, c2, t0, t1, t2, s0, s1, s2, s3, rMod5Fix;
   long unsigned int not7;
   int counter = 0;
   
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
      unsigned long long compensate = rdtsc();
      createTable(T, B, TRows, n, t, numberOfTables, i, 0, 0);
//      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
      startTime += (rdtsc() - compensate);
      //      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
      for (int j = 0; j < m; j++){
         bits0 = A[(j*la + i)*3 + 0];
         bits1 = A[(j*la + i)*3 + 1];
         bits2 = A[(j*la + i)*3 + 2];
//         printGF5(bits0, bits1, bits2);
         for (int k = numberOfTables - 1; k >= 0; k--){
            
            index0 = bits0 >> k*2;
            index1 = bits1 >> k*2;
            index2 = bits2 >> k*2;
//            printGF5(index0, index1, index2);
            index = ((index2 & 2) << 3) + ((index2 & 2) << 2) + ((index2 & 2) << 1)
            + ((index1 & 2) << 2) + ((index1 & 2) << 1) + ((index1 & 2))
            + ((index0 & 2) << 1) + ((index0 & 2)) + ((index0 & 2) >> 1)
            + ((index2 & 1) << 2) + ((index1 & 1) << 1) + (index0 & 1);
//            printf("bits\n");
//            printBinary(index2);
//            printBinary(index1);
//            printBinary(index0);
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
               s3 = (s2 & ((t1 & c1) | (s1 & (t0 & c0)))) ^ (t2 & c2);
               s2 ^= (s1 & (t0 & c0)) ^ t1 & c1;
               s1 ^= t0 & c0;
//               not7 = (s0 & s1 & s2);
               C[(j*n + l)*3 + 2] = (s2 ^ (s1 & (s0 & s3)));// & ~not7; //my meathod
               C[(j*n + l)*3 + 1] = (s1 ^ (s0 & s3));// & ~not7;// maybe try (s3 & ...) | (~s3 & ...)
               C[(j*n + l)*3 + 0] = (s0 ^ s3);// & ~not7;
//               printf("r: ");printGF5((~s2 & s0) | (s2 & s1 & ~s0) | (s3 & ~s2), (s3) | (s1 & s0) | (~s2 & s1), (s2 & ~s1 & ~s0));
            }
         }
         //                     printf("C\n");
         //                     printGF5Matrix(C, m, n);
      }
   }
   for (int j = 0; j < m; j++){
      for (int l = 0; l < n; l++){
         c0 = C[(j*n + l)*3 + 0];
         c1 = C[(j*n + l)*3 + 1];
         c2 = C[(j*n + l)*3 + 2];
         not7 = (c0 & c1 & c2);
         C[(j*n + l)*3 + 0] = c0 & ~not7;
         C[(j*n + l)*3 + 1] = c1 & ~not7;
         C[(j*n + l)*3 + 2] = c2 & ~not7;
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
         Cf[i * nActual + j] %= 7;
      }
   }
   endTime = rdtsc();
   //   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
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
   
   free(A);
   free(B);
   free(C);
   free(T);
}



