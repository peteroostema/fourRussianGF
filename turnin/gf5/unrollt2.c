/* unrollt 1, 2, 4.c are unrolled bitpacked four Russian matrix multiplies.
 * Since the compiler really likes t*fourrusgf5.c and uses
 * instrinsics and unrolls automatically, unrollt*.c runs slower.
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static inline void createTable(long unsigned int* T, long unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){
      *(T + i*TRows) = 0;
   }
   int index = 0;
   
   for (int i = 0; i < numOfTables; i++){
      long unsigned int b0, b1, s, s0, s1, s2, s3, s4, s5, s10, s11, s12, s20, s21, s22, s30, s31, s32;
      long unsigned int s00, s01, s02, s40, s41, s42;
      long unsigned int d, d0, d1, d2;
      long unsigned int r, r0, r1, r2;
      r = 0;
      d = 0;
      long unsigned int rMod5Fix1, rMod5Fix2, rMod5Fix3, rMod5Fix4, rMod5Fix0, rMod5Fix;
//      printf("subtable: %d,   index: %d\n", i, index);
      
      for (int o = 0; o < n; o++){
         index = o + i*TRows*n;
         b0 = Bt[tableNumber*t*n*numOfTables + i*t*n + 1*n + o];
         b1 = Bt[tableNumber*t*n*numOfTables + i*t*n + 0 + o];
         //         for (int j = 0; j < 5; j++){
         //            for (int k = 0; k < 5; k++){
         
         
         s = b0;
//               index += n;
         s0 = s & 0x1111111111111111;
         s1 = (s & 0x2222222222222222) >> 1;
         s2 = (s & 0x4444444444444444) >> 2;
//         rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
//         rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//         T[index + n] = (s - rMod5Fix);
         T[index + n] = b0;
         
         r |= (~s1 & s0) << 2; //my meathod
         r |= s1 << 1;
         r |= (s2 | s1) & ~s0;
         
         T[index + 4*n] = r;
         
         d |= (s1 & ~s0) << 2; //my meathod
         d |= ((~s1 & s0) | s2) << 1;
         d |= s2 | (s1 & s0);
         T[index + 2*n] = d;
         d0 = d & 0x1111111111111111;
         d1 = (d & 0x2222222222222222) >> 1;
         d2 = (d & 0x4444444444444444) >> 2;
         
         r = (~d1 & d0) << 2; //my meathod
         r |= d1 << 1;
         r |= (d2 | d1) & ~d0;
         
         T[index + 3*n] = r;
         
         
         
//         T[index + 5*n] = T[index + 0*n];
//         T[index + 6*n] = T[index + 1*n];
//         T[index + 7*n] = T[index + 2*n];
//         T[index + 8*n] = T[index + 3*n];
//         T[index + 9*n] = T[index + 4*n];
//         T[index + 10*n] = T[index + 0*n];
//         T[index + 11*n] = T[index + 1*n];
//         T[index + 12*n] = T[index + 2*n];
//         T[index + 13*n] = T[index + 3*n];
//         T[index + 14*n] = T[index + 4*n];
//         T[index + 15*n] = T[index + 0*n];
//         T[index + 16*n] = T[index + 1*n];
//         T[index + 17*n] = T[index + 2*n];
//         T[index + 18*n] = T[index + 3*n];
//         T[index + 19*n] = T[index + 4*n];
//         T[index + 10*n] = T[index + 0*n];
//         T[index + 21*n] = T[index + 1*n];
//         T[index + 22*n] = T[index + 2*n];
//         T[index + 23*n] = T[index + 3*n];
//         T[index + 24*n] = T[index + 4*n];
         
         for (int j = 0; j < 4; j++){
            s0 = T[index] + b1;
            s1 = T[index + 1*n] + b1;
            s2 = T[index + 2*n] + b1;
            s3 = T[index + 3*n] + b1;
            s4 = T[index + 4*n] + b1;
            index += 5*n;
            s00 = s0 & 0x1111111111111111;
            s01 = (s0 & 0x2222222222222222) >> 1;
            s02 = (s0 & 0x4444444444444444) >> 2;
            s10 = s1 & 0x1111111111111111;
            s11 = (s1 & 0x2222222222222222) >> 1;
            s12 = (s1 & 0x4444444444444444) >> 2;
            s20 = s2 & 0x1111111111111111;
            s21 = (s2 & 0x2222222222222222) >> 1;
            s22 = (s2 & 0x4444444444444444) >> 2;
            s30 = s3 & 0x1111111111111111;
            s31 = (s3 & 0x2222222222222222) >> 1;
            s32 = (s3 & 0x4444444444444444) >> 2;
            s40 = s4 & 0x1111111111111111;
            s41 = (s4 & 0x2222222222222222) >> 1;
            s42 = (s4 & 0x4444444444444444) >> 2;
            rMod5Fix0 = ((s0 & 0x8888888888888888) >> 3) | (s02 & s01) | (s02 & s00);
            rMod5Fix1 = ((s1 & 0x8888888888888888) >> 3) | (s12 & s11) | (s12 & s10);
            rMod5Fix2 = ((s2 & 0x8888888888888888) >> 3) | (s22 & s21) | (s22 & s20);
            rMod5Fix3 = ((s3 & 0x8888888888888888) >> 3) | (s32 & s31) | (s32 & s30);
            rMod5Fix4 = ((s4 & 0x8888888888888888) >> 3) | (s42 & s41) | (s42 & s40);
            rMod5Fix0 = (rMod5Fix0 << 2) | (rMod5Fix0);
            rMod5Fix1 = (rMod5Fix1 << 2) | (rMod5Fix1);
            rMod5Fix2 = (rMod5Fix2 << 2) | (rMod5Fix2);
            rMod5Fix3 = (rMod5Fix3 << 2) | (rMod5Fix3);
            rMod5Fix4 = (rMod5Fix4 << 2) | (rMod5Fix4);
            T[index + 0*n] = (s0 - rMod5Fix0);
            T[index + 1*n] = (s1 - rMod5Fix1);
            T[index + 2*n] = (s2 - rMod5Fix2);
            T[index + 3*n] = (s3 - rMod5Fix3);
            T[index + 4*n] = (s4 - rMod5Fix4);
         }
               
               
               
               
//            s = T[index] + b1;
//            s0 = s & 0x1111111111111111;
//            s1 = (s & 0x2222222222222222) >> 1;
//            s2 = (s & 0x4444444444444444) >> 2;
//            rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
//            rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//            T[index] = (s - rMod5Fix);
         }
         //               s = T[index] + b2;
         //               s0 = s & 0x1111111111111111;
         //               s1 = (s & 0x2222222222222222) >> 1;
         //               s2 = (s & 0x4444444444444444) >> 2;
         //               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
         //               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
         //               T[index] = (s - rMod5Fix);
         //            }
         //            s = T[index] + b3;
         //            s0 = s & 0x1111111111111111;
         //            s1 = (s & 0x2222222222222222) >> 1;
         //            s2 = (s & 0x4444444444444444) >> 2;
         //            rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
         //            rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
         //            T[index] = (s - rMod5Fix);
         //         }
         //   //      index++;
      }
//   printf("\nfinal index: %d\n", index);
}

static inline void createTableOld(long unsigned int* T, long unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256){
   for (int i = 0; i < numOfTables; i++){
      *(T + i*TRows) = 0;
   }
   int index = 0;
   
   for (int i = 0; i < numOfTables; i++){
      long unsigned int b0, b1, b2, b3, s0, s1, s2, s3, s, rMod5Fix;
      //      printf("subtable: %d,   index: %d\n", i, index);
      
      for (int o = 0; o < n; o++){
         index = o + i*TRows*n;
         b0 = Bt[tableNumber*t*n*numOfTables + i*t*n + 1*n + o];
         b1 = Bt[tableNumber*t*n*numOfTables + i*t*n + 0 + o];
         //         for (int j = 0; j < 5; j++){
         //            for (int k = 0; k < 5; k++){
         for (int l = 0; l < 5; l++){
            for (int m = 0; m < 5; m++){
               s = T[index] + b0;
               index += n;
               s0 = s & 0x1111111111111111;
               s1 = (s & 0x2222222222222222) >> 1;
               s2 = (s & 0x4444444444444444) >> 2;
               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
               T[index] = (s - rMod5Fix);
            }
            s = T[index] + b1;
            s0 = s & 0x1111111111111111;
            s1 = (s & 0x2222222222222222) >> 1;
            s2 = (s & 0x4444444444444444) >> 2;
            rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
            rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
            T[index] = (s - rMod5Fix);
         }
         //               s = T[index] + b2;
         //               s0 = s & 0x1111111111111111;
         //               s1 = (s & 0x2222222222222222) >> 1;
         //               s2 = (s & 0x4444444444444444) >> 2;
         //               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
         //               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
         //               T[index] = (s - rMod5Fix);
         //            }
         //            s = T[index] + b3;
         //            s0 = s & 0x1111111111111111;
         //            s1 = (s & 0x2222222222222222) >> 1;
         //            s2 = (s & 0x4444444444444444) >> 2;
         //            rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
         //            rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
         //            T[index] = (s - rMod5Fix);
         //         }
         //   //      index++;
      }
   }
   //   printf("\nfinal index: %d\n", index);
}

void print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns){
   long unsigned int tmp = 0;
   int rowCounter = 0;
   int subtableSize = rows*(majorColumns*innerColumns + outerColumns);
   for (int l = 0; l < depth; l++){
      printf("subTable: %d\nsubtableSize: %d\n", l, subtableSize);
      
      for (int i = 0; i < rows; i++){
         for (int j = 0; j < innerColumns; j++){
            tmp = matrix[l*subtableSize + j + i*innerColumns];
            for (int k = 15; k >= 0; k--){
               long unsigned tmp2 = (tmp >> k*4) & 0xf;
               //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
               printf("%lu", tmp2); //"   %lu"
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
   long unsigned int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         tmp = matrix[i*columns + j];
         for (int i = 15; i >= 0; i--){
            long unsigned tmp2 = (tmp >> i*4) & 0xf;
            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
            printf("%lu", tmp2); //"   %lu"
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

void printGF5(long unsigned int r){
   for (int i = 15; i >= 0; i--){
      long unsigned int tmp = (r >> i*4) & 0xf;
      //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
      printf("   %lu", tmp);
   }
   printf("\n");
}

void convertGF5MatrixToFullSize(unsigned int* fullSizeMatrix, long unsigned int* gf5Matrix, int rows, int columns, int integerColumns){
   long unsigned int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp = gf5Matrix[i*integerColumns + j];
         for (int k = 15; k >= 0; k--){
            fullSizeMatrix[i*columns + j*16 + (15 - k)] = (tmp >> k*4) & 0xf;
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
   long unsigned int tmp2;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < integerColumns; j++){
         tmp = C[i*integerColumns + j];
         for (int k = 15; k >= 0; k--){
            if (j*16 + (16 - k) >= columns){
               break;
            }
            tmp2 = (tmp >> k*4) & 0xf;
            //      long unsigned int printingNumber = tmp2*4 + tmp1*2 + tmp0;
            if (!(tmp2 == Cf[i*columns + j * 16 + (15 - k)])){
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
   int m = 1024;  // rows in A
   int lb = 1024; // rows in B
   int la = (lb - 1) / 16 + 1;  // integer columns in A
   int nActual = 1024;
   int n = (nActual - 1) / 16 + 1;  // integer columns in B
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
   int lPadding = (16 - (lb % 16)) % 16;
   int nPadding = (16 - (nActual % 16)) % 16;
   printf("Lpadding: %d\n", lPadding);
   printf("Npadding: %d\n", nPadding);
   
   long unsigned int* A = malloc(m*la*(sizeof(long unsigned int)));
   long unsigned int* B = malloc((lb + lPadding)*n*(sizeof(long unsigned int)));
   long unsigned int* C = calloc(m*n, (sizeof(long unsigned int)));
   
   for (int i = 0; i < m; i++){
      for (int j = 0; j < la; j++){
         long unsigned int a = 0;
         if (j < la - 1){
            for (int i = 0; i < 16; i++){
               long unsigned int tmp = rand() % 5;
               a <<= 4;
               a |= tmp;
            }
         } else{
            for (int i = 0; i < 16 - lPadding; i++){
               long unsigned int tmp = rand() % 5;
               a <<= 4;
               a |= tmp;
            }
            a <<= lPadding*4;
         }
         A[i*la + j] = a;
      }
   }
//   printf("A\n");
//   printGF5Matrix(A, m, la);
   
   for (int i = 0; i < lb; i++){
      for (int j = 0; j < n; j++){
         long unsigned int b = 0;
         if (j < n - 1){
            for (int i = 0; i < 16; i++){
               long unsigned int tmp = rand() % 5;
               b <<= 4;
               b |= tmp;
            }
         } else{
            for (int i = 0; i < 16 - nPadding; i++){
               long unsigned int tmp = rand() % 5;
               b <<= 4;
               b |= tmp;
            }
            b <<= nPadding*4;
         }
         B[i*n + j] = b;
      }
   }
   for (int i = 0; i < lPadding; i++){
      for (int j = 0; j < n; j++){
         B[(i+lb)*n + j] = 0;
      }
   }
//   printf("B\n");
//   printGF5Matrix(B, lb, n);
   
   int TRows = 1;
   for (int i = 0; i < t; i++){
      TRows *= 5;
   }
   int TSize = TRows*n;
   int numberOfTables = 16 / t;
   printf("TRows: %d\n# of tables: %d\n", TRows, numberOfTables);
   long unsigned int* T = malloc((TSize*(numberOfTables)+n)*sizeof(long int));
//   long unsigned int* TOld = malloc((TSize*(numberOfTables)+n)*sizeof(long int));
   unsigned long long startTime = rdtsc();
//   createTable(T, B, TRows, n, t, numberOfTables, 0, 0, 0);
//   print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
   //print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns)
   
   long unsigned int bits, index;
   long unsigned int s, s0, s1, s2, rMod5Fix;
   int counter = 0;
   
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i, 0, 0);
//      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
//      createTableOld(TOld, B, TRows, n, t, numberOfTables, i, 0, 0);
//      print4Dgf5Matrix(TOld, TRows, 1, numberOfTables, n, 0);
      //      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
      for (int j = 0; j < m; j++){
         for (int k = numberOfTables - 1; k >= 0; k--){
            bits = A[j*la + i];
            bits >>= k*8;
            index =  (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            //            printf("id: %lu\t%d\n", index, counter++);
            for (int l = 0; l < n; l++){
               s = C[j*n + l] + T[index*n + TSize*(numberOfTables - k - 1) + l];
               s0 = s & 0x1111111111111111;
               s1 = (s & 0x2222222222222222) >> 1;
               s2 = (s & 0x4444444444444444) >> 2;
               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
               C[j*n + l] = (s - rMod5Fix);
            }
            bits = A[j*la + i];
            index7 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index6 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index5 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index4 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index3 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index2 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index1 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            bits >>= 8;
            index0 = (bits & 0xf) + ((bits >> 4) & 0xf)*5;
            for (int l = 0; l < n; l++){
               s0 = C[j*n + l] + T[index0*n + TSize*0 + l];
               s1 = C[j*n + l] + T[index1*n + TSize*1 + l];
               s2 = C[j*n + l] + T[index2*n + TSize*2 + l];
               s3 = C[j*n + l] + T[index3*n + TSize*3 + l];
               s4 = C[j*n + l] + T[index4*n + TSize*4 + l];
               s5 = C[j*n + l] + T[index5*n + TSize*5 + l];
               s6 = C[j*n + l] + T[index6*n + TSize*6 + l];
               s7 = C[j*n + l] + T[index7*n + TSize*7 + l];
               s00 = s0 & 0x1111111111111111;
               s01 = (s0 & 0x2222222222222222) >> 1;
               s02 = (s0 & 0x4444444444444444) >> 2;
               s10 = s1 & 0x1111111111111111;
               s11 = (s1 & 0x2222222222222222) >> 1;
               s12 = (s1 & 0x4444444444444444) >> 2;
               s20 = s2 & 0x1111111111111111;
               s21 = (s2 & 0x2222222222222222) >> 1;
               s22 = (s2 & 0x4444444444444444) >> 2;
               s30 = s3 & 0x1111111111111111;
               s31 = (s3 & 0x2222222222222222) >> 1;
               s32 = (s3 & 0x4444444444444444) >> 2;
               s40 = s4 & 0x1111111111111111;
               s41 = (s4 & 0x2222222222222222) >> 1;
               s42 = (s4 & 0x4444444444444444) >> 2;
               s50 = s5 & 0x1111111111111111;
               s51 = (s5 & 0x2222222222222222) >> 1;
               s52 = (s5 & 0x4444444444444444) >> 2;
               s60 = s6 & 0x1111111111111111;
               s61 = (s6 & 0x2222222222222222) >> 1;
               s62 = (s6 & 0x4444444444444444) >> 2;
               s70 = s7 & 0x1111111111111111;
               s71 = (s7 & 0x2222222222222222) >> 1;
               s72 = (s7 & 0x4444444444444444) >> 2;
               rMod5Fix0 = ((s0 & 0x8888888888888888) >> 3) | (s02 & s01) | (s02 & s00);// | (s0 & s1 & s2)
               rMod5Fix1 = ((s1 & 0x8888888888888888) >> 3) | (s12 & s11) | (s12 & s10);
               rMod5Fix2 = ((s2 & 0x8888888888888888) >> 3) | (s22 & s21) | (s22 & s20);
               rMod5Fix3 = ((s3 & 0x8888888888888888) >> 3) | (s32 & s31) | (s32 & s30);
               rMod5Fix4 = ((s4 & 0x8888888888888888) >> 3) | (s42 & s41) | (s42 & s40);
               rMod5Fix5 = ((s5 & 0x8888888888888888) >> 3) | (s52 & s51) | (s52 & s50);
               rMod5Fix6 = ((s6 & 0x8888888888888888) >> 3) | (s62 & s61) | (s62 & s60);
               rMod5Fix7 = ((s7 & 0x8888888888888888) >> 3) | (s72 & s71) | (s72 & s70);
               rMod5Fix0 = (rMod5Fix0 << 2) | (rMod5Fix0);
               rMod5Fix1 = (rMod5Fix1 << 2) | (rMod5Fix1);
               rMod5Fix2 = (rMod5Fix2 << 2) | (rMod5Fix2);
               rMod5Fix3 = (rMod5Fix3 << 2) | (rMod5Fix3);
               rMod5Fix4 = (rMod5Fix4 << 2) | (rMod5Fix4);
               rMod5Fix5 = (rMod5Fix5 << 2) | (rMod5Fix5);
               rMod5Fix6 = (rMod5Fix6 << 2) | (rMod5Fix6);
               rMod5Fix7 = (rMod5Fix7 << 2) | (rMod5Fix7);
               C[j*n + l] = (s - rMod5Fix);
            }
         }
         //                     printf("C\n");
         //                     printGF5Matrix(C, m, n);
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



