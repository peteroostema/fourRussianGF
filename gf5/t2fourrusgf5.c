/* t2fourrusgf5.c is the same implementation as fourrusgf5.c
 * but with t = 2
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
   int m = 2048;  // rows in A
   int lb = 2048; // rows in B
   int la = (lb - 1) / 16 + 1;  // integer columns in A
   int nActual = 2048;
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
   
//   createTable(T, B, TRows, n, t, numberOfTables, 0, 0, 0);
//   print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
   //print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns)
   
   long unsigned int bits, index;
   long unsigned int s, s0, s1, s2, rMod5Fix;
   int counter = 0;
   
   unsigned long long startTime = rdtsc();
   for (int i = 0; i < la; i++){
//      unsigned long long compensate = rdtsc();
      createTable(T, B, TRows, n, t, numberOfTables, i, 0, 0);
//      startTime += (rdtsc() - compensate);
//      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
//      startTime = rdtsc();
      //      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
      for (int j = 0; j < m; j++){
         bits = A[j*la + i];
         for (int k = numberOfTables - 1; k >= 0; k--){
            
            index = bits >> k*8;
//            printf("idHigher: %lu\tidLower: %lu\t%d\n", (index & 0xf0) >> 4, index & 0xf, counter++);
//            index =  (index & 0xf) + ((index >> 4) & 0xf)*5; //(((index >> 4) & 0xf) << 2) ^ ((index >> 4) & 0xf); // *5
            index =  (index & 0xf) + (((index >> 4) & 0xf) << 2) + ((index >> 4) & 0xf);
//                        printf("id: %lu\t%d\n", index, counter++);
            for (int l = 0; l < n; l++){
               s = C[j*n + l] + T[index*n + TSize*(numberOfTables - k - 1) + l];
               s0 = s & 0x1111111111111111;
               s1 = (s & 0x2222222222222222) >> 1;
               s2 = (s & 0x4444444444444444) >> 2;
               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & (s1 | s0));//(s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
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
   
//   startTime = rdtsc();
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            Cf[i * nActual + j] += (Af[i * lb + k] * Bf[j + k * nActual]);
         }
         Cf[i * nActual + j] %= 5;
      }
   }
//   endTime = rdtsc();
//   //   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
//   printf("cycles: %llu\n", endTime - startTime);
   
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



