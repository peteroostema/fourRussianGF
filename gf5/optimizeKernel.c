/* optimizeKernel.c is an extension of kernalt2.c,
 * but doens't run much quicker
 *
 * Again, this is for a bitpacked operation and slicing
 * would be ~2.5x faster. A bitsliced kernel was never made
 *
 * Author: Peter Oostema
 * Date: Summer 2017
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
   int m = 4;  // rows in A
   int lb = 16; // rows in B
   int la = (lb - 1) / 16 + 1;  // integer columns in A
   int nActual = 64;
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
   printf("A\n");
   printGF5Matrix(A, m, la);
   
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
   printf("B\n");
   printGF5Matrix(B, lb, n);
   
   int TRows = 1;
   for (int i = 0; i < t; i++){
      TRows *= 5;
   }
   int TSize = TRows*n;
   int numberOfTables = 16 / t;
   printf("TRows: %d\n# of tables: %d\n", TRows, numberOfTables);
   long unsigned int* T = malloc((TSize*(numberOfTables)+n)*sizeof(long int));
   long unsigned int* subDiv = T + 700; //+ (numberOfTables - 1)*TSize;
   
   unsigned long long startTime;//rdtsc()
   createTable(T, B, TRows, n, t, numberOfTables, 0, 0, 0);
   print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
   //print4Dgf5Matrix(long unsigned int* matrix, int rows, int majorColumns, int depth, int innerColumns, int outerColumns)
   long unsigned int* cIntermediateTest = calloc(m*n, (sizeof(long unsigned int)));
   
   long unsigned int a0, a1, a2, a3;
   long unsigned int index00, index01, index10, index11, index20, index21, index30, index31, index0, index1;
   long unsigned int index40, index41, index50, index51, index60, index61, index70, index71;
   long unsigned int index02, index03, index12, index13, index22, index23, index32, index33, index2, index3;
   long unsigned int index42, index43, index52, index53, index62, index63, index72, index73;
   __m256i xmm0, xmm1, ymm0, ymm1, xmm00, xmm01, xmm02, xmm03, xmm10, xmm11, xmm12, xmm13;
   __m256i xmm2, xmm3, ymm2, ymm3, xmm20, xmm21, xmm22, xmm23, xmm30, xmm31, xmm32, xmm33;
   __m256i mask0, mask1, mask2, mask3;
   mask0 = _mm256_set_epi64x (0x1111111111111111, 0x1111111111111111, 0x1111111111111111, 0x1111111111111111);
   mask1 = _mm256_set_epi64x (0x2222222222222222, 0x2222222222222222, 0x2222222222222222, 0x2222222222222222);
   mask2 = _mm256_set_epi64x (0x4444444444444444, 0x4444444444444444, 0x4444444444444444, 0x4444444444444444);
   mask3 = _mm256_set_epi64x (0x8888888888888888, 0x8888888888888888, 0x8888888888888888, 0x8888888888888888);
   
   startTime = rdtsc();
   
   
//   for (int i = 0; i < la; i++){
//      createTable(T, B, TRows, n, t, numberOfTables, i, 0, 0);
////      print4Dgf5Matrix(T, TRows, 1, numberOfTables, n, 0);
//      //      createTable(T, B, TRows, n, t, numberOfTables, i); // expand by n
//      //      print3DBinaryMatrix(T, TRows, n, numberOfTables);
//      for (int j = 0; j < m; j++){
   a0 = *(A);
   a1 = *(A + 1);
   a2 = *(A + 2);
   a3 = *(A + 3);
   
   
   index00 = ((a0 & 0xf) + ((a0 >> 4) & 0xf)*5) << 2;
   index01 = ((a1 & 0xf) + ((a1 >> 4) & 0xf)*5) << 2;
   index02 = ((a2 & 0xf) + ((a2 >> 4) & 0xf)*5) << 2;
   index03 = ((a3 & 0xf) + ((a3 >> 4) & 0xf)*5) << 2;
   VLOAD((subDiv + index00), xmm0);
   VLOAD((subDiv + index01), xmm1);
   VLOAD((subDiv + index02), xmm2);
   VLOAD((subDiv + index03), xmm3);

   index10 = a0 >> 8;
   index11 = a1 >> 8;
   index20 = a0 >> 16;
   index21 = a1 >> 16;
   index30 = a0 >> 24;
   index31 = a1 >> 24;
   index40 = a0 >> 32;
   index41 = a1 >> 32;
   index50 = a0 >> 40;
   index51 = a1 >> 40;
   index60 = a0 >> 48;
   index61 = a1 >> 48;
   index70 = a0 >> 56;
   index71 = a1 >> 56;
   
   index12 = a2 >> 8;
   index13 = a3 >> 8;
   index22 = a2 >> 16;
   index23 = a3 >> 16;
   index32 = a2 >> 24;
   index33 = a3 >> 24;
   index42 = a2 >> 32;
   index43 = a3 >> 32;
   index52 = a2 >> 40;
   index53 = a3 >> 40;
   index62 = a2 >> 48;
   index63 = a3 >> 48;
   index72 = a2 >> 56;
   index73 = a3 >> 56;
   
   index10 = ((index10 & 0xf) + ((index10 >> 4) & 0xf)*5) << 2;
   index11 = ((index11 & 0xf) + ((index11 >> 4) & 0xf)*5) << 2;
   index20 = ((index20 & 0xf) + ((index20 >> 4) & 0xf)*5) << 2;
   index21 = ((index21 & 0xf) + ((index21 >> 4) & 0xf)*5) << 2;
   index30 = ((index30 & 0xf) + ((index30 >> 4) & 0xf)*5) << 2;
   index31 = ((index31 & 0xf) + ((index31 >> 4) & 0xf)*5) << 2;
   index40 = ((index40 & 0xf) + ((index40 >> 4) & 0xf)*5) << 2;
   index41 = ((index41 & 0xf) + ((index41 >> 4) & 0xf)*5) << 2;
   index50 = ((index50 & 0xf) + ((index50 >> 4) & 0xf)*5) << 2;
   index51 = ((index51 & 0xf) + ((index51 >> 4) & 0xf)*5) << 2;
   index60 = ((index60 & 0xf) + ((index60 >> 4) & 0xf)*5) << 2;
   index61 = ((index61 & 0xf) + ((index61 >> 4) & 0xf)*5) << 2;
   index70 = ((index70 & 0xf) + ((index70 >> 4) & 0xf)*5) << 2;
   index71 = ((index71 & 0xf) + ((index71 >> 4) & 0xf)*5) << 2;
   
   index12 = ((index12 & 0xf) + ((index12 >> 4) & 0xf)*5) << 2;
   index13 = ((index13 & 0xf) + ((index13 >> 4) & 0xf)*5) << 2;
   index22 = ((index22 & 0xf) + ((index22 >> 4) & 0xf)*5) << 2;
   index23 = ((index23 & 0xf) + ((index23 >> 4) & 0xf)*5) << 2;
   index32 = ((index32 & 0xf) + ((index32 >> 4) & 0xf)*5) << 2;
   index33 = ((index33 & 0xf) + ((index33 >> 4) & 0xf)*5) << 2;
   index42 = ((index42 & 0xf) + ((index42 >> 4) & 0xf)*5) << 2;
   index43 = ((index43 & 0xf) + ((index43 >> 4) & 0xf)*5) << 2;
   index52 = ((index52 & 0xf) + ((index52 >> 4) & 0xf)*5) << 2;
   index53 = ((index53 & 0xf) + ((index53 >> 4) & 0xf)*5) << 2;
   index62 = ((index62 & 0xf) + ((index62 >> 4) & 0xf)*5) << 2;
   index63 = ((index63 & 0xf) + ((index63 >> 4) & 0xf)*5) << 2;
   index72 = ((index72 & 0xf) + ((index72 >> 4) & 0xf)*5) << 2;
   index73 = ((index73 & 0xf) + ((index73 >> 4) & 0xf)*5) << 2;
   VLOAD((T + 600 + index10), ymm0);
   VLOAD((T + 600 + index11), ymm1);
   VLOAD((T + 600 + index12), ymm2);
   VLOAD((T + 600 + index13), ymm3);// T[index*n + TSize*(numberOfTables - k - 1) + l];
//      printf("id0: %lu, full %lu\n", index10/4, index10);
//      printf("id1: %lu, full %lu\n", index11/4, index11);
//      printf("id2: %lu, full %lu\n", index12/4, index12);
//      printf("id3: %lu, full %lu\n", index13/4, index13);
//      _mm256_storeu_si256(cIntermediateTest, ymm0);
//      _mm256_storeu_si256(cIntermediateTest + n, ymm1);
//      _mm256_storeu_si256(cIntermediateTest + 2*n, ymm2);
//      _mm256_storeu_si256(cIntermediateTest + 3*n, ymm3);
//      printf("SecondLoad\n");
//      printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);//s = C[j*n + l] + T[index*n + TSize*(numberOfTables - k - 1) + l];
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);//s0 = s & 0x1111111111111111;
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);//s1 = (s & 0x2222222222222222) >> 1;
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);//2 = (s & 0x4444444444444444) >> 2;
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);//((s & 0x8888888888888888) >> 3)
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));//rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);//rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);//C[j*n + l] = (s - rMod5Fix);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   VLOAD((T + 500 + index20), ymm0);
   VLOAD((T + 500 + index21), ymm1);
   VLOAD((T + 500 + index22), ymm2);
   VLOAD((T + 500 + index23), ymm3);
   //   printf("id0: %lu, full %lu\n", index20/4, index20);
   //   printf("id1: %lu, full %lu\n", index21/4, index21);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("ThirdLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   VLOAD((T + 400 + index30), ymm0);
   VLOAD((T + 400 + index31), ymm1);
   VLOAD((T + 400 + index32), ymm2);
   VLOAD((T + 400 + index33), ymm3);
   //   printf("id0: %lu, full %lu\n", index0/4, index0);
   //   printf("id1: %lu, full %lu\n", index1/4, index1);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("ForthLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   
   VLOAD((T + 300 + index40), ymm0);
   VLOAD((T + 300 + index41), ymm1);
   VLOAD((T + 300 + index42), ymm2);
   VLOAD((T + 300 + index43), ymm3);
   //   printf("id0: %lu, full %lu\n", index00/4, index00);
   //   printf("id1: %lu, full %lu\n", index01/4, index01);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("FifthLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   VLOAD((T + 200 + index50), ymm0);
   VLOAD((T + 200 + index51), ymm1);
   VLOAD((T + 200 + index52), ymm2);
   VLOAD((T + 200 + index53), ymm3);
   //   printf("id0: %lu, full %lu\n", index10/4, index10);
   //   printf("id1: %lu, full %lu\n", index11/4, index11);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("SixthLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   VLOAD((T + 100 + index60), ymm0);
   VLOAD((T + 100 + index61), ymm1);
   VLOAD((T + 100 + index62), ymm2);
   VLOAD((T + 100 + index63), ymm3);
   //   printf("id0: %lu, full %lu\n", index20/4, index20);
   //   printf("id1: %lu, full %lu\n", index21/4, index21);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("SeventhLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
   
   
   VLOAD((T + 00 + index70), ymm0);
   VLOAD((T + 00 + index71), ymm1);
   VLOAD((T + 00 + index72), ymm2);
   VLOAD((T + 00 + index73), ymm3);
   //   printf("id0: %lu, full %lu\n", index30/4, index30);
   //   printf("id1: %lu, full %lu\n", index31/4, index31);
   //   _mm256_storeu_si256(cIntermediateTest, ymm0);
   //   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
   //   printf("LastLoad\n");
   //   printGF5Matrix(cIntermediateTest, m, n);
   xmm0 = _mm256_add_epi16 (xmm0, ymm0);
   xmm1 = _mm256_add_epi16 (xmm1, ymm1);
   xmm2 = _mm256_add_epi16 (xmm2, ymm2);
   xmm3 = _mm256_add_epi16 (xmm3, ymm3);
   xmm00 = _mm256_and_si256 (xmm0, mask0);
   xmm10 = _mm256_and_si256 (xmm1, mask0);
   xmm20 = _mm256_and_si256 (xmm2, mask0);
   xmm30 = _mm256_and_si256 (xmm3, mask0);
   xmm01 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask1), 1);
   xmm11 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask1), 1);
   xmm21 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask1), 1);
   xmm31 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask1), 1);
   xmm02 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask2), 2);
   xmm12 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask2), 2);
   xmm22 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask2), 2);
   xmm32 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask2), 2);
   xmm03 = _mm256_srli_epi16 (_mm256_and_si256 (xmm0, mask3), 3);
   xmm13 = _mm256_srli_epi16 (_mm256_and_si256 (xmm1, mask3), 3);
   xmm23 = _mm256_srli_epi16 (_mm256_and_si256 (xmm2, mask3), 3);
   xmm33 = _mm256_srli_epi16 (_mm256_and_si256 (xmm3, mask3), 3);
   ymm0 = _mm256_or_si256 (xmm03, _mm256_and_si256 (xmm02, _mm256_or_si256 (xmm00, xmm01)));
   ymm1 = _mm256_or_si256 (xmm13, _mm256_and_si256 (xmm12, _mm256_or_si256 (xmm10, xmm11)));
   ymm2 = _mm256_or_si256 (xmm23, _mm256_and_si256 (xmm22, _mm256_or_si256 (xmm20, xmm21)));
   ymm3 = _mm256_or_si256 (xmm33, _mm256_and_si256 (xmm32, _mm256_or_si256 (xmm30, xmm31)));
   ymm0 = _mm256_or_si256 (_mm256_slli_epi16 (ymm0, 2), ymm0);
   ymm1 = _mm256_or_si256 (_mm256_slli_epi16 (ymm1, 2), ymm1);
   ymm2 = _mm256_or_si256 (_mm256_slli_epi16 (ymm2, 2), ymm2);
   ymm3 = _mm256_or_si256 (_mm256_slli_epi16 (ymm3, 2), ymm3);
   xmm0 = _mm256_sub_epi16 (xmm0, ymm0);
   xmm1 = _mm256_sub_epi16 (xmm1, ymm1);
   xmm2 = _mm256_sub_epi16 (xmm2, ymm2);
   xmm3 = _mm256_sub_epi16 (xmm3, ymm3);
//         index2 = ((a2 & 0xf) + ((a2 >> 4) & 0xf)*5) << 3;
//         index3 = ((a3 & 0xf) + ((a3 >> 4) & 0xf)*5) << 3;
//   printf("a0: "); printGF5(a0); printf("\n");
//   printf("a1: "); printGF5(a1); printf("\n");
   
   
//         VLOAD((subDiv + index2), xmm2);
//         VLOAD((subDiv + index3), xmm3);
//   subDiv -= TSize;
   
//               printf("id0: %lu, full %lu\n", index00/4, index00);
//               printf("id1: %lu, full %lu\n", index01/4, index01);
//               _mm256_storeu_si256(cIntermediateTest, xmm0);
//               _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//               printf("FirstLoad\n");
//               printGF5Matrix(cIntermediateTest, m, n);
   
//   index0 = a0 >> 8;
//   index1 = a1 >> 8;
//         index2 = a2 >> 8;
//         index3 = a3 >> 8;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
//         index2 = ((index2 & 0xf) + ((index2 >> 4) & 0xf)*5) << 3;
//         index3 = ((index3 & 0xf) + ((index3 >> 4) & 0xf)*5) << 3;
   
   

   
//         VLOAD((subDiv + index2), ymm2);
//         VLOAD((subDiv + index3), ymm3);
  
   
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("Added\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
   
   
//   _mm256_storeu_si256(cIntermediateTest, xmm03);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm13);
//   printf("ForthBit\n");
//   printGF5Matrix(cIntermediateTest, m, n);
//   
//   _mm256_storeu_si256(cIntermediateTest, xmm02);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm12);
//   printf("ThirdBit\n");
//   printGF5Matrix(cIntermediateTest, m, n);
//   
//   _mm256_storeu_si256(cIntermediateTest, xmm01);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm11);
//   printf("SecondBit\n");
//   printGF5Matrix(cIntermediateTest, m, n);
//   
//   _mm256_storeu_si256(cIntermediateTest, xmm00);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm10);
//   printf("FirstBit\n");
//   printGF5Matrix(cIntermediateTest, m, n);

   
   
   
//   _mm256_storeu_si256(cIntermediateTest, ymm0);
//   _mm256_storeu_si256(cIntermediateTest + n, ymm1);
//   printf("Mod5Fix\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
  
//   subDiv -= TSize;
   
   
   
   
//   index0 = a0 >> 16;
//   index1 = a1 >> 16;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
//   subDiv -= TSize;
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
   
//   index0 = a0 >> 24;
//   index1 = a1 >> 24;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
//   subDiv -= TSize;
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
   
   
//   index0 = a0 >> 32;
//   index1 = a1 >> 32;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
//   subDiv -= TSize;
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
//   index0 = a0 >> 40;
//   index1 = a1 >> 40;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
//   subDiv -= TSize;
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
//   index0 = a0 >> 48;
//   index1 = a1 >> 48;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
//   subDiv -= TSize;
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);
   
//   index0 = a0 >> 56;
//   index1 = a1 >> 56;
//   index0 = ((index0 & 0xf) + ((index0 >> 4) & 0xf)*5) << 2;
//   index1 = ((index1 & 0xf) + ((index1 >> 4) & 0xf)*5) << 2;
   
   
//   uint8_t *val = (uint8_t*) &xmm0;
//   printf("xmm0: ");
//   for (int k = 0; k < 32; k++){
//      printf("%i", val[k] >> 4);
//      printf("%i", val[k] & 4);
//   }
//   printf("\n");
   
//   printf("id0: %lu, full %lu\n", index0/4, index0);
//   printf("id1: %lu, full %lu\n", index1/4, index1);
//   _mm256_storeu_si256(cIntermediateTest, xmm0);
//   _mm256_storeu_si256(cIntermediateTest + n, xmm1);
//   printf("CIntermediate\n");
//   printGF5Matrix(cIntermediateTest, m, n);

//   _mm256_and_si256 (ymm00, __m256i b)
//   _mm256_or_si256 (__m256i a, __m256i b)
//   _mm256_sub_epi16 (__m256i a, __m256i b)
//   s = C[j*n + l] + T[index*n + TSize*(numberOfTables - k - 1) + l];
//   s0 = s & 0x1111111111111111;
//   s1 = (s & 0x2222222222222222) >> 1;
//   s2 = (s & 0x4444444444444444) >> 2;
//   rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2) // (s2) & (s0 | s1)
//   rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//   C[j*n + l] = (s - rMod5Fix);
   
   _mm256_storeu_si256(C, xmm0);
   _mm256_storeu_si256(C + n, xmm1);
   _mm256_storeu_si256(C + 2*n, xmm2);
   _mm256_storeu_si256(C + 3*n, xmm3);
   
   
          
   
//         VXOR(xmm2, ymm2, xmm2);
//         VXOR(xmm3, ymm3, xmm3);
   
//         sibDiv -+ TSize;
//         
//         for (int k = numberOfTables - 1; k >= 0; k--){
//            
//            index = bits >> k*8;
//            index =  (index & 0xf) + ((index >> 4) & 0xf)*5;
//            //            printf("id: %lu\t%d\n", index, counter++);
//            for (int l = 0; l < n; l++){
//               s = C[j*n + l] + T[index*n + TSize*(numberOfTables - k - 1) + l];
//               s0 = s & 0x1111111111111111;
//               s1 = (s & 0x2222222222222222) >> 1;
//               s2 = (s & 0x4444444444444444) >> 2;
//               rMod5Fix = ((s & 0x8888888888888888) >> 3) | (s2 & s1) | (s2 & s0);// | (s0 & s1 & s2)
//               rMod5Fix = (rMod5Fix << 2) | (rMod5Fix);
//               C[j*n + l] = (s - rMod5Fix);
//            }
//         }
//         //                     printf("C\n");
//         //                     printGF5Matrix(C, m, n);
//      }
//   }
   unsigned long long endTime = rdtsc();
      printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   printf("cycles: %llu\n", endTime - startTime);
   
   printf("C\n");
   printGF5Matrix(C, m, n);
   
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
   
   printf("Cf\n");
   int rowCounter = 0;
   for (int i = 0; i < m; i++){
      for (int j = 0; j < nActual; j++){
         printf("%d", Cf[i*nActual + j]);
      }
      printf("\t");
      printf(" %d", rowCounter);
      printf("\n");
      rowCounter++;
   }
   printf("\n");
   
   printf("Correct = 1, Wrong = 0: %d\n", checkCorrectResult(C, Cf, m, nActual, n));
   
   free(A);
   free(B);
   free(C);
   free(T);
}



