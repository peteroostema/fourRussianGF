/* assembly2.c is an older version of assebly.c
 * except the createTable function doesn't compile.
 * But there might be something of use here so I
 * haven't deleted it
 * Author: Peter Oostema
 * Date: Summer 2017
 */

/* *
#ifndef _ISA_INTRINSICS_AVX_4X32F_H
#define _ISA_INTRINSICS_AVX_4X32F_H

#include <immintrin.h>

#define INTTYPE       unsigned long
#define DATATYPE      double
#define VECTSZ        4
#define VECTORTYPE    __m256d


#define PREFETCH_NONE(x)

// Source --> destination
#define VLOAD_IA(srcaddr,dest)           \
asm volatile(                            \
"vmovapd %[addr], %[rvdest]\n"       \
: [rvdest] "=x"(dest)                \
: [addr] "m"(*srcaddr)               \
);

#define PREFETCH_T0_IA(srcaddr)             \
asm volatile(                           \
"prefetcht0 %[addr]\n"              \
:                                   \
: [addr] "m"(*srcaddr)              \
);


#define PREFETCH_T1_IA(srcaddr)             \
asm volatile(                           \
"prefetcht1 %[addr]\n"              \
:                                   \
: [addr] "m"(*srcaddr)              \
);

#define PREFETCH_T2_IA(srcaddr)             \
asm volatile(                           \
"prefetcht2 %[addr]\n"              \
:                                   \
: [addr] "m"(*srcaddr)              \
);

#define VBROADCAST_IA(srcaddr,dest)        \
asm volatile(                              \
"vbroadcastsd %[addr], %[rvdest]\n"    \
: [rvdest] "=x"(dest)                  \
: [addr] "m"(*srcaddr)                 \
);

#define VMUL_IA(srca,srcb,dest)             \
asm volatile(                               \
"vmulpd %[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
: [vsrca] "x"(srca),                    \
[vsrcb] "x"(srcb)                     \
);

#define VADD_IA(srca,srcb,dest)             \
asm volatile(                               \
"vaddpd %[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
: [vsrca] "x"(srca),                    \
[vsrcb] "x"(srcb)                     \
);


#define VZERO_IA(dest)                      \
asm volatile(                               \
"vxorpd %[vdest],%[vdest],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
);

#define VPERM2F128_IA(shf,srca,srcb,dest)                   \
asm volatile(                                           \
"vperm2f128 %[ishf],%[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                                \
: [ishf]  "i"(shf),                                 \
[vsrca] "x"(srca),                                \
[vsrcb] "x"(srcb)                                 \
);

#define VSHUFFLE_IA(shf,srca,srcb,dest)                     \
asm volatile(                                           \
"vshufpd %[ishf],%[vsrca],%[vsrcb],%[vdest]\n"      \
: [vdest] "=x"(dest)                                \
: [ishf]  "i"(shf),                                 \
[vsrca] "x"(srca),                                \
[vsrcb] "x"(srcb)                                 \
);


#define DECREMENT(dest)                                   \
asm volatile(                                         \
"sub $1, %[_dest]\n"                                     \
: [_dest] "+r"(dest)                                  \
);


#define VLOAD_IN(srcaddr,dest)      dest = _mm256_load_pd(srcaddr)
#define VBROADCAST_IN(srcaddr,dest) dest = _mm256_set1_pd(*(srcaddr))
#define VMUL_IN(srca,srcb,dest)     dest = _mm256_mul_pd(srca,srcb)
#define VADD_IN(srca,srcb,dest)     dest = _mm256_add_pd(srca,srcb)

#endif // _ISA_INTRINSICS_AVX_4X32F_H

vmovups (%rax,%rcx,4), %ymm0
vmovaps %ymm0, 0x240(%rsp)
movq    0x348(%rsp), %rax
movq    0x320(%rsp), %rcx
movl    0x18(%rcx), %edx
movl    %edx, %ecx
leaq    (%rax,%rcx,4), %rsi
movq    %rsi, 0x1788(%rsp)
vmovups (%rax,%rcx,4), %ymm0
vmovaps %ymm0, 0x220(%rsp)
movq    0x348(%rsp), %rax
movq    0x320(%rsp), %rcx
movl    0x1c(%rcx), %edx
movl    %edx, %ecx
leaq    (%rax,%rcx,4), %rsi
movq    %rsi, 0x1780(%rsp)
vmovups (%rax,%rcx,4), %ymm0
vmovaps %ymm0, 0x200(%rsp)
vmovaps 0x2e0(%rsp), %ymm0
movq    0x340(%rsp), %rax
movq    0x318(%rsp), %rcx
movl    (%rcx), %edx
movl    %edx, %ecx
leaq    (%rax,%rcx,4), %rsi
movq    %rsi, 0x1778(%rsp)
vmovups (%rax,%rcx,4), %ymm1
vmovaps %ymm0, 0x1740(%rsp)
vmovaps %ymm1, 0x1720(%rsp)
vmovaps 0x1740(%rsp), %ymm0
vpxor   %ymm1, %ymm0, %ymm0
vmovaps %ymm0, 0x2e0(%rsp)
vmovaps 0x2c0(%rsp), %ymm0
movq    0x340(%rsp), %rax
movq    0x318(%rsp), %rcx
movl    0x4(%rcx), %edx
movl    %edx, %ecx
leaq    (%rax,%rcx,4), %rsi
movq    %rsi, 0x1718(%rsp)
vmovups (%rax,%rcx,4), %ymm1
vmovaps %ymm0, 0x16e0(%rsp)
vmovaps %ymm1, 0x16c0(%rsp)
vmovaps 0x16e0(%rsp), %ymm0
vpxor   %ymm1, %ymm0, %ymm0
vmovaps %ymm0, 0x2c0(%rsp)
vmovaps 0x2a0(%rsp), %ymm0
movq    0x340(%rsp), %rax
movq    0x318(%rsp), %rcx
movl    0x8(%rcx), %edx
movl    %edx, %ecx
leaq    (%rax,%rcx,4), %rsi
movq    %rsi, 0x16b8(%rsp)
vmovups (%rax,%rcx,4), %ymm1
vmovaps %ymm0, 0x1680(%rsp)
vmovaps %ymm1, 0x1660(%rsp)
vmovaps 0x1680(%rsp), %ymm0

/*   */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <immintrin.h>

//static inline void FIRST_TRY(){
#define FIRST_TRY()
asm volatile(
             "movq    0x340(%rsp), %rax\n\t"
             "movq    0x318(%rsp), %rcx\n\t"
             "movl    (%rcx), %edx\n\t"
             "movl    %edx, %ecx\n\t"
             "leaq    (%rax,%rcx,4), %rsi\n\t"
             "movq    %rsi, 0x1778(%rsp)\n\t"
             "vmovups (%rax,%rcx,4), %ymm1\n\t"
             "vmovaps %ymm0, 0x1740(%rsp)\n\t"
             "vmovaps %ymm1, 0x1720(%rsp)\n\t"
             "vmovaps 0x1740(%rsp), %ymm0\n\t"
             "vpxor   %ymm1, %ymm0, %ymm0\n\t"
             "vmovaps %ymm0, 0x2e0(%rsp)\n\t"
             "vmovaps 0x2c0(%rsp), %ymm0\n\t"
             );

#define VXOR_IA(srca,srcb,dest)             \
asm volatile(                               \
"vpxor %[vsrca],%[vsrcb],%[vdest]\n"   \
: [vdest] "=x"(dest)                    \
: [vsrca] "x"(srca),                    \
[vsrcb] "x"(srcb)                     \
);

#define VLOAD_IA(srcaddr,dest)           \
asm volatile(                            \
"vmovdqu %[addr], %[rvdest]\n"       \
: [rvdest] "=x"(dest)                \
: [addr] "m"(*srcaddr)               \
);


static inline unsigned long long rdtsc(){
   unsigned long long int x;
   unsigned a, d;
   
   __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));
   
   return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

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

static inline void createTable(unsigned int* T, unsigned int* Bt, int TRows, int n, int t, int numOfTables, int tableNumber, int nq256, int nr256, int section){// get rid of n store column majorish
   int start = section*8;
   if (section < nq256){
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
         for (int k = t - 1; k >= rowsToXOR; k--){
            for (int l = 0; l < numOfTables; l++){
               if (k == t - 1){
                  _mm256_storeu_si256(T + l*TRows*8 + i*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*8 + (i - 1)*8), _mm256_load_si256 (Bt + (l*t + k)*n + start + tableNumber*32*n)));
               } else{
                  _mm256_storeu_si256(T + l*TRows*8 + i*8, _mm256_xor_si256 (_mm256_load_si256 (T + l*TRows*8 + i*8), _mm256_load_si256 (Bt + (l*t + k)*n + start + tableNumber*32*n)));
               }
            }
         }
      }
   } else {
      for (int i = 0; i < numOfTables; i++){
         for (int j = 0; j < nr256; j++){
            T[i*TRows*nr256 + j*TRows] = 0;
         }
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
         for (int k = t - 1; k >= rowsToXOR; k--){
            for (int j = 0; j < nr256; j++){
               for (int l = 0; l < numOfTables; l++){
                  if (k == t - 1){
                     T[l*TRows*8 + i + j*TRows] = T[l*TRows*8 + (i-1) + j*TRows] ^ Bt[(l*t + k)*n + j + tableNumber*32*n + nq256*8];
                  } else{
                     T[l*TRows*8 + i + j*TRows] ^= Bt[(l*t + k)*n + j + tableNumber*32*n + nq256*8];
                  }
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

void print4DBinaryMatrix(unsigned int* matrix, int rows, int majorColumns, int depth, innerColumns, int outerColumns){
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
   
   int m = 1024;  // rows in A
   int lb = 1024; // rows in B
   int la = (lb - 1) / 32 + 1;  // columns in A
   int nActual = 1024;
   int n = (nActual - 1) / 32 + 1;  // columns in B
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
//   printf("m: %d\n", m);
//   printf("lb: %d\n", lb);
//   printf("la: %d\n", la);
//   printf("n: %d\n", n);
//   printf("ra: %d\n", ra);
//   printf("rb: %d\n", rb);
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
   printf("256 bit reads: %d\n", readsPerColumn);
   printf("leftOver: %d\n", rowsLeftover);
   printf("B is %d * 256bits + %d ints\n", nq256, nr256);
   unsigned long long startTime = rdtsc();
   {
   //assume tables have been created
   
      uint64_t *aptr = (uint64_t*)A;
      uint64_t *bptr = (uint64_t*)T;
      
      uint64_t a0 = *aptr;
      uint64_t a1 = *(aptr + 1);
      uint64_t a2 = *(aptr + 2);
      uint64_t a3 = *(aptr + 3);
   
   uint64_t mask = 0xF;
   
   __m256i xmm0, xmm1, xmm2, xmm3, ymm0, ymm1, ymm2, ymm3;
   
   xmm0 = _mm256_setzero_si256();
   xmm1 = _mm256_setzero_si256();
   xmm2 = _mm256_setzero_si256();
   xmm3 = _mm256_setzero_si256();
   
   for (int j = 0; j < n/8; ++j){
                  for (int k = 0; k < m/4; ++k){
                     
          createTable(T, B, 16, n, 4, 16, i, nq256, nr256, j);
                     
         for (int i = 0; i < la/2; i++){
            //create 16 new tables of 2^4 entries for new 64 bit integer of A
            //assuming 4 bit chunks  (t=4) => 64/4 = 16
            for (int ii = 0; ii < 16; ++ii){
               uint64_t index0 = ((a0 >> (60-ii*4)) &  mask) << 2;
               VLOAD((bptr + index0), ymm0);
               VXOR(xmm0, ymm0, xmm0);
         
               uint64_t index1 = ((a1 >> (60-ii*4)) &  mask) << 2;
               VLOAD((bptr + index1), ymm1);
               VXOR(xmm1, ymm1, xmm1);
         
               uint64_t index2 = ((a2 >> (60-ii*4)) &  mask) << 2;
               VLOAD((bptr + index2), ymm2);
               VXOR(xmm2, ymm2, xmm2);
         
               uint64_t index3 = ((a3 >> (60-ii*4)) &  mask) << 2;
               VLOAD((bptr + index3), ymm3);
               VXOR(xmm3, ymm3, xmm3);
            }
      
            _mm256_storeu_si256(C + 0*n + j*8+k*4*n,
                                _mm256_xor_si256 (xmm0, _mm256_load_si256(C +0*n+ j*8+k*4*n)));
            _mm256_storeu_si256(C + 1*n + j*8+k*4*n,
                                _mm256_xor_si256 (xmm0, _mm256_load_si256(C +1*n+ j*8+k*4*n)));
            _mm256_storeu_si256(C + 2*n + j*8+k*4*n,
                                _mm256_xor_si256 (xmm0, _mm256_load_si256(C +2*n+ j*8+k*4*n)));
            _mm256_storeu_si256(C + 3*n + j*8+k*4*n,
                                _mm256_xor_si256 (xmm0, _mm256_load_si256(C +3*n+ j*8+k*4*n)));
         }
      }
   }
      
}

   
   for (int i = 0; i < la; i++){
      createTable(T, B, TRows, n, t, numberOfTables, i, nq256, nr256);
//      printBinaryMatrix(T, 8, 1);
//      printBinaryMatrix(T, TRows, n);
//      print4DBinaryMatrix(T, TRows, nq256, numberOfTables, 8, nr256);
      
      register unsigned int* subDiv1 = T;
      register unsigned int* subDiv2 = T + TSize;
      register unsigned int* subDiv3 = T + 2*TSize;
      register unsigned int* subDiv4 = T + 3*TSize;
      for (int j = 0; j < readsPerColumn; j++){
         id4 = _mm256_load_si256(A + j*8 + i*m);
         //         unsigned int* idInt = (unsigned int*) &id4;
         //         printf("Integers: \n");
         //         printBinary(idInt[0]); printBinary(idInt[1]); printBinary(idInt[2]); printBinary(idInt[3]); printBinary(idInt[4]); printBinary(idInt[5]); printBinary(idInt[6]); printBinary(idInt[7]);
         //         printf("IntegersAgain killme: %u %u %u %u %u %u %u %u\n", idInt[0], idInt[1], idInt[2], idInt[3], idInt[4], idInt[5], idInt[6], idInt[7]);
         //__m256i id = _mm256_load_ps(A + j/8 + i*m);
         id1 = _mm256_srli_epi32 (id4, 24);
         id1 = _mm256_slli_epi32(id1, 3);
         register unsigned int* id1Int = (unsigned int*) &(id1);
         id4 = _mm256_slli_epi32(id4, 8);
         id2 = _mm256_srli_epi32 (id4, 24);
         id2 = _mm256_slli_epi32(id2, 3);
         register unsigned int* id2Int = (unsigned int*) &(id2);
         id4 = _mm256_slli_epi32(id4, 8);
         id3 = _mm256_srli_epi32 (id4, 24);
         id3 = _mm256_slli_epi32(id3, 3);
         register unsigned int* id3Int = (unsigned int*) &(id3);
         id4 = _mm256_slli_epi32(id4, 8);
         id4 = _mm256_srli_epi32 (id4, 24);
         id4 = _mm256_slli_epi32(id4, 3);
         register unsigned int* id4Int = (unsigned int*) &(id4);
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
         register __m256i xmm0 = _mm256_setzero_si256();
         register __m256i xmm1 = _mm256_setzero_si256();
         register __m256i xmm2 = _mm256_setzero_si256();
         register __m256i xmm3 = _mm256_setzero_si256();
         register __m256i xmm4 = _mm256_setzero_si256();
         register __m256i xmm5 = _mm256_setzero_si256();
         register __m256i xmm6 = _mm256_setzero_si256();
         register __m256i xmm7 = _mm256_setzero_si256();
         

         
         for (int k = 0; k < nq256; k++){ //nActual / 256
            //            id = _mm256_slli_epi32 (id, int imm8);
            
            //prefetch C
//            xmm0 = _mm256_setzero_si256();
//            xmm1 = _mm256_setzero_si256();
//            xmm2 = _mm256_setzero_si256();
//            xmm3 = _mm256_setzero_si256();
//            xmm4 = _mm256_setzero_si256();
//            xmm5 = _mm256_setzero_si256();
//            xmm6 = _mm256_setzero_si256();
//            xmm7 = _mm256_setzero_si256();
         
//            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv1 + (*(id1Int + 0))*n));
//            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv1 + (*(id1Int + 1))*n));
//            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv1 + (*(id1Int + 2))*n));
//            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv1 + (*(id1Int + 3))*n));
//            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv1 + (*(id1Int + 4))*n));
//            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv1 + (*(id1Int + 5))*n));
//            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv1 + (*(id1Int + 6))*n));
//            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv1 + (*(id1Int + 7))*n));

            
//            xmm0 = _mm256_load_si256 (subDiv1 + (*(id1Int + 0)));
//            xmm1 = _mm256_load_si256 (subDiv1 + (*(id1Int + 1)));
//            xmm2 = _mm256_load_si256 (subDiv1 + (*(id1Int + 2)));
//            xmm3 = _mm256_load_si256 (subDiv1 + (*(id1Int + 3)));
//            xmm4 = _mm256_load_si256 (subDiv1 + (*(id1Int + 4)));
//            xmm5 = _mm256_load_si256 (subDiv1 + (*(id1Int + 5)));
//            xmm6 = _mm256_load_si256 (subDiv1 + (*(id1Int + 6)));
//            xmm7 = _mm256_load_si256 (subDiv1 + (*(id1Int + 7)));//id1Int++
//            id1Int -= 8;
            //            printf("C\n");
            //            printBinaryMatrix(C, m, n);
            
            __m256i a0, a1, a2, a3, a4, a5, a6, a7;
            VLOAD_IA((subDiv1 + (*(id1Int + 0))), xmm0);
            VLOAD_IA((subDiv1 + (*(id1Int + 1))), xmm1);
            VLOAD_IA((subDiv1 + (*(id1Int + 2))), xmm2);
            VLOAD_IA((subDiv1 + (*(id1Int + 3))), xmm3);
            VLOAD_IA((subDiv1 + (*(id1Int + 4))), xmm4);
            VLOAD_IA((subDiv1 + (*(id1Int + 5))), xmm5);
            VLOAD_IA((subDiv1 + (*(id1Int + 6))), xmm6);
            VLOAD_IA((subDiv1 + (*(id1Int + 7))), xmm7);
            
            VLOAD_IA((subDiv2 + (*(id2Int + 0))), a0);
            VXOR_IA(xmm0, a0, xmm0);
            VLOAD_IA((subDiv2 + (*(id2Int + 1))), a1);
            VXOR_IA(xmm1, a1, xmm1);
            VLOAD_IA((subDiv2 + (*(id2Int + 2))), a2);
            VXOR_IA(xmm2, a2, xmm2);
            VLOAD_IA((subDiv2 + (*(id2Int + 3))), a3);
            VXOR_IA(xmm3, a3, xmm3);
            VLOAD_IA((subDiv2 + (*(id2Int + 4))), a4);
            VXOR_IA(xmm4, a4, xmm4);
            VLOAD_IA((subDiv2 + (*(id2Int + 5))), a5);
            VXOR_IA(xmm5, a5, xmm5)
            VLOAD_IA((subDiv2 + (*(id2Int + 6))), a6);
            VXOR_IA(xmm6, a6, xmm6);
            VLOAD_IA((subDiv2 + (*(id2Int + 7))), a7);
            VXOR_IA(xmm7, a7, xmm7);
            
            VLOAD_IA((subDiv3 + (*(id3Int + 0))), a0);
            VXOR_IA(xmm0, a0, xmm0);
            VLOAD_IA((subDiv3 + (*(id3Int + 1))), a1);
            VXOR_IA(xmm1, a1, xmm1);
            VLOAD_IA((subDiv3 + (*(id3Int + 2))), a2);
            VXOR_IA(xmm2, a2, xmm2);
            VLOAD_IA((subDiv3 + (*(id3Int + 3))), a3);
            VXOR_IA(xmm3, a3, xmm3);
            VLOAD_IA((subDiv3 + (*(id3Int + 4))), a4);
            VXOR_IA(xmm4, a4, xmm4);
            VLOAD_IA((subDiv3 + (*(id3Int + 5))), a5);
            VXOR_IA(xmm5, a5, xmm5);
            VLOAD_IA((subDiv3 + (*(id3Int + 6))), a6);
            VXOR_IA(xmm6, a6, xmm6);
            VLOAD_IA((subDiv3 + (*(id3Int + 7))), a7);
            VXOR_IA(xmm7, a7, xmm7);
            
            VLOAD_IA((subDiv4 + (*(id4Int + 0))), a0);
            VXOR_IA(xmm0, a0, xmm0);
            VLOAD_IA((subDiv4 + (*(id4Int + 1))), a1);
            VXOR_IA(xmm1, a1, xmm1);
            VLOAD_IA((subDiv4 + (*(id4Int + 2))), a2);
            VXOR_IA(xmm2, a2, xmm2);
            VLOAD_IA((subDiv4 + (*(id4Int + 3))), a3);
            VXOR_IA(xmm3, a3, xmm3);
            VLOAD_IA((subDiv4 + (*(id4Int + 4))), a4);
            VXOR_IA(xmm4, a4, xmm4);
            VLOAD_IA((subDiv4 + (*(id4Int + 5))), a5);
            VXOR_IA(xmm5, a5, xmm5);
            VLOAD_IA((subDiv4 + (*(id4Int + 6))), a6);
            VXOR_IA(xmm6, a6, xmm6);
            VLOAD_IA((subDiv4 + (*(id4Int + 7))), a7);
            VXOR_IA(xmm7, a7, xmm7);

            
//            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv2 + (*(id2Int + 0))));
//            xmm1 = _mm256_xor_si256 (xmm1, _mm256_load_si256 (subDiv2 + (*(id2Int + 1))));
//            xmm2 = _mm256_xor_si256 (xmm2, _mm256_load_si256 (subDiv2 + (*(id2Int + 2))));
//            xmm3 = _mm256_xor_si256 (xmm3, _mm256_load_si256 (subDiv2 + (*(id2Int + 3))));
//            xmm4 = _mm256_xor_si256 (xmm4, _mm256_load_si256 (subDiv2 + (*(id2Int + 4))));
//            xmm5 = _mm256_xor_si256 (xmm5, _mm256_load_si256 (subDiv2 + (*(id2Int + 5))));
//            xmm6 = _mm256_xor_si256 (xmm6, _mm256_load_si256 (subDiv2 + (*(id2Int + 6))));
//            xmm7 = _mm256_xor_si256 (xmm7, _mm256_load_si256 (subDiv2 + (*(id2Int + 7))));
            
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
            //            C[8*n*j + n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[1] ));//wrong idInt
            //            C[8*n*j + 2*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[2] ));//wrong idInt
            //            C[8*n*j + 3*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[3] ));//wrong idInt
            //            C[8*n*j + 4*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[4] ));//wrong idInt
            //            C[8*n*j + 5*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[5] ));//wrong idInt
            //            C[8*n*j + 6*n + k] ^= (unsigned int*) (_mm256_load_si256(T + idInt[6] ));//wrong idInt
            //            C[8*n*j + 7*n + k] ^=  (_mm256_load_si256(T + idInt[7] ));//wrong idInt
            
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VXOR_IA(xmm0, a0, xmm0);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + n + 8*k), a1);
//            VXOR_IA(xmm1, a1, xmm1);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 2*n + 8*k), a2);
//            VXOR_IA(xmm2, a2, xmm2);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 3*n + 8*k), a3);
//            VXOR_IA(xmm3, a3, xmm3);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 4*n + 8*k), a4);
//            VXOR_IA(xmm4, a0, xmm4);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 5*n + 8*k), a5);
//            VXOR_IA(xmm5, a1, xmm5);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 6*n + 8*k), a6);
//            VXOR_IA(xmm6, a2, xmm6);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
//            VLOAD_IA((C + 8*n*j + 7*n + 8*k), a7);
//            VXOR_IA(xmm7, a3, xmm7);
//            VLOAD_IA((C + 8*n*j + 8*k), a0);
            
            _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + 8*n*j + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (xmm1, _mm256_load_si256 (C + 8*n*j + n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (xmm2, _mm256_load_si256 (C + 8*n*j + 2*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (xmm3, _mm256_load_si256 (C + 8*n*j + 3*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (xmm4, _mm256_load_si256 (C + 8*n*j + 4*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (xmm5, _mm256_load_si256 (C + 8*n*j + 5*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (xmm6, _mm256_load_si256 (C + 8*n*j + 6*n + 8*k)));
            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (xmm7, _mm256_load_si256 (C + 8*n*j + 7*n + 8*k)));
            
//            _mm256_storeu_si256(C + 8*n*j + 8*k, xmm0);
//            _mm256_storeu_si256(C + 8*n*j + n + 8*k, xmm1);
//            _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, xmm2);
//            _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, xmm3);
//            _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, xmm4);
//            _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, xmm5);
//            _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, xmm6);
//            _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, xmm7);
            
            subDiv1 += 8*TRows;
            subDiv2 += 8*TRows;
            subDiv3 += 8*TRows;
            subDiv4 += 8*TRows;
         }
//         subDiv1 = T;
//         subDiv2 = T + TSize;
//         subDiv3 = T + 2*TSize;
//         subDiv4 = T + 3*TSize;
         
         id1 = _mm256_srli_epi32(id1, 3);
         id1Int = (unsigned int*) &(id1);
         id2 = _mm256_srli_epi32(id2, 3);
         id2Int = (unsigned int*) &(id2);
         id3 = _mm256_srli_epi32(id3, 3);
         id3Int = (unsigned int*) &(id3);
         id4 = _mm256_srli_epi32(id4, 3);
         id4Int = (unsigned int*) &(id4);
         
//         _mm256_storeu_si256(C + 8*n*j + 8*k, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + 8*n*j + 8*k)));// sequential C read write
//         _mm256_storeu_si256(C + 8*n*j + n + 8*k, _mm256_xor_si256 (xmm1, _mm256_load_si256 (C + 8*n*j + n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 2*n + 8*k, _mm256_xor_si256 (xmm2, _mm256_load_si256 (C + 8*n*j + 2*n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 3*n + 8*k, _mm256_xor_si256 (xmm3, _mm256_load_si256 (C + 8*n*j + 3*n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 4*n + 8*k, _mm256_xor_si256 (xmm4, _mm256_load_si256 (C + 8*n*j + 4*n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 5*n + 8*k, _mm256_xor_si256 (xmm5, _mm256_load_si256 (C + 8*n*j + 5*n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 6*n + 8*k, _mm256_xor_si256 (xmm6, _mm256_load_si256 (C + 8*n*j + 6*n + 8*k)));
//         _mm256_storeu_si256(C + 8*n*j + 7*n + 8*k, _mm256_xor_si256 (xmm7, _mm256_load_si256 (C + 8*n*j + 7*n + 8*k)));

         for (int l = 0; l < nr256; l++){
            C[8*n*j + 8*nq256 + l] ^= *(subDiv1 + *(id1Int + 0)); //T[id1Int[0]*n + l + 8*nq256];
            C[8*n*j + 8*nq256 + n + l] ^= *(subDiv1 + *(id1Int + 1));
            C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv1 + *(id1Int + 2));
            C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv1 + *(id1Int + 3));
            C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv1 + *(id1Int + 4));
            C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv1 + *(id1Int + 5));
            C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv1 + *(id1Int + 6));
            C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv1 + *(id1Int + 7));
            C[8*n*j + 8*nq256 + l] ^= *(subDiv2 + *(id2Int + 0));
            C[8*n*j + 8*nq256 + n + l] ^= *(subDiv2 + *(id2Int + 1));
            C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv2 + *(id2Int + 2));
            C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv2 + *(id2Int + 3));
            C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv2 + *(id2Int + 4));
            C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv2 + *(id2Int + 5));
            C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv2 + *(id2Int + 6));
            C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv2 + *(id2Int + 7));
            C[8*n*j + 8*nq256 + l] ^= *(subDiv3 + *(id3Int + 0));
            C[8*n*j + 8*nq256 + n + l] ^= *(subDiv3 + *(id3Int + 1));
            C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv3 + *(id3Int + 2));
            C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv3 + *(id3Int + 3));
            C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv3 + *(id3Int + 4));
            C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv3 + *(id3Int + 5));
            C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv3 + *(id3Int + 6));
            C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv3 + *(id3Int + 7));
            C[8*n*j + 8*nq256 + l] ^= *(subDiv4 + *(id4Int + 0));
            C[8*n*j + 8*nq256 + n + l] ^= *(subDiv4 + *(id4Int + 1));
            C[8*n*j + 8*nq256 + 2*n + l] ^= *(subDiv4 + *(id4Int + 2));
            C[8*n*j + 8*nq256 + 3*n + l] ^= *(subDiv4 + *(id4Int + 3));
            C[8*n*j + 8*nq256 + 4*n + l] ^= *(subDiv4 + *(id4Int + 4));
            C[8*n*j + 8*nq256 + 5*n + l] ^= *(subDiv4 + *(id4Int + 5));
            C[8*n*j + 8*nq256 + 6*n + l] ^= *(subDiv4 + *(id4Int + 6));
            C[8*n*j + 8*nq256 + 7*n + l] ^= *(subDiv4 + *(id4Int + 7)); //T[id4Int[7]*n + l + 8*nq256 + 3*TSize]
            
//            T[id2*n + TSize + l];
//            C[(j+readsPerColumn)*8*n + l] ^= T[id3*n + TSize*2 + l];
//            C[(j+readsPerColumn)*8*n + l] ^= T[id4*n + TSize*3 + l];
            
            subDiv1 += TRows;
            subDiv2 += TRows;
            subDiv3 += TRows;
            subDiv4 += TRows;
         }
         subDiv1 = T;
         subDiv2 = T + TSize;
         subDiv3 = T + 2*TSize;
         subDiv4 = T + 3*TSize;
//         printf("C\n");
//         printBinaryMatrix(C, m, n);
      }
      for (int j = 0; j < rowsLeftover; j++){
         unsigned int id1 = readBits(A, j*32 + i*32*m + 256*readsPerColumn, t); //  + i*32*m//
         unsigned int id2 = readBits(A, t + j*32 + i*32*m + 256*readsPerColumn, t);
         unsigned int id3 = readBits(A, 2*t + j*32 + i*32*m + 256*readsPerColumn, t);
         unsigned int id4 = readBits(A, 3*t + j*32 + i*32*m + 256*readsPerColumn, t);
//         printf("id1: %d\n", id1);
//         printf("id2: %d\n", id2);
//         printf("id3: %d\n", id3);
//         printf("id4: %d\n", id4);
//                        printf("APostition: %d\n", k*t + j*32 + i*32*m);
//                        printf("id: %d, Table#: %d\n", id, k);
//                        printf("CPosition: %d\n", j*n + i);
//                        printf("TValue: ");
//                        printBinary(T[id + TSize*k]);
//                        printf("\n");
         register __m256i xmm0 = _mm256_setzero_si256();
         
         
         for (int k = 0; k < nq256; k++){
            xmm0 = _mm256_load_si256 (subDiv1 + k*TRows*8 + id1*8);
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv2 + k*TRows*8 + id2*8));
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv3 + k*TRows*8 + id3*8));
            xmm0 = _mm256_xor_si256 (xmm0, _mm256_load_si256 (subDiv4 + k*TRows*8 + id4*8));
            _mm256_storeu_si256(C + (j+8*readsPerColumn)*n + k*8, _mm256_xor_si256 (xmm0, _mm256_load_si256 (C + (j+8*readsPerColumn)*n + k*8)));
         }
         for (int k = 0; k < nr256; k++){
            C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv1 + k*TRows + nq256*8*TRows + id1);//T[id1 + l];
            C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv2 + k*TRows + nq256*8*TRows + id2);//T[id2 + TSize + l];
            C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv3 + k*TRows + nq256*8*TRows + id3);//T[id3 + TSize*2 + l];
            C[(j+8*readsPerColumn)*n + nq256*8 + k] ^= *(subDiv4 + k*TRows + nq256*8*TRows + id4);//T[id4 + TSize*3 + l];
         }
         
            
            
//         for (int l = 0; l < n; l++){
//            C[(j+8*readsPerColumn)*n + l] ^= *(subDiv1 + l*TRows + id1);//T[id1 + l];
//            C[(j+8*readsPerColumn)*n + l] ^= *(subDiv2 + l*TRows + id2);//T[id2 + TSize + l];
//            C[(j+8*readsPerColumn)*n + l] ^= *(subDiv3 + l*TRows + id3);//T[id3 + TSize*2 + l];
//            C[(j+8*readsPerColumn)*n + l] ^= *(subDiv4 + l*TRows + id4);//T[id4 + TSize*3 + l];
//         }
      }
      //      printf("\n");
      
   }
   unsigned long long endTime = rdtsc();
   printf("time: %f s\n", (float)(endTime - startTime) / 2800000000.0);
   
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
//   printf("A\n");
//   printBinaryColumnMatrix(A, m, la);
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