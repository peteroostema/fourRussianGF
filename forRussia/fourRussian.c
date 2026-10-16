/* fourRussina.c is my first try at implementing the
 * method of four Russians on GF2. It is the slowest
 * and most general result I have since it accepts any
 * value for t, but doesn't consider any HPC optimizations
 *
 * Although this program is one of my slowest most of the 
 * concepts outside of optimization are best shown here,
 * so I've included full documentation throughout this file.
 *
 * Author: Peter Oostema
 * Date: Summer 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* converts a GF2 matrix stored in row major order in 32 bit integers
 *  other versions contain a convertBinaryColumnMatrixToFullSize to convert column major matrixes
 * @param: fullSizeMatrix, matrix where each elements stored in a seperate unsignedint
 * @param: binaryMatrix, matrix with packed elements
 * @param: rows, # of rows in the matrixes
 * @param: columns, # of columns in the matrixes
 * @param: binaryColumns, # of 32 bit unsignedints to store a row of binaryMatrix
 * result: fullSizeMatrix contains every element from binaryMatrix stored in a seperate unsignedint
 */
void convertBinaryMatrixToFullSize(unsigned int* fullSizeMatrix, unsigned int* binaryMatrix, int rows, int columns, int binaryColumns){
   int tmp = 0;
   int rowCounter = 0;
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < binaryColumns; j++){
         tmp = binaryMatrix[i*binaryColumns + j];
         for (int k = 0; k < 32; k++){ // reads through a 32 bit integer
            if (j*32 + k >= columns){
               break;
            }
            if (tmp & 0x80000000){ // check last bit
               fullSizeMatrix[i*columns + j * 32 + k] = 1;
            } else {
               fullSizeMatrix[i*columns + j * 32 + k] = 0;
            }
            tmp <<= 1;
         }
      }
      //rowCounter++;
   }
   printf("size: %d\n", rows*columns);
   for (int i = 0; i < rows; i++){
      for (int j = 0; j < columns; j++){
         printf("%d", fullSizeMatrix[i*columns + j]);
      }
      printf("\t");
      printf(" %d", rowCounter);
      printf("\n");
      rowCounter++;
   }
   printf("\n");
}

/* checks that matrix C matches the result from a matrix multiplication by definition
 * @param: C, resulting matrix from four Russian matrix multiply
 * @param: Cf, matrix of a elementwise matrix multiplication
 * @param: rows, # of rows in the matrixes
 * @param: columns, # of columns in the matrixes
 * @param: binaryColumns, # of 32 bit unsignedints to store a row of binaryMatrix
 * result: returns 1 if it worked, 0 if it didn't
 */
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

/* print out a matrix of GF2 packed elements
 *  other versions have a print4DBinaryMatrix for printing out the T matrix since it has a few more dimensions to it
 */
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

/* power computes the result from a base to an exponent
 *  in other versions this isn't used as t is resricted
 *  to a power of 2 and this computation can be acconplished
 *  by one left shift
 */
int power(int base, int exp){
   int result = 1;
   for (int i = 0; i < exp; i++){
      result *= base;
   }
   return result;
}

/* finds the next row that should be added to
 * complete a table stored in greycode order
 *
 * This was never used since the Table can be
 * created in less than twice the time storing
 * in binary order, and has less overhead
 * when the table needs to be read from.
 */
int findNextGreyRow (int position, int depth){
   printf("position: %d\t", position);
   if (position % 2){
      return depth;
   } else {
      return findNextGreyRow(position/2, depth + 1);
   }
   return -1; //error
}

/* fill lookup table T with all row combinations of B
 * @param: T, lookup table
 * @param: Bt, B
 * @param: TRow, # of rows needed in T to store all combinations of t rows of B
 * @param: startRow, row to start reading from B
 *    (other version change this to numberOfTables and TableNumber as construction got more complicated)
 * @param: n, matrix dimension n
 * @param: t, four Russian divisor
 * result: T's 2^t rows are filled with all the combinations of t rows of B
 */
void createTable(unsigned int* T, unsigned int* Bt, int TRows, int startRow, int n, int t){
   for (int j = 0; j < n; j++){
      T[j] = 0; // initialize the first row to zero
   }
   for (int i = 1; i < TRows; i++){ // find what row to add to the next row of T
      int rowsToXOR = t - 1; //        takes O(lg(i)) time
      unsigned int position = i;
      while (rowsToXOR >= 0){
         //printf("position: %d\t", position);
         if (position % 2){
            break;
         }
         position /= 2;
         rowsToXOR--;
      }
      //printf("rowtoXor: %d\n", rowsToXOR);
      for (int j = 0; j < n; j++){
         for (int k = t - 1; k >= rowsToXOR; k--){
            if (k == t - 1){ // if its the first addition read add to the previous row in T
               T[i*n + j] = T[(i - 1)*n + j] ^ Bt[(startRow + k)*n + j];
            } else{  // otherwise add to the same row
               T[i*n + j] ^= Bt[(startRow + k)*n + j];
            }
            
         }
      }
   }
}

/* reads a # of bits from anywher in A
 * (in other versions it is limited to a multiple of 2 bits not crossing an integer
 *  and in most version an entire integer of A is just read in without this function)
 * @param: A, matrix to read from
 * @param: start, element # one wants to start with
 * @param: length, # of elements to read
 * result: return length number of elements starting at start
 */
unsigned int readBits(unsigned int* A, int start, int length){
   unsigned int matrixStartID = start / 32; // integer id to start from
   printf("matrixID: %d\n", matrixStartID);
   unsigned int startIntPos = start % 32; // position in integer to read from
   printf("intID: %d\n", startIntPos);
   printf("length: %d\n", length);
   if (startIntPos + length > 32){ // if you need to read across integer do a bunch of extar stuff
      unsigned int leftSide = A[matrixStartID] << startIntPos;
      printf("A[0]: %d\n", A[0]);
      printf("leftSide: %d\n", leftSide);
      leftSide >>= (32 - length);
      printf("leftSide: %d\n", leftSide);
      unsigned int rightSide = A[matrixStartID + 1] >> (64 - startIntPos - length);
      printf("rightSide: %d\n", rightSide);
      return leftSide | rightSide;
   } else{ // otherwise just read it out
      unsigned result = A[matrixStartID] << startIntPos;
      return result >> (32 - length);
   }
}

/* main gets its name from m being two letters before p,
 * so main is a little less of a pain, just as four Russian
 * matrix multiplication is a little less of a pain than normal
 * matrix multiplication \s
 */
int main(){
   int m = 6;  // rows in A, rows in C
   int lb = 6; // rows in B, columns in A
   int la = lb / 32 + 1;  // # of integer columns in A
   int nActual = 6; // columns in B and C
   int n = (nActual - 1) / 32 + 1;  // integer columns in B and C
   int ra = 32 - (lb % 32); // zeros to pad on A
   int rb = (32 - (nActual % 32)) % 32; // zeros to pad on B
   printf("m: %d\n", m);
   printf("lb: %d\n", lb);
   printf("la: %d\n", la);
   printf("n: %d\n", n);
   printf("ra: %d\n", ra);
   printf("rb: %d\n", rb);
   int t = 4;  // four Russian divisor
   unsigned int* A = malloc(m*la*(sizeof(int)));
   
   int lilExtra = (t - ((lb*n) % t)) % t; // since A is padded to 32x columns, B needs a lilExtra rows
   unsigned int* B = malloc((lb + lilExtra)*n*(sizeof(int))); // pad rows for createTable
   printf("SizeB: %d\n", (lb + lilExtra)*n);
   unsigned int* C = calloc(m*n,(sizeof(int)));
   srand(time(NULL)); // initialize random number generator
   
//   unsigned int* greyTable = malloc(power(2, t)*sizeof(int));
//   unsigned int currentPower = 2;
   //   greyTable[0] = 0;
   //   greyTable[1] = 1;
   //   for (int i = 1; i < t; i++){
   //      for (int j = 0; j < currentPower; j++){
   //         greyTable[currentPower + j] = greyTable[currentPower - j - 1];
   //         greyTable[currentPower + j] |= currentPower;
   //      }
   //      currentPower *= 2;
   //   }
   //   printf("greyTable\n");
   //   printBinaryMatrix(greyTable, power(2, t), 1);
   
   for (int i = 0; i < m; i++){// init A
      for (int j = 0; j < la; j++){
         if (j < la - 1){
            A[i*la + j] = rand();
         } else {
            A[i*la + j] = rand() >> ra;
            A[i*la + j] <<= ra;
         }
      }
   }
   printf("A\n");
   printBinaryMatrix(A, m, la);
   
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
   printf("B\n");
   printBinaryMatrix(B, lb + lilExtra, n);
   
   // read through A, t columns at a time so lb/t loops
   // other version read though 32 rows at a time so 32/t there
   int tLoops = lb / t;
   if (lb % t){
      tLoops++;
   }
   int TRows = power(2, t);
   printf("2^t: %d\n", TRows);
   unsigned int* T = malloc(TRows*n*sizeof(int));
   // time
   for (int i = 0; i < tLoops; i++){
      // create a lookup table for t columns/ t rows of B
      createTable(T, B, TRows, i*t, n, t);
      printBinaryMatrix(T, TRows, n);
      for (int j = 0; j < m; j++){
         // read t elements from A at once
         unsigned int id = readBits(A, i*t + 32*j*la, t);
         //id = (TRows - id) % TRows;
         printf("id: %d\n", id);
         for (int k = 0; k < n; k++){
            // lookup where the combination of rows of B are at the index read from A
            C[j*n + k] ^= T[id*n + k];
         }
      }
      printf("C\n");
      printBinaryMatrix(C, m, n);
      
   }
   // time
   printf("C\n");
   printBinaryMatrix(C, m, n);
   
   ///////////////////////////////////////////////////
   // testing                                        /
   //   matrix multiplication by definition          /
   //    take column of A, vector multiply with B    /
   //    add appropriatly                            /
   ///////////////////////////////////////////////////
   unsigned int* Af = malloc(m*lb*sizeof(int));
   unsigned int* Bf = malloc(lb*nActual*sizeof(int));
   unsigned int* Cf = calloc(m*nActual,sizeof(int));
   printf("Af\n");
   convertBinaryMatrixToFullSize(Af, A, m, lb, la);
   printf("Bf\n");
   convertBinaryMatrixToFullSize(Bf, B, lb, nActual, n);
   
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < nActual; j++){ // column
         for (int k = 0; k < lb; k++){ // interate
            Cf[i * nActual + j] ^= (Af[i * lb + k] & Bf[j + k * nActual]);
         }
      }
   }
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
//   free(greyTable);
}
