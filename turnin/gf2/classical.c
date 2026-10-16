/* This file was my starting point and just initializes
 * three matrixies and multiplies them without modulus
 * then prints them out
 * Author: Peter Oostema
 * Date: Summer 2017
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
   int m = 7;  // rows in A
   int l = 2;  // columns in A and rows in B
   int n = 4;  // columns in B
   int t = 2;  // four Russian divisor
   int* A = malloc(m*l*(sizeof(int)));
   int* B = malloc(l*n*(sizeof(int)));
   int* C = calloc(m*n,(sizeof(int)));
   int b = 0; // counter for printing
   srand(time(NULL));
   for (int i = 0; i < m*l; i++){
      A[i] = i;
      //A[i] = rand();
      
      printf("%d", A[i]);
      if (b < l - 1){
         printf("\t");
         b++;
      } else{
         printf("\n");
         b = 0;
      }
   }
   printf("\n");
   for (int i = 0; i < l*n; i++){
      B[i] = i;
      
      printf("%d", B[i]);
      if (b < n - 1){
         printf("\t");
         b++;
      } else{
         printf("\n");
         b = 0;
      }
   }
   printf("\n");
   
//   for (int i = 0; i < m*n; i++){
//      C[i] = 0;
//   }
   //C[1] = 0;
   //C[3] = 0;
   
//   for (int i = 0; i < m*n; i++){
//      
//      printf("%d", C[i]);
//      if (b < n - 1){
//         printf("\t");
//         b++;
//      } else{
//         printf("\n");
//         b = 0;
//      }
//   }
//   printf("\n");
   
   //printf("idk\n");
   //printf("size of char:%lu \n", sizeof(char));
   //printf("size of long:%lu \n", sizeof(long));
   //printf("size of long long:%lu \n", sizeof(long long));
   //int a = 2000000000;
   //printf("two billion %d\n", a);
   
   for (int i = 0; i < m; i++){ // row
      for (int j = 0; j < n; j++){ // column
         for (int k = 0; k < l; k++){ // interate
            C[i * n + j] = A[i * l + k] * B[j + k * n] + C[i * n + j];
         }
         printf("%d", C[i * n + j]);
         
         if (b < n - 1){
            printf("\t");
            b++;
         } else{
            printf("\n");
            b = 0;
         }
      }
   }
   free(A);
   free(B);
   free(C);
}
