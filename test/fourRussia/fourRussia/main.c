/* RCEREF.C == Referecne imoplementation of RC5-32/12/16 in C.
 * Copyright (C) 1995 RSA Data Security Inc.
 */

//SIEC keys 751 Bytes, 188 words

#include <stdio.h>
typedef unsigned int WORD; // should be 32-bit = 4 bytes
#define w 32                    // word size in bits
#define r 12                    // number of rounds
#define b 256                    // number of bytes in key
#define c 64                     // number of words in key, c = max(1,ceil(8*b/w))
#define t 26                    // size of table S = 2*(r+1) words
WORD S[t];                      // expended key table
WORD P = 0xb7e15163, Q = 0x9e3779b9; // magic constatns

// Rotation operators, x must be unsigned, to get logical right shift
#define ROTL(x,y) (((x)<<(y&(w-1))) | ((x)>>(w-(y&(w-1)))))
#define ROTR(x,y) (((x)>>(y&(w-1))) | ((x)<<(w-(y&(w-1)))))
void RC5_SETUP(unsigned char *K)
{
   // w = 32, r = 12, b = 16
   // c = max(1, ceil(8 * b/w))
   // t = 2 * (r+1)
   WORD i, j, k, u = w/8, A, B, L[c];
   
   for(i = b-1, L[c-1] = 0; i != -1; i--)
      L[i/u] = (L[i/u] << 8) + K[i];
   
   for(S[0] = P, i = 1; i < t; i++)
      S[i] = S[i-1] + Q;
   
   for(A = B = i = j = k = 0; k < 3 * t; k++, i = (i+1) % t, j = (j+1) % c) {
      A = S[i] = ROTL(S[i] + (A + B), 3);
      B = L[j] = ROTL(L[j] + (A + B), (A + B));
   }
}

void RC5_ENCRYPT(WORD *pt, WORD *ct) {
   WORD i, A = pt[0] + S[0], B = pt[1] + S[1];
   
   for(i = 1; i <= r; i++) {
      A = ROTL(A ^ B, B) + S[2*i];
      B = ROTL(B ^ A, A) + S[2*i + 1];
   }
   ct[0] = A; ct[1] = B;
}

void RC5_DECRYPT(WORD *ct, WORD *pt) {
   WORD i, B=ct[1], A=ct[0];
   
   for(i = r; i > 0; i--) {
      B = ROTR(B - S[2*i + 1], A) ^ A;
      A = ROTR(A - S[2*i], B) ^ B;
   }
   
   pt[1] = B - S[1]; pt[0] = A - S[0];
}

void printword(WORD A) {
   WORD k;
   for (k = 0; k < w; k+=8)
      printf("%x", (A >> k) & 0xFF);
}

int main(){
   WORD i, j, k, pt1[2], pt2[2], ct[2] = {0, 0};
   unsigned char key[b];
   if (sizeof(WORD) !=4)
      printf("RC5 error: WORD has %d byes.\n", sizeof(WORD));
   printf("RC5 examples:\n");
   for (i = 1; i < 6; i++){
      //Initialize pt  and key pseudorandomly based on previous ct
      pt1[0] = ct[0]; pt1[1] = ct[2];
      for (j = 0; j < b; j++) key[j] = ct[0]%(768-j);
      // Setup, encrypt, and decrypt
      RC5_SETUP(key);
      RC5_ENCRYPT(pt1, ct);
      RC5_DECRYPT(ct, pt2);
      // Print out results, checking for decryption failure
      printf("\n%d. key = ", i);
      for (j = 0; j < b; j++) printf("%x", key[j]);
      printf("\n     plaintext "); printword(pt1[0]); printword(pt1[1]);
      printf(" ---> ciphertext "); printword(ct[0]); printword(ct[1]);
      printf("\ndecipheredtext "); printword(pt2[0]); printword(pt2[1]);
      printf("\n");
      if (pt1[0]!=pt2[0] || pt1[1]!=pt2[1])
         printf("Dcryption Error");
   }
}



