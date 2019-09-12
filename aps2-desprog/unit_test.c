#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fourier.h"

#define NUM_TESTS 90

int main(void) {
  char path[33];
  double complex s[MAX_SIZE], nt[MAX_SIZE], ns[MAX_SIZE], ft[MAX_SIZE], fs[MAX_SIZE];

  srand(0);

  int num_passed = 0;

  for (int n = 1; n <= MAX_SIZE; n *= 2) {
    for (int t = 1; t <= NUM_TESTS; t++) {
      int j;

      for (j = 0; j < n; j++) {
        double level = (rand() / (double)RAND_MAX) * 256;

        if (level >= 256) {
          s[j] = 255;
        } else {
          s[j] = floor(level);
        }
      }

      nft_forward(s, nt, n);
      nft_inverse(nt, ns, n);

      fft_forward(s, ft, n);
      fft_inverse(ft, fs, n);

      int passed = 1;

      for (j = 0; j < n; j++) {
        if ((int)round(creal(ns[j])) != (int)round(creal(fs[j]))) {
          passed = 0;
          break;
        }
      }

      sprintf(path, "unit_outputs/%03d-%02d.txt", n, t);

      FILE *file = fopen(path, "w");

      fprintf(file, "SINAL     NFT                         FFT                         NFT INVERSA     FFT INVERSA\n");
      fprintf(file, "-----     -----------------------     -----------------------     -----------     -----------\n");

      for (j = 0; j < n; j++) {
        fprintf(file, " %4.0f     %11.3f %11.3f     %11.3f %11.3f            %4.0f            %4.0f\n", round(creal(s[j])), creal(nt[j]), cimag(nt[j]), creal(ft[j]), cimag(ft[j]), round(creal(ns[j])), round(creal(fs[j])));
      }

      fclose(file);

      if (passed) {
        num_passed++;
      } else {
        printf("NÃO PASSOU NO TESTE %d-%d: examine o texto %s\n", n, t, path);
      }
    }
  }

  printf("passou em %d de %d testes\n", num_passed, ((int)log2(MAX_SIZE) + 1) * NUM_TESTS);

  return EXIT_SUCCESS;
}
