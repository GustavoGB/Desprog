#include <math.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  for (int k = 0; k < n; k++) {
    t[k] = 0;

    for (int j = 0; j < n; j++) {
      t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
    }
  }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  nft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  double complex sp[MAX_SIZE / 2];
  double complex si[MAX_SIZE / 2];
  double complex tp[MAX_SIZE / 2];
  double complex ti[MAX_SIZE / 2];
  int k = 0; // iterador para sp
  int h = 0; // iterador para si

  if (n == 1) {
    t[0] = s[0];
    return; //se n for 1
  }

  int metade = n / 2;

  for (int n1 = 0; n1 < metade; n1++) { // percorre para dividir os indices em pares e impares

    sp[k] = s[n1 * 2]; // pares
    k++;

    si[h] = s[n1 * 2 + 1]; // impares
    h++;

  } // depois disso já temos as listas

  for (int j = 0; j < metade; j++) {
    ti[j] += si[j] * cexp(sign * 2 * PI * j * I / n);
  }

  for (int o = 0; o < metade; o++) {
    tp[o] += sp[o] * cexp(sign * 2 * PI * o * I / n);
  }

  fft(si, ti, metade, sign);
  fft(sp, tp, metade, sign); // depois disso teremos tp e ti preenchidas

  for (int g = 0; g < metade; g++) {
    t[g] = tp[g] + ti[g] * cexp(sign * 2 * PI * g * I / n);
    t[g + (n / 2)] = tp[g] - ti[g] * cexp(sign * 2 * PI * g * I / n);
  }
}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  fft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  //row é height e col é width
  /*
para cada linha l da matriz
    aplique a transformada unidimensional normal sobre l
para cada coluna c da matriz
    aplique a transformada unidimensional normal sobre c
  */
  //para gravar os valores da transposta
  double complex tLinha[width];
  double complex sLinha[width];
  double complex tColuna[height];
  double complex sColuna[height];

  for (int l = 0; l < height; l++) { //percorre-se as linhas
    for (int c = 0; c < width; c++)  //e as colunas
    {
      sLinha[c] = matrix[l][c]; // Guarda em um "temp" para depois armazenar na matrix original
    }

    fft_forward(sLinha, tLinha, width); // Aplicamos a fft para as linhas
    for (int c = 0; c < width; c++) {
      matrix[l][c] = tLinha[c]; //guarda estes valores com a fft em nossa variavel temporaria tLinha
    }
  }

  // Aplicar a mesma logica porem invertendo linhas e colunas
  for (int c = 0; c < width; c++) {
    for (int l = 0; l < height; l++) {
      sColuna[l] = matrix[l][c];
    }
    fft_forward(sColuna, tColuna, height);
    for (int l = 0; l < height; l++) {
      matrix[l][c] = tColuna[l];
    }
  }
  return;
}
// Mesma logica utilizando fft inversa
void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  double complex tLinha[width];
  double complex sLinha[width];
  double complex tColuna[height];
  double complex sColuna[height];

  for (int l = 0; l < height; l++) { //percorre-se as linhas
    for (int c = 0; c < width; c++)  //e as colunas
    {
      sLinha[c] = matrix[l][c]; // Guarda em um "temp" para depois armazenar na matrix original
    }
    fft_inverse(sLinha, tLinha, width); // Aplicamos a fft para as linhas
    for (int c = 0; c < width; c++) {
      matrix[l][c] = tLinha[c]; //guarda estes valores com a fft em nossa variavel temporaria tLinha
    }
  }
  // Aplicar a mesma logica porem invertendo linhas e colunas
  for (int c = 0; c < width; c++) {
    for (int l = 0; l < height; l++) {
      sColuna[l] = matrix[l][c];
    }

    fft_inverse(sColuna, tColuna, height);
    for (int l = 0; l < height; l++) {
      matrix[l][c] = tColuna[l];
    }
  }
  return;
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
  int center_x = width / 2;
  int center_y = height / 2;

  double variance = -2 * SIGMA * SIGMA;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int dx = center_x - (x + center_x) % width;
      int dy = center_y - (y + center_y) % height;

      double d = dx * dx + dy * dy;

      double g = exp(d / variance);

      if (flip) {
        g = 1 - g;
      }

      output[y][x] = g * input[y][x];
    }
  }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 1);
}