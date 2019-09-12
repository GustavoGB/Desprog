#ifndef FOURIER_H
#define FOURIER_H

#include <complex.h>

#define MAX_SIZE 512

#define PI 3.141592

#define SIGMA 16

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n);
void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n);

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n);
void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n);

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height);
void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height);

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height);
void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height);

#endif
