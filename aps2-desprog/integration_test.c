#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fourier.h"
#include "image.h"

#define NUM_TESTS 9

int main(void) {
  static char path[39];
  static double complex input[MAX_SIZE][MAX_SIZE], output[MAX_SIZE][MAX_SIZE];
  static image_t image;

  for (int t = 1; t <= NUM_TESTS; t++) {
    int x, y;

    sprintf(path, "integration_inputs/%d.pgm", t);
    load_image(path, &image);

    int width = pow(2, ceil(log2(image.width)));
    int height = pow(2, ceil(log2(image.height)));

    for (y = 0; y < image.height; y++) {
      for (x = 0; x < image.width; x++) {
        input[y][x] = image.levels[y][x];
      }
      for (; x < width; x++) {
        input[y][x] = 0;
      }
    }
    for (; y < height; y++) {
      for (x = 0; x < width; x++) {
        input[y][x] = 0;
      }
    }

    fft_forward_2d(input, width, height);

    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        output[y][x] = input[y][x];
      }
    }
    fft_inverse_2d(output, width, height);
    convert(output, &image);
    sprintf(path, "integration_outputs/%d.pgm", t);
    save_image(&image, path);

    filter_lp(input, output, width, height);
    fft_inverse_2d(output, width, height);
    convert(output, &image);
    sprintf(path, "integration_outputs/%d-lp.pgm", t);
    save_image(&image, path);

    filter_hp(input, output, width, height);
    fft_inverse_2d(output, width, height);
    convert(output, &image);
    sprintf(path, "integration_outputs/%d-hp.pgm", t);
    save_image(&image, path);
  }

  printf("tudo pronto para rodar integration_test.py\n");

  return EXIT_SUCCESS;
}
