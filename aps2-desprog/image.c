#include <math.h>
#include <stdio.h>

#include "image.h"

void load_image(char path[], image_t *image) {
  FILE *file = fopen(path, "r");

  getc(file);
  getc(file);

  if (fscanf(file, "%d", &image->width) != 1) {
    fprintf(stderr, "erro ao ler largura da imagem %s\n", path);
  }
  if (fscanf(file, "%d", &image->height) != 1) {
    fprintf(stderr, "erro ao ler altura da imagem %s\n", path);
  }
  if (fscanf(file, "%*d") != 0) {
    fprintf(stderr, "erro ao ler 255 da imagem %s\n", path);
  }

  for (int y = 0; y < image->height; y++) {
    for (int x = 0; x < image->width; x++) {
      if (fscanf(file, "%d", &image->levels[y][x]) != 1) {
        fprintf(stderr, "erro ao ler nivel da imagem %s\n", path);
      }
    }
  }

  fclose(file);
}

void save_image(image_t *image, char path[]) {
  FILE *file = fopen(path, "w");

  fprintf(file, "P2\n%d %d\n255\n", image->width, image->height);

  for (int y = 0; y < image->height; y++) {
    for (int x = 0; x < image->width; x++) {
      fprintf(file, "%d\n", image->levels[y][x]);
    }
  }

  fclose(file);
}

void convert(double complex matrix[MAX_SIZE][MAX_SIZE], image_t *image) {
  for (int y = 0; y < image->height; y++) {
    for (int x = 0; x < image->width; x++) {
      image->levels[y][x] = fmin(fmax(0, round(creal(matrix[y][x]))), 255);
    }
  }
}
