#ifndef IMAGE_H
#define IMAGE_H

#include "fourier.h"

typedef struct {
  int levels[MAX_SIZE][MAX_SIZE];
  int width;
  int height;
} image_t;

void load_image(char path[], image_t *image);
void save_image(image_t *image, char path[]);

void convert(double complex matrix[MAX_SIZE][MAX_SIZE], image_t *image);

#endif
