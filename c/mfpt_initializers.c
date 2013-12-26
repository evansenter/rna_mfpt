#include <stdio.h>
#include <stdlib.h>
#include "initializers.h"

KLP_MATRIX init_klp_matrix(int length) {
  KLP_MATRIX klp_matrix = {
    .k      = (int*)malloc(length * sizeof(int)),
    .l      = (int*)malloc(length * sizeof(int)),
    .p      = (double*)malloc(length * sizeof(double)),
    .length = length
  };
  return klp_matrix;
}

void free_klp_matrix(KLP_MATRIX klp_matrix) {
  free(klp_matrix.k);
  free(klp_matrix.l);
  free(klp_matrix.p);
}

void print_klp_matrix(KLP_MATRIX klp_matrix) {
  int i;
  printf("\nk/l/p matrix:\n");

  for (i = 0; i < klp_matrix.length; ++i) {
    printf("%d\t%d\t%.8f\n", klp_matrix.k[i], klp_matrix.l[i], klp_matrix.p[i]);
  }
}

double** init_transition_matrix(int length) {
  int i;
  double** transition_matrix;
  transition_matrix = malloc(length * sizeof(double*));

  for (i = 0; i < length; ++i) {
    transition_matrix[i] = calloc(length, sizeof(double));
  }

  return transition_matrix;
}

void free_transition_matrix(double** transition_matrix, int length) {
  int i;

  for (i = 0; i < length; ++i) {
    free(transition_matrix[i]);
  }

  free(transition_matrix);
}

void print_transition_matrix(KLP_MATRIX klp_matrix, double** transition_matrix) {
  int i;
  printf("Transition matrix:\n");
  printf("(x)\t(y)\tp(x, y)\n");

  for (i = 0; i < klp_matrix.length; ++i) {
    printf("(%d)\t=>\t(%d)\t%.8f\n", klp_matrix.k[i], klp_matrix.l[i], transition_matrix[klp_matrix.k[i]][klp_matrix.l[i]]);
  }
}