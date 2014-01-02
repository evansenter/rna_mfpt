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
  printf("\nk/l/p matrix (index, k, l, p):\n");

  for (i = 0; i < klp_matrix.length; ++i) {
    printf("%d\t%d\t%d\t%+.8f\n", i, klp_matrix.k[i], klp_matrix.l[i], klp_matrix.p[i]);
  }
  
  printf("\n");
}

double* init_transition_matrix(int length) {
  return calloc(length * length, sizeof(double));;
}

double* transpose_matrix(double* matrix, int length) {
  int i, j;
  double* transposed_matrix = init_transition_matrix(length);
  
  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) {
      COL_ORDER(transposed_matrix, i, j, length) = ROW_ORDER(matrix, i, j, length);
    }
  }
  
  free_transition_matrix(matrix);
  return transposed_matrix;
}

void free_transition_matrix(double* transition_matrix) {
  free(transition_matrix);
}

void print_transition_matrix(KLP_MATRIX klp_matrix, double* transition_matrix) {
  int i, j;
  printf("Transition matrix:\n");
  printf("(from)\t(to)\tp(to | from)\n");

  for (i = 0; i < klp_matrix.length; ++i) {
    for (j = 0; j < klp_matrix.length; ++j) {
      printf(
        "(%d, %d)\t=>\t(%d, %d)\t%.8f\n", 
        klp_matrix.k[i], 
        klp_matrix.l[i], 
        klp_matrix.k[j], 
        klp_matrix.l[j], 
        ROW_ORDER(transition_matrix, i, j, klp_matrix.length)
      );
    }
  }
  
  printf("\n");
}