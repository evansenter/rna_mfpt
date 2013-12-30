#ifndef MFPT_INITIALIZERS_H
#define MFPT_INITIALIZERS_H

#include "data_structures.h"

#define ROW_ORDER(matrix, i, j, n) ((matrix)[((i) * (n)) + (j)])
#define COL_ORDER(matrix, i, j, n) ((matrix)[((j) * (n)) + (i)])

KLP_MATRIX init_klp_matrix(int);
void free_klp_matrix(KLP_MATRIX);
void print_klp_matrix(KLP_MATRIX);
double* init_transition_matrix(int);
double* transpose_matrix(double*, int);
void free_transition_matrix(double*);
void print_transition_matrix(KLP_MATRIX, double*);

#endif