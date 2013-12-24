#ifndef MFPT_INITIALIZERS_H
#define MFPT_INITIALIZERS_H

#include "data_structures.h"

KLP_MATRIX init_klp_matrix(int);
void free_klp_matrix(KLP_MATRIX);
void print_klp_matrix(KLP_MATRIX);
double** init_transition_matrix(int);
void free_transition_matrix(double**, int);
void print_transition_matrix(KLP_MATRIX, double**);

#endif