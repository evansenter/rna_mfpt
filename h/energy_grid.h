#ifndef MFPT_ENERGY_GRID_H
#define MFPT_ENERGY_GRID_H

#include "data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);  
int dgels_(char *t, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

#ifdef __cplusplus
}
#endif

double** convert_energy_grid_to_transition_matrix(KLP_MATRIX*, MFPT_PARAMETERS);
double compute_mfpt(KLP_MATRIX, MFPT_PARAMETERS, double**);
double* inverse(double*, int);
double* pseudoinverse(double*, int);
int number_of_permissible_single_bp_moves(int, int, KLP_MATRIX);
int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX, MFPT_PARAMETERS, int*, int*);
double transition_rate_from_probabilities(double, double, double);
double transition_rate_from_energies(double, double, double);
double transition_rate_from_probabilities_with_hastings(double, double, double, double);
double transition_rate_from_energies_with_hastings(double, double, double, double);

#endif
