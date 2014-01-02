#ifndef MFPT_ENERGY_GRID_H
#define MFPT_ENERGY_GRID_H

#include "data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
int dgels_(char* t, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info);

#ifdef __cplusplus
}
#endif

double* convert_klp_matrix_to_transition_matrix(KLP_MATRIX*, MFPT_PARAMS*);
double compute_mfpt(KLP_MATRIX*, const MFPT_PARAMS, const double*);
double* inverse(double*, int);
double* pseudoinverse(double*, int);
int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX*, MFPT_PARAMS*);
void set_bp_dist_from_start_and_end_positions(const KLP_MATRIX, MFPT_PARAMS*, int);
int extend_klp_matrix_to_all_possible_positions(KLP_MATRIX*, const MFPT_PARAMS);
void populate_remaining_probabilities_in_klp_matrix(KLP_MATRIX*, const MFPT_PARAMS, int);
double* populate_number_of_adjacent_moves(const KLP_MATRIX, const MFPT_PARAMS);
int number_of_permissible_single_bp_moves(const KLP_MATRIX, int);
double* populate_transition_matrix_from_stationary_matrix(const KLP_MATRIX, const MFPT_PARAMS, const double*, transition_probability);
double radial_probability(const KLP_MATRIX, int, int);
double transition_rate_from_probabilities(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_probabilities_with_hastings(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_energies_with_hastings(const KLP_MATRIX, const double*, int, int, short);
double transition_rate_from_radial_probability(const KLP_MATRIX, const double*, int, int, short);

#endif
