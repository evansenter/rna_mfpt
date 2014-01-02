#ifndef MFPT_DATA_STRUCTURES_H
#define MFPT_DATA_STRUCTURES_H

typedef struct {
  int start_state;
  int end_state;
  int max_dist;
  int bp_dist;
  double epsilon;
  short energy_based;
  short transition_matrix_input;
  short pseudoinverse;
  short fully_connected;
  short single_bp_moves_only;
  short hastings;
  short radial_probability;
  short rate_matrix;
  short verbose;
} MFPT_PARAMS;

typedef struct {
  int* k;
  int* l;
  double* p;
  int length;
} KLP_MATRIX;

typedef double (*transition_probability)(const KLP_MATRIX, const double*, int, int, short);

#endif
