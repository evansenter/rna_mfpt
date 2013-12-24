#ifndef MFPT_DATA_STRUCTURES_H
#define MFPT_DATA_STRUCTURES_H

typedef struct {
  int start_state;
  int end_state;
  int sequence_length;
  int bp_dist;
  short energy_based;
  short transition_matrix_input;
  short pseudoinverse;
  short single_bp_moves_only;
  short hastings;
  short verbose;
  double additive_epsilon;
  double distributed_epsilon;
} MFPT_PARAMETERS;

typedef struct {
  int* k;
  int* l;
  double* p;
  int length;
} KLP_MATRIX;

#endif
