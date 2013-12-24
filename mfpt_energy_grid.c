#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "mfpt_params.h"
#include "mfpt_energy_grid.h"
#include "mfpt_initializers.h"

double** convert_energy_grid_to_transition_matrix(KLP_MATRIX* klp_matrix, MFPT_PARAMETERS parameters) {
  int i, j, m, bp_distance, input_data_index, start_index, end_index, resolved, distance_from_start = -1, distance_from_end = -1, pointer = 0, validPositions = 0;
  int* old_k;
  int* old_l;
  double row_sum, epsilon;
  double* old_p;
  double* number_of_adjacent_moves;
  double** transition_probabilities;
  
  number_of_adjacent_moves = malloc(klp_matrix->length * sizeof(double));
  
  if (parameters.single_bp_moves_only) {
    if (parameters.bp_dist) {
      bp_distance = parameters.bp_dist;
    } else {
      resolved = find_start_and_end_positions_in_klp_matrix(*klp_matrix, parameters, &start_index, &end_index);
    
      if (start_index > 0) {
        distance_from_start = klp_matrix->l[start_index];
      }
      
      if (end_index > 0) {
        distance_from_end = klp_matrix->k[end_index];
      }
    
      if (distance_from_start == distance_from_end && resolved) {
        bp_distance = distance_from_start;
      } else if (distance_from_start >= 0 && distance_from_end == -1) {
        bp_distance = distance_from_start;
      } else if (distance_from_end >= 0 && distance_from_start == -1) {
        bp_distance = distance_from_end;
      } else {
        fprintf(stderr, "Can't infer the input structure distances for the energy grid. We found (0, %d) and (%d, 0). Consider using the -D flag to manually set the base pair distance between the two structures.\n", distance_from_end, distance_from_start);
        printf("-3\n");
        exit(0);
      }
    }
    
    if (parameters.sequence_length) {
      #ifdef DEBUG
        printf("\nAccessible positions (top-left is [0, 0]):\n");
      #endif
      
      for (i = 0; i <= parameters.sequence_length; ++i) {
        for (j = 0; j <= parameters.sequence_length; ++j) {
          if (
            i + j >= bp_distance &&
            i + bp_distance >= j &&
            j + bp_distance >= i &&
            (i + j) % 2 == bp_distance % 2
          ) {
            #ifdef DEBUG
              printf(" ");
            #endif
            validPositions++;
          } else {
            #ifdef DEBUG
              printf("X");
            #endif
          }
        }
        #ifdef DEBUG
          printf("\n");      
        #endif
      }
      
      old_k = malloc(klp_matrix->length * sizeof(int));
      old_l = malloc(klp_matrix->length * sizeof(int));
      old_p = malloc(klp_matrix->length * sizeof(double));
      
      for (i = 0; i <= klp_matrix->length; ++i) {
        old_k[i] = klp_matrix->k[i];
        old_l[i] = klp_matrix->l[i];
        old_p[i] = klp_matrix->p[i];
      }
      
      klp_matrix->k =    (int*)realloc(klp_matrix->k, validPositions * sizeof(int));
      klp_matrix->l =    (int*)realloc(klp_matrix->l, validPositions * sizeof(int));
      klp_matrix->p = (double*)realloc(klp_matrix->p, validPositions * sizeof(double));
      
      epsilon = (parameters.additive_epsilon ? parameters.additive_epsilon : parameters.distributed_epsilon / validPositions);
      
      for (i = 0; i <= parameters.sequence_length; ++i) {
        for (j = 0; j <= parameters.sequence_length; ++j) {
          if (
            i + j >= bp_distance &&
            i + bp_distance >= j &&
            j + bp_distance >= i &&
            (i + j) % 2 == bp_distance % 2
          ) {
            input_data_index = -1;
            
            for (m = 0; m < klp_matrix->length && input_data_index == -1; ++m) {
              if (old_k[m] == i && old_l[m] == j) {
                input_data_index = m;
              }
            }
            
            klp_matrix->k[pointer] = i;
            klp_matrix->l[pointer] = j;
            klp_matrix->p[pointer] = (input_data_index == -1 ? 0. : old_p[input_data_index]) + epsilon;
            
            if (!parameters.energy_based) {
              klp_matrix->p[pointer] /= 1. + epsilon * validPositions;
            }
            
            pointer++;
          }
        }
      }
      
      free(old_k);
      free(old_l);
      free(old_p);
      
      klp_matrix->length = validPositions;
    }
    
    #ifdef DEBUG
      printf("\nFull dataset:\n");
    #endif
  
    for (i = 0; i < klp_matrix->length; ++i) {
      number_of_adjacent_moves[i] = (double)number_of_permissible_single_bp_moves(klp_matrix->k[i], klp_matrix->l[i], *klp_matrix);
    
      #ifdef DEBUG
        printf("%d\t%d\t%.15f\t%d possible moves\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i], (int)number_of_adjacent_moves[i]);
      #endif
    }
  }
  
  transition_probabilities = init_transition_matrix(klp_matrix->length);
  
  for (i = 0; i < klp_matrix->length; ++i) {
    row_sum = 0.;
      
    for (j = 0; j < klp_matrix->length; ++j) {
      if (i != j) {
        if (parameters.single_bp_moves_only) {
          if ((int)abs(klp_matrix->k[i] - klp_matrix->k[j]) == 1 && (int)abs(klp_matrix->l[i] - klp_matrix->l[j]) == 1) {
            if (parameters.hastings) {
              if (parameters.energy_based) {
                transition_probabilities[i][j] = transition_rate_from_energies_with_hastings(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i], number_of_adjacent_moves[j]);
              } else {
                transition_probabilities[i][j] = transition_rate_from_probabilities_with_hastings(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i], number_of_adjacent_moves[j]);
              }
            } else {
              if (parameters.energy_based) {
                transition_probabilities[i][j] = transition_rate_from_energies(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i]);
              } else {
                transition_probabilities[i][j] = transition_rate_from_probabilities(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i]);
              }
            }
          }
        } else {
          if (parameters.energy_based) {
            transition_probabilities[i][j] = transition_rate_from_energies(klp_matrix->p[i], klp_matrix->p[j], (double)(klp_matrix->length - 1));
          } else {
            transition_probabilities[i][j] = transition_rate_from_probabilities(klp_matrix->p[i], klp_matrix->p[j], (double)(klp_matrix->length - 1));
          }
        }
        
        row_sum += transition_probabilities[i][j];
      }
    }
    
    transition_probabilities[i][i] = 1 - row_sum;
  }
  
  return transition_probabilities;
}

double compute_mfpt(KLP_MATRIX klp_matrix, MFPT_PARAMETERS parameters, double **transition_probabilities) {
  int i, j, x, y, start_index, end_index, resolved, inversion_matrix_row_length = klp_matrix.length - 1;
  double mfpt_from_start;
  
  resolved = find_start_and_end_positions_in_klp_matrix(klp_matrix, parameters, &start_index, &end_index);
  
  if (resolved != 2) {
    if (start_index < 0) {
      #ifdef DEBUG
        fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
      #endif
      return -1;
    }
  
    if (end_index < 0) {
      #ifdef DEBUG
        fprintf(stderr, "We can not find any position in the energy grid correspondent to the stopping state.\n");
      #endif
      return -2;
    }
  }
  
  // If start_index > end_index, we need to shift to the left by one because the end_index row / column is being removed.
  if (start_index > end_index) {
    start_index--;
  }
  
  double *mfpt             = calloc(inversion_matrix_row_length, sizeof(double));
  double *inversion_matrix = malloc((int)pow((double)inversion_matrix_row_length, 2.) * sizeof(double));
  
  for (i = 0; i < klp_matrix.length; ++i) {
    for (j = 0; j < klp_matrix.length; ++j) { 
      if (i != end_index && j != end_index) {
        x = (i > end_index ? i - 1 : i);
        y = (j > end_index ? j - 1 : j);
        
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (end_index) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversion_matrix and i, j are used for indexing into transition_probabilities.
        inversion_matrix[x * inversion_matrix_row_length + y] = (i == j ? 1 - transition_probabilities[i][j] : -transition_probabilities[i][j]);
      }
    }
  }
  
  inversion_matrix = parameters.pseudoinverse ? pseudoinverse(inversion_matrix, inversion_matrix_row_length) : inverse(inversion_matrix, inversion_matrix_row_length);
  
  #ifdef DEBUG
    printf("\nMFPT values for indices into the full dataset:\n");
  #endif
  
  for (i = 0; i < inversion_matrix_row_length; ++i) {
    for (j = 0; j < inversion_matrix_row_length; ++j) {
      mfpt[i] += inversion_matrix[i * inversion_matrix_row_length + j];
    }
    
    #ifdef DEBUG
      // The business with this i < end_index stuff is inorder to ensure that the output MFPT debug indices are representative of the input data.
      printf("%d:\t%f\n", i < end_index ? i : i + 1, mfpt[i]);
    #endif
  }
    
  mfpt_from_start = mfpt[start_index];
  free(mfpt);
  free(inversion_matrix);
  
  return mfpt_from_start;
}

double* inverse(double* a, int size) {
  // http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
  int* ipiv    = malloc((size + 1) * sizeof(int));
  int lwork    = size * size;
  double* work = malloc(lwork * sizeof(double));
  int info;

  dgetrf_(&size, &size, a, &size, ipiv, &info);
  
  #ifdef DEBUG
    printf("dgetrf_(&size, &size, a, &size, ipiv, &info) info:\t%d\n", info);
  #endif
  
  dgetri_(&size, a, &size, ipiv, work, &lwork, &info);
  
  #ifdef DEBUG
    printf("dgetri_(&size, a, &size, ipiv, work, &lwork, &info) info:\t%d\n", info);
  #endif

  free(ipiv);
  free(work);
  
  return a;
}

double* pseudoinverse(double* a, int size) {
  // Least-squares fit solution to B - Ax, where (in this case) A is square and B is the identity matrix.
  char trans;
  int i, m, n, nrhs, lda, ldb, lwork, info;

  trans = 'N';
  m     = size;
  n     = m;
  nrhs  = m;
  lda   = m;
  ldb   = m;
  
  double* b = calloc(ldb * nrhs, sizeof(double));
  for (i = 0; i < ldb; ++i) {
    b[i * nrhs + i] = 1.;
  }
  
  #ifdef SUPER_HEAVY_DEBUG
    printf("dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info)\n\n");
  #endif
    
  lwork        = -1;
  double* work = malloc(MAX(1, lwork) * sizeof(double));
  
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  
  lwork = (int)work[0];
  
  #ifdef SUPER_HEAVY_DEBUG
    printf("workspace query (lwork):\t%d\n", lwork);
  #endif
    
  free(work);
  
  work = malloc(MAX(1, lwork) * sizeof(double));
  
  dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  
  #ifdef DEBUG
    printf("info:\t%d\n", info);
  #endif
  
  free(work);
  
  return(b);
}

int number_of_permissible_single_bp_moves(int x, int y, KLP_MATRIX klp_matrix) {
  int j, a, b, num_moves = 0;
  
  for (j = 0; j < klp_matrix.length; ++j) {
    a = klp_matrix.k[j];
    b = klp_matrix.l[j];
    
    if (
      // Because N(x, y) is restricted to entries in *k and *l, we *assume* the input data satisfies the triangle inequality and bounds.
      (int)abs(x - a) == 1 && (int)abs(y - b) == 1
    ) {
      num_moves++;
    }
  }
  
  return num_moves;
}

int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX klp_matrix, MFPT_PARAMETERS parameters, int* start_index, int* end_index) {
  int i, resolved = 0;
  
  if (parameters.start_state == -1) {
    for (i = 0, *start_index = -1; i < klp_matrix.length && *start_index < 0; ++i) {
      if (klp_matrix.k[i] == 0) {
        *start_index = i;
        resolved++;
      }
    }
  } else {
    *start_index = parameters.start_state;
    resolved++;
  }
  
  if (parameters.end_state == -1) {
    for (i = 0, *end_index = -1; i < klp_matrix.length && *end_index < 0; ++i) {
      if (klp_matrix.l[i] == 0) {
        *end_index = i;
        resolved++;
      }
    }
  } else {
    *end_index = parameters.end_state;
    resolved++;
  }
  
  #ifdef DEBUG
    printf("\nstart_index:\t%d\n", *start_index);
    printf("end_index:\t%d\n", *end_index);
    printf("resolved:\t%d\n", resolved);
  #endif
    
  return resolved;
}

double transition_rate_from_probabilities(double from, double to, double num_from) {
  return MIN(1., to / from) / num_from;
}

double transition_rate_from_energies(double from, double to, double num_from) {
  return MIN(1., exp(-(to - from) / RT)) / num_from;
}

double transition_rate_from_probabilities_with_hastings(double from, double to, double num_from, double num_to) {
  return MIN(1., (num_from / num_to) * (to / from)) / num_from;
}

double transition_rate_from_energies_with_hastings(double from, double to, double num_from, double num_to) {
  return MIN(1., (num_from / num_to) * exp(-(to - from) / RT)) / num_from;
}