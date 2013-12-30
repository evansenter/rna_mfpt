#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "constants.h"
#include "initializers.h"
#include "energy_grid.h"

#define DEBUG 1;

double* convert_energy_grid_to_transition_matrix(KLP_MATRIX* klp_matrix, MFPT_PARAMETERS* parameters) {
  int i, j, start_index, end_index, resolved;
  double row_sum;
  double* number_of_adjacent_moves;
  double* transition_probabilities;
  
  // This is just to get gcc to shut up, it's initialized with the proper (klp_matrix->length is increased in extend_klp_matrix_to_all_possible_positions) below.
  number_of_adjacent_moves = calloc(klp_matrix->length, sizeof(double));
  
  if (parameters->single_bp_moves_only) {
    if (!parameters->bp_dist) {
      resolved = find_start_and_end_positions_in_klp_matrix(klp_matrix, *parameters, &start_index, &end_index);
      set_bp_dist_from_start_and_end_positions(*klp_matrix, parameters, start_index, end_index, resolved);
    }
    
    if (parameters->sequence_length) {
      extend_klp_matrix_to_all_possible_positions(klp_matrix, *parameters);
    }
    
    number_of_adjacent_moves = realloc(number_of_adjacent_moves, klp_matrix->length * sizeof(double));
    for (i = 0; i < klp_matrix->length; ++i) {
      number_of_adjacent_moves[i] = (double)number_of_permissible_single_bp_moves(klp_matrix, i);
    }
    
    #ifdef DEBUG
    printf("\nFull dataset:\n");
    for (i = 0; i < klp_matrix->length; ++i) {
      printf("%d\t%d\t%f\t%d possible moves\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i], (int)number_of_adjacent_moves[i]);
    }
    #endif
  }
  
  transition_probabilities = init_transition_matrix(klp_matrix->length);
    
  for (i = 0; i < klp_matrix->length; ++i) {
    row_sum = 0.;
    
    for (j = 0; j < klp_matrix->length; ++j) {
      if (i != j) {        
        if (parameters->single_bp_moves_only) {
          if ((int)abs(klp_matrix->k[i] - klp_matrix->k[j]) == 1 && (int)abs(klp_matrix->l[i] - klp_matrix->l[j]) == 1) {
            // The problem here is that the position we're transitioning to has to be non-zero. So what happens is that if the starting position has a 0 probability
            // and it's adjacent positions have zero probabilities as well, this code won't do shit. So what we really need to do is allow arbitrary transitions on the
            // grid, which means that we no longe need to use the Hastings adjustment for detailed balance. In this way, it's always possible to leave a 0-probability
            // position, but I'm not sure it makes complete sense.
            printf("%d, %d => %d, %d\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->k[j], klp_matrix->l[j]);
            
            if (parameters->sequence_length && klp_matrix->p[i] == 0 && klp_matrix->p[j] > 0) {
              printf("Transitioning from p(%d, %d) = 0 to p(%d, %d) = %f\n", klp_matrix->k[i], klp_matrix->l[i], klp_matrix->k[j], klp_matrix->l[j], klp_matrix->p[j]);
              radial_probability(klp_matrix, i, parameters->sequence_length);
            }
            
            if (parameters->hastings) {
              if (parameters->energy_based) {
                ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                    transition_rate_from_energies_with_hastings(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i], number_of_adjacent_moves[j]);
              } else {
                ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                    transition_rate_from_probabilities_with_hastings(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i], number_of_adjacent_moves[j]);
              }
            } else {
              if (parameters->energy_based) {
                ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                    transition_rate_from_energies(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i]);
              } else {
                ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                    transition_rate_from_probabilities(klp_matrix->p[i], klp_matrix->p[j], number_of_adjacent_moves[i]);
              }
            }
          }
        } else {
          if (parameters->energy_based) {
            ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                transition_rate_from_energies(klp_matrix->p[i], klp_matrix->p[j], (double)(klp_matrix->length - 1));
          } else {
            ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) = \
                transition_rate_from_probabilities(klp_matrix->p[i], klp_matrix->p[j], (double)(klp_matrix->length - 1));
          }
        }
        
        row_sum += ROW_ORDER(transition_probabilities, i, j, klp_matrix->length);
      }
    }
    
    ROW_ORDER(transition_probabilities, i, i, klp_matrix->length) = 1 - row_sum;
  }
  
  return transition_probabilities;
}

double compute_mfpt(KLP_MATRIX* klp_matrix, MFPT_PARAMETERS parameters, double* transition_probabilities) {
  int i, j, x, y, start_index, end_index, resolved, inversion_matrix_row_length = klp_matrix->length - 1;
  double mfpt_from_start;
  
  printf("\nFull dataset:\n");
  for (i = 0; i < klp_matrix->length; ++i) {
    printf("%d\t%d\t%d\t%f\n", i, klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i]);
  }
  
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
  
  double* mfpt             = calloc(inversion_matrix_row_length, sizeof(double));
  double* inversion_matrix = malloc((int)pow((double)inversion_matrix_row_length, 2.) * sizeof(double));
  
  for (i = 0; i < klp_matrix->length; ++i) {
    for (j = 0; j < klp_matrix->length; ++j) {
      if (i != end_index && j != end_index) {
        x = (i > end_index ? i - 1 : i);
        y = (j > end_index ? j - 1 : j);
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (end_index) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversion_matrix and i, j are used for indexing into transition_probabilities.
        inversion_matrix[x * inversion_matrix_row_length + y] = \
            (i == j ? 1 - ROW_ORDER(transition_probabilities, i, j, klp_matrix->length) : -ROW_ORDER(transition_probabilities, i, j, klp_matrix->length));
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
  printf("dgetrf_(&size, &size, a, &size, ipiv, &info)\ninfo:\t%d\n", info);
  #endif
  dgetri_(&size, a, &size, ipiv, work, &lwork, &info);
  #ifdef DEBUG
  printf("dgetri_(&size, a, &size, ipiv, work, &lwork, &info)\ninfo:\t%d\n", info);
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
  return (b);
}

int find_start_and_end_positions_in_klp_matrix(KLP_MATRIX* klp_matrix, MFPT_PARAMETERS parameters, int* start_index, int* end_index) {
  int i, resolved = 0;
  
  if (parameters.start_state == -1) {
    for (i = 0, *start_index = -1; i < klp_matrix->length && *start_index < 0; ++i) {
      if (klp_matrix->k[i] == 0) {
        *start_index = i;
        resolved++;
      }
    }
  } else {
    *start_index = parameters.start_state;
    resolved++;
  }
  
  if (parameters.end_state == -1) {
    for (i = 0, *end_index = -1; i < klp_matrix->length && *end_index < 0; ++i) {
      if (klp_matrix->l[i] == 0) {
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

void set_bp_dist_from_start_and_end_positions(KLP_MATRIX klp_matrix, MFPT_PARAMETERS* parameters, int start_index, int end_index, int resolved) {
  int distance_from_start, distance_from_end;
  
  distance_from_start = distance_from_end = -1;
  
  if (start_index > 0) {
    distance_from_start = klp_matrix.l[start_index];
  }
  
  if (end_index > 0) {
    distance_from_end = klp_matrix.k[end_index];
  }
  
  if (distance_from_start == distance_from_end && resolved) {
    parameters->bp_dist = distance_from_start;
  } else if (distance_from_start >= 0 && distance_from_end == -1) {
    parameters->bp_dist = distance_from_start;
  } else if (distance_from_end >= 0 && distance_from_start == -1) {
    parameters->bp_dist = distance_from_end;
  } else {
    fprintf(stderr, "Can't infer the input structure distances for the energy grid. We found (0, %d) and (%d, 0). Consider using the -d flag to manually set the base pair distance between the two structures.\n", distance_from_end, distance_from_start);
    printf("-3\n");
    exit(0);
  }
}

void extend_klp_matrix_to_all_possible_positions(KLP_MATRIX* klp_matrix, MFPT_PARAMETERS parameters) {
  int i, j, m, position_in_input_data, pointer, valid_positions = 0;
  
  #ifdef DEBUG
  printf("\nAccessible positions (top-left is [0, 0]):\n");
  #endif

  for (i = 0; i <= parameters.sequence_length; ++i) {
    for (j = 0; j <= parameters.sequence_length; ++j) {
      if (
        i + j >= parameters.bp_dist &&
        i + parameters.bp_dist >= j &&
        j + parameters.bp_dist >= i &&
        (i + j) % 2 == parameters.bp_dist % 2
      ) {
        #ifdef DEBUG
        position_in_input_data = -1;
        for (m = 0; m < klp_matrix->length && position_in_input_data == -1; ++m) {
          if (klp_matrix->k[m] == i && klp_matrix->l[m] == j) {
            position_in_input_data = m;
          }
        }
        printf(position_in_input_data == -1 ? "X" : "O");
        #endif
        valid_positions++;
      } else {
        #ifdef DEBUG
        printf(" ");
        #endif
      }
    }
    #ifdef DEBUG
    printf("\n");
    #endif
  }

  klp_matrix->k = realloc(klp_matrix->k, valid_positions * sizeof(int));
  klp_matrix->l = realloc(klp_matrix->l, valid_positions * sizeof(int));
  klp_matrix->p = realloc(klp_matrix->p, valid_positions * sizeof(double));

  pointer = klp_matrix->length;

  #ifdef DEBUG
  printf("\nInput dataset:\n");
  for (i = 0; i < klp_matrix->length; ++i) {
    printf("%d\t%d\t%d\t%f\n", i, klp_matrix->k[i], klp_matrix->l[i], klp_matrix->p[i]);
  }
  #endif

  for (i = 0; i <= parameters.sequence_length; ++i) {
    for (j = 0; j <= parameters.sequence_length; ++j) {
      if (
        i + j >= parameters.bp_dist &&
        i + parameters.bp_dist >= j &&
        j + parameters.bp_dist >= i &&
        (i + j) % 2 == parameters.bp_dist % 2
      ) {
        position_in_input_data = -1;
      
        for (m = 0; m < klp_matrix->length && position_in_input_data == -1; ++m) {
          if (klp_matrix->k[m] == i && klp_matrix->l[m] == j) {
            position_in_input_data = m;
          }
        }
      
        if (position_in_input_data < 0) {
          klp_matrix->k[pointer] = i;
          klp_matrix->l[pointer] = j;
          klp_matrix->p[pointer] = 0.;
          pointer++;
        }
      }
    }
  }

  klp_matrix->length = valid_positions;
}

int number_of_permissible_single_bp_moves(KLP_MATRIX* klp_matrix, int i) {
  int j, x, y, a, b, num_moves = 0;
  
  x = klp_matrix->k[i];
  y = klp_matrix->l[i];
  
  for (j = 0; j < klp_matrix->length; ++j) {
    a = klp_matrix->k[j];
    b = klp_matrix->l[j];
    
    if (
      // Because N(x, y) is restricted to entries in *k and *l, we *assume* the input data satisfies the triangle inequality and bounds.
      (int)abs(x - a) == 1 && (int)abs(y - b) == 1
    ) {
      num_moves++;
    }
  }
  
  return num_moves;
}

double radial_probability(KLP_MATRIX* klp_matrix, int klp_index, int row_size) {
  int distance, i, j, k, x, y, a, b, found_klp_index, found_nonzero_probability = 0;
  double probability_sum = 0;
  
  x = klp_matrix->k[klp_index];
  y = klp_matrix->l[klp_index];
  
  for (distance = 1; distance < row_size && !found_nonzero_probability; ++distance) {
    for (i = 0; i <= distance; ++i) {
      for (j = 0; j <= distance; ++j) {
        // If the position is on the perimeter...
        if (!(i > 0 && i < distance && j > 0 && j < distance)) {
          a = x + i * 2 - distance;
          b = y + j * 2 - distance;
          
          // Look for the position in the stationary probability data structure (assumes klp_matrix is populated with exactly all possible positions)
          found_klp_index = -1;
          for (k = 0; k < klp_matrix->length && found_klp_index < 0; ++k) {
            if (klp_matrix->k[k] == a && klp_matrix->l[k] == b) {
              found_klp_index = 1;
              
              if (klp_matrix->p[k] > 0) {
                printf("%d, %d: (distance %d (%d, %d) = %f)\n", x, y, distance, a, b, klp_matrix->p[k]);
                
                probability_sum          += klp_matrix->p[k];
                found_nonzero_probability = 1;
              }
            }
          }
        }
      }
    }
  }
  
  printf("probability_sum: %f\n", probability_sum);
  
  return probability_sum;
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
