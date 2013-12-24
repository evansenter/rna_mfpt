#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "constants.h"
#include "mfpt_params.h"
#include "mfpt_parser.h"
#include "mfpt_initializers.h"
#include "mfpt_energy_grid.h"

int main(int argc, char* argv[]) {
  int i, line_count, row_length;
  double mfpt;
  double** transition_matrix;
  MFPT_PARAMETERS parameters;
  KLP_MATRIX klp_matrix;
  
  parameters = parse_mfpt_args(argc, argv);
  line_count = count_lines(argv[argc - 1]);
  
  if (!line_count) {
    fprintf(stderr, "%s appears to have no data.\n", argv[argc - 1]);
    return 0;
  }
  
  klp_matrix = init_klp_matrix(line_count);
  populate_arrays(argv[argc - 1], klp_matrix, parameters);
  
  #ifdef SUPER_HEAVY_DEBUG
    print_klp_matrix(klp_matrix);
  #endif
  
  // We already have a transition matrix, this is the easy case. Just need to find MFPT.
  if (parameters.transition_matrix_input) {
    // We need to infer the dimensions of the transition matrix.
    row_length = 0;
    for (i = 0; i < line_count; ++i) {
      row_length = klp_matrix.k[i] > row_length ? klp_matrix.k[i] : row_length;
      row_length = klp_matrix.l[i] > row_length ? klp_matrix.l[i] : row_length;
    }
    row_length++;
    
    transition_matrix = init_transition_matrix(row_length);
    for (i = 0; i < line_count; ++i) {
      transition_matrix[klp_matrix.k[i]][klp_matrix.l[i]] = klp_matrix.p[i];
    }
  // We have an energy grid, this requires converting the energy grid into a transition matrix data structure before finding MFPT.
  } else {
    row_length        = line_count;
    transition_matrix = convert_energy_grid_to_transition_matrix(&klp_matrix, parameters);
  }
  
  #if SUPER_HEAVY_DEBUG
    print_transition_matrix(klp_matrix, transition_matrix);
  #endif
  
  mfpt = compute_mfpt(klp_matrix, parameters, transition_matrix);
  printf("%+.8f\n", mfpt);
  
  free_klp_matrix(klp_matrix);
  free_transition_matrix(transition_matrix, row_length);
    
  return 0;
}
