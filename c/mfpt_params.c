#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "params.h"

MFPT_PARAMETERS init_mfpt_params() {
  MFPT_PARAMETERS parameters = {
    .start_state             = -1,
    .end_state               = -1,
    .sequence_length         = 0,
    .bp_dist                 = 0,
    .energy_based            = 0,
    .transition_matrix_input = 0,
    .pseudoinverse           = 0,
    .single_bp_moves_only    = 0,
    .hastings                = 0,
    .verbose                 = 0,
    .additive_epsilon        = 0,
    .distributed_epsilon     = 0
  };
  return parameters;
}

MFPT_PARAMETERS parse_mfpt_args(int argc, char* argv[]) {
  int c;
  MFPT_PARAMETERS parameters;
  parameters = init_mfpt_params();
  
  while ((c = getopt(argc, argv, "EeTtPpXxHhVvA:a:Z:z:N:n:D:d:O:o:Q:q:")) != -1) {
    switch (c) {
      case 'E':
      case 'e':
        parameters.energy_based = 1;
        break;
        
      case 'T':
      case 't':
        parameters.transition_matrix_input = 1;
        break;
        
      case 'P':
      case 'p':
        parameters.pseudoinverse = 1;
        break;
        
      case 'X':
      case 'x':
        parameters.single_bp_moves_only = 1;
        break;
        
      case 'H':
      case 'h':
        parameters.hastings = 1;
        break;
        
      case 'V':
      case 'v':
        parameters.verbose = 1;
        break;
        
      case 'A':
      case 'a':
        if (!sscanf(optarg, "%d", &(parameters.start_state))) {
          mfpt_usage();
        } else if (parameters.start_state < 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'Z':
      case 'z':
        if (!sscanf(optarg, "%d", &(parameters.end_state))) {
          mfpt_usage();
        } else if (parameters.end_state < 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'N':
      case 'n':
        if (!sscanf(optarg, "%d", &(parameters.sequence_length))) {
          mfpt_usage();
        } else if (parameters.sequence_length <= 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'D':
      case 'd':
        if (!sscanf(optarg, "%d", &(parameters.bp_dist))) {
          mfpt_usage();
        } else if (parameters.bp_dist <= 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'O':
      case 'o':
        if (!sscanf(optarg, "%lf", &(parameters.additive_epsilon))) {
          mfpt_usage();
        } else if (parameters.additive_epsilon <= 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'Q':
      case 'q':
        if (!sscanf(optarg, "%lf", &(parameters.distributed_epsilon))) {
          mfpt_usage();
        } else if (parameters.distributed_epsilon <= 0) {
          mfpt_usage();
        }
        
        break;
        
      case '?':
        switch (optopt) {
          case 'A':
          case 'a':
          case 'Z':
          case 'z':
          case 'N':
          case 'n':
          case 'D':
          case 'd':
          case 'O':
          case 'o':
          case 'Q':
          case 'q':
            fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            break;
            
          default:
            if (isprint(optopt)) {
              fprintf(stderr, "Unknown option `-%c'.\n", optopt);
            } else {
              fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            }
        }
        
        mfpt_usage();
        
      default:
        mfpt_usage();
    }
  }
  
  if (optind + 1 != argc) {
    mfpt_usage();
  }
  
  if (parameters.verbose) {
    debug_mfpt_parameters(parameters);
  }
  
  if (mfpt_error_handling(parameters)) {
    mfpt_usage();
  }
  
  return parameters;
}

int mfpt_error_handling(MFPT_PARAMETERS parameters) {
  int error = 0;
  
  if (parameters.distributed_epsilon < 1e-15) {
    parameters.distributed_epsilon = 0;
  }
  
  if (parameters.additive_epsilon < 1e-15) {
    parameters.additive_epsilon = 0;
  }
  
  if (parameters.start_state == parameters.end_state && parameters.start_state >= 0) {
    fprintf(stderr, "Error: If the -a and -z flags are identical the MFPT is 0!\n");
    error++;
  }
  
  if (parameters.transition_matrix_input && (
        parameters.hastings ||
        parameters.sequence_length ||
        parameters.additive_epsilon ||
        parameters.single_bp_moves_only)) {
    fprintf(stderr, "Error: If the -t flag is provided, -h, -n, -o and -x are not permitted!\n");
    error++;
  }
  
  if (parameters.transition_matrix_input && !(parameters.start_state >= 0 && parameters.end_state >= 0)) {
    fprintf(stderr, "Error: If the -t flag is provided, -a and -z must be explicitly set!\n");
    error++;
  }
  
  if (parameters.single_bp_moves_only && parameters.energy_based && (parameters.sequence_length || parameters.additive_epsilon)) {
    fprintf(stderr, "Error: If the -x and -e flags are provided, -n and -o are not permitted!\n");
    error++;
  }
  
  if (parameters.hastings && !parameters.single_bp_moves_only) {
    fprintf(stderr, "Error: If the -h flag is provided, -x must be explicitly set!\n");
    error++;
  }
  
  if ((parameters.sequence_length && !(parameters.additive_epsilon || parameters.distributed_epsilon)) ||
      (!parameters.sequence_length && (parameters.additive_epsilon || parameters.distributed_epsilon))) {
    fprintf(stderr, "Error: The -n flag must be used with exactly one of the -o or -q flags!\n");
    error++;
  }
  
  if (parameters.additive_epsilon && parameters.distributed_epsilon) {
    fprintf(stderr, "Error: The -o and -q flags are not permitted together!\n");
    error++;
  }
  
  if (parameters.sequence_length && (parameters.energy_based || parameters.start_state >= 0 || parameters.end_state >= 0)) {
    fprintf(stderr, "Error: If the -n flag is provided, -e, -a and -z are not permitted!\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_mfpt_parameters(MFPT_PARAMETERS parameters) {
  char buffer[128];
  printf("(e) energy_based\t\t%s\n",          parameters.energy_based            ? "Yes" : "No");
  printf("(t) transition_matrix_input\t%s\n", parameters.transition_matrix_input ? "Yes" : "No");
  printf("(p) pseudoinverse\t\t%s\n",         parameters.pseudoinverse           ? "Yes" : "No");
  printf("(x) single_bp_moves_only\t%s\n",    parameters.single_bp_moves_only    ? "Yes" : "No");
  printf("(h) hastings\t\t\t%s\n",            parameters.hastings                ? "Yes" : "No");
  sprintf(buffer, "%d", parameters.sequence_length);
  printf("(n) sequence_length\t\t%s\n", parameters.sequence_length ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.bp_dist);
  printf("(d) bp_dist\t\t\t%s\n", parameters.bp_dist ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%.2e", parameters.additive_epsilon);
  printf("(o) additive_epsilon\t\t%s\n", parameters.additive_epsilon ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%.2e", parameters.distributed_epsilon);
  printf("(q) distributed_epsilon\t\t%s\n", parameters.distributed_epsilon ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.start_state);
  printf("(a) start_state\t\t\t%s\n", parameters.start_state >= 0 ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.end_state);
  printf("(z) end_state\t\t\t%s\n", parameters.end_state >= 0 ? buffer : "N/A");
  memset(buffer, ' ', 128 * sizeof(char));
}

void mfpt_usage() {
  fprintf(stderr, "RNAmfpt [options] input_csv\n\n");
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-A/a\tstart state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the first column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the start state.\n");
  fprintf(stderr, "-D/d\tstart/end distance, the default is disabled. When provided, indicates the base pair distance between the starting / ending structures. This flag is used in conjunction with the -Q flag, and is needed in cases when the base pair distance between the two structures can't be inferred from the input grid.\n");
  fprintf(stderr, "-E/e\tenergy-based transitions, the default is disabled. If this flag is provided, the transition from state a to b will be calculated as (min(1, exp(-(E_b - E_a) / RT) / n) rather than (min(1, p_b / p_a) / n).\n");
  fprintf(stderr, "-H/h\tHastings adjustment, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted (in the all-to-all transition case, N(X) / N(Y) == 1). Calculating N(X) and N(Y) will respect grid boundaries and the triangle equality, and the basepair distance between the two structures for kinetics is inferred from the energy grid.\n");
  fprintf(stderr, "-N/n\tsequence length, the default is disabled. This flag represents the sequence length of the sequence on which kinetics is being performed. It is used in conjunction with the -O flag to ensure that the graph is fully connected. If -O is not explicitly set, the -Q flag must be set.\n");
  fprintf(stderr, "-O/o\tepsilon, the default is disabled. This flag should be a %%f-parseable epsilon value added to *all* positions in the energy grid (including valid positions not present in the input), which is then renormalized and used to ensure that the graph is fully connected.\n");
  fprintf(stderr, "-P/p\tpseudoinverse, the default is disabled. If this flag is provided, the Moore-Penrose pseudoinverse is computed for the transition probability matrix, rather than the true inverse.\n")
  ;
  fprintf(stderr, "-Q/q\talternative epsilon, the default is disabled. If this flag is provided, each accessible position is increased by .01 / num_accessible_positions, and then renormalized to ensure that all valid positions have a small non-zero probability. This flag is mutually exclusive with -O.\n");
  fprintf(stderr, "-T/t\ttransition matrix input, the default is disabled. If this flag is provided, the input is expected to be a transition probability matrix, rather than a 2D energy grid. In this case, the first two columns in the CSV file are row-order indices into the transition probability matrix, and the third (final) column is the transition probability of that cell.\n");
  fprintf(stderr, "-V/v\tverbose, the default is disabled. If this flag is provided, light debug data will be printed. To enable heavy debugging, use the flags in mfpt_constants.h\n");
  fprintf(stderr, "-X/x\tsingle basepair moves, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted. This option makes the assumption that the input is *not* a transition probability matrix already, and the input energy grid already satisfies the triangle inequality / parity condition.\n");
  fprintf(stderr, "-Z/z\tend state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the second column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the end state.\n");
  fprintf(stderr, "\nProgram returns -1 (resp. -2) if the start state (resp. end state) probability is 0. -3 is returned if the distance between the two input structures could not be inferred from the input data (usually also means that one of the states has a 0-probability). Otherwise returns the MFPT as predicted by matrix inversion.\n");
  abort();
}
