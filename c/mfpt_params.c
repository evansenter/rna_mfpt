#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include "params.h"

MFPT_PARAMS init_mfpt_params() {
  MFPT_PARAMS parameters = {
    .start_state             = -1,
    .end_state               = -1,
    .max_dist                = 0,
    .bp_dist                 = 0,
    .epsilon                 = 0.,
    .energy_based            = 0,
    .transition_matrix_input = 0,
    .pseudoinverse           = 0,
    .fully_connected         = 0,
    .single_bp_moves_only    = 0,
    .hastings                = 0,
    .radial_probability      = 0,
    .rate_matrix             = 0,
    .verbose                 = 0
  };
  return parameters;
}

MFPT_PARAMS parse_mfpt_args(int argc, char** argv) {
  int c;
  MFPT_PARAMS parameters;
  parameters = init_mfpt_params();
  
  while ((c = getopt(argc, argv, "EeTtPpXxHhRrFfQqVvA:a:Z:z:N:n:D:d:O:o:")) != -1) {
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
        
      case 'R':
      case 'r':
        parameters.rate_matrix = 1;
        break;
        
      case 'F':
      case 'f':
        parameters.fully_connected = 1;
        break;
        
      case 'Q':
      case 'q':
        parameters.radial_probability = 1;
        break;
        
      case 'V':
      case 'v':
        parameters.verbose = 1;
        break;
        
      case 'A':
      case 'a':
        if (!sscanf(optarg, "%d", &parameters.start_state)) {
          mfpt_usage();
        } else if (parameters.start_state < 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'Z':
      case 'z':
        if (!sscanf(optarg, "%d", &parameters.end_state)) {
          mfpt_usage();
        } else if (parameters.end_state < 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'N':
      case 'n':
        if (!sscanf(optarg, "%d", &parameters.max_dist)) {
          mfpt_usage();
        } else if (parameters.max_dist <= 0) {
          mfpt_usage();
        }
        break;
        
      case 'D':
      case 'd':
        if (!sscanf(optarg, "%d", &parameters.bp_dist)) {
          mfpt_usage();
        } else if (parameters.bp_dist <= 0) {
          mfpt_usage();
        }
        
        break;
        
      case 'O':
      case 'o':
        if (!sscanf(optarg, "%lf", &parameters.epsilon)) {
          mfpt_usage();
        } else if (parameters.epsilon <= 1e-16) {
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

int mfpt_error_handling(MFPT_PARAMS parameters) {
  int error = 0;
  
  // Type of run check.
  if (parameters.transition_matrix_input + parameters.single_bp_moves_only + parameters.fully_connected != 1) {
    fprintf(stderr, "Error: Exactly one of -t, -x or -f must be provided!\n");
    error++;
  }
  
  // Transition matrix input requirements.
  if (parameters.transition_matrix_input && !(parameters.start_state >= 0 && parameters.end_state >= 0)) {
    fprintf(stderr, "Error: If the -t flag is provided, -a and -z must be explicitly set!\n");
    error++;
  }
  
  // Transition matrix input restrictions.
  if (parameters.transition_matrix_input && (parameters.hastings || parameters.energy_based)) {
    fprintf(stderr, "Error: If the -t flag is provided, -h and -e are not permitted!\n");
    error++;
  }
  
  // Fully connected graph restrictions.
  if (parameters.fully_connected && (parameters.hastings || parameters.energy_based)) {
    fprintf(stderr, "Error: If the -f flag is provided, -h and -e are not permitted (or -a and -z without -d)!\n");
    error++;
  }
  
  // Extending k/l/p restrictions.
  if (parameters.max_dist && parameters.energy_based) {
    fprintf(stderr, "Error: If the -n flag is provided, -e is not permitted!\n");
    error++;
  }
  
  // If we're extending k/l/p we need to know how the user wants the zero-positions filled.
  if (parameters.max_dist && !(parameters.epsilon || parameters.radial_probability)) {
    fprintf(stderr, "Error: If using the full grid (bounded by -n), either -o or -q needs to be specified for populating 0-probability positions!\n");
    error++;
  }
  
  // Can't populate the zero entries in k/l/p matrix with both of these techniques.
  if (parameters.epsilon > 0 && parameters.radial_probability) {
    fprintf(stderr, "Error: -o and -q are mutually exclusive!\n");
    error++;
  }
  
  // MFPT will be 0.
  if (parameters.start_state == parameters.end_state && parameters.start_state >= 0) {
    fprintf(stderr, "Error: If the -a and -z flags are identical the MFPT is 0!\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_mfpt_parameters(MFPT_PARAMS parameters) {
  char* buffer = calloc(128, sizeof(char));
  
  printf("(e) energy_based\t\t%s\n",          parameters.energy_based            ? "Yes" : "No");
  printf("(t) transition_matrix_input\t%s\n", parameters.transition_matrix_input ? "Yes" : "No");
  printf("(p) pseudoinverse\t\t%s\n",         parameters.pseudoinverse           ? "Yes" : "No");
  printf("(x) single_bp_moves_only\t%s\n",    parameters.single_bp_moves_only    ? "Yes" : "No");
  printf("(f) fully_connected\t\t%s\n",       parameters.fully_connected         ? "Yes" : "No");
  printf("(h) hastings\t\t\t%s\n",            parameters.hastings                ? "Yes" : "No");
  printf("(q) radial_probability\t\t%s\n",    parameters.radial_probability      ? "Yes" : "No");
  printf("(r) rate_matrix\t\t\t%s\n",         parameters.rate_matrix             ? "Yes" : "No");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.max_dist);
  printf("(n) max_dist\t\t\t%s\n", parameters.max_dist ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.bp_dist);
  printf("(d) bp_dist\t\t\t%s\n", parameters.bp_dist ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.start_state);
  printf("(a) start_state\t\t\t%s\n", parameters.start_state >= 0 ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%d", parameters.end_state);
  printf("(z) end_state\t\t\t%s\n", parameters.end_state >= 0 ? buffer : "N/A");
  
  memset(buffer, ' ', 128 * sizeof(char));
  sprintf(buffer, "%.2e", parameters.epsilon);
  printf("(o) epsilon\t\t\t%s\n", parameters.epsilon >= 0 ? buffer : "N/A");
}

void mfpt_usage() {
  fprintf(stderr, "RNAmfpt [options] input_csv\n\n");
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");
  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-A/a\tstart state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the first column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the start state.\n");
  fprintf(stderr, "-D/d\tstart/end distance, the default is disabled. When provided, indicates the base pair distance between the starting / ending structures. This flag is used in conjunction with the -n flag, and is needed in cases when the base pair distance between the two structures can't be inferred from the input grid.\n");
  fprintf(stderr, "-E/e\tenergy-based transitions, the default is disabled. If this flag is provided, the transition from state a to b will be calculated as (min(1, exp(-(E_b - E_a) / RT) / n) rather than (min(1, p_b / p_a) / n).\n");
  fprintf(stderr, "-F/f\tfully connected, if the graph is fully connected we permit transitions between arbitrary positions in the 2D grid.\n");
  fprintf(stderr, "-H/h\tHastings adjustment, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted (in the all-to-all transition case, N(X) / N(Y) == 1). Calculating N(X) and N(Y) will respect grid boundaries and the triangle equality, and the basepair distance between the two structures for kinetics is inferred from the energy grid.\n");
  fprintf(stderr, "-N/n\tsequence length, the default is disabled. This flag represents the sequence length of the sequence on which kinetics is being performed. It is used to ensure that the graph is fully connected.\n");
  fprintf(stderr, "-O/o\tepsilon. If the graph is going to be populated with all possible moves (via the -n flag), this will inflate all 0-probability positions.\n");
  fprintf(stderr, "-P/p\tpseudoinverse, the default is disabled. If this flag is provided, the Moore-Penrose pseudoinverse is computed for the transition probability matrix, rather than the true inverse.\n");
  fprintf(stderr, "-Q/q\tradial probability. If the graph is going to be populated with all possible moves (via the -n flag), this will adjust all transition probabilities from 0-probability positions.\n");
  fprintf(stderr, "-R/r\trate matrix, the default is disabled. If this flag is provided, the transition rate matrix is computed rather than the transition probability matrix.\n");
  fprintf(stderr, "-T/t\ttransition matrix input, the default is disabled. If this flag is provided, the input is expected to be a transition probability matrix, rather than a 2D energy grid. In this case, the first two columns in the CSV file are row-order indices into the transition probability matrix, and the third (final) column is the transition probability of that cell.\n");
  fprintf(stderr, "-V/v\tverbose, the default is disabled. If this flag is provided, light debug data will be printed. To enable heavy debugging, use the flags in mfpt_constants.h\n");
  fprintf(stderr, "-X/x\tsingle basepair moves, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted. This option makes the assumption that the input is *not* a transition probability matrix already, and the input energy grid already satisfies the triangle inequality / parity condition.\n");
  fprintf(stderr, "-Z/z\tend state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the second column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the end state.\n");
  fprintf(stderr, "\nProgram returns -1 (resp. -2) if the start state (resp. end state) probability is 0. -3 is returned if the distance between the two input structures could not be inferred from the input data (usually also means that one of the states has a 0-probability). Otherwise returns the MFPT as predicted by matrix inversion.\n");
  abort();
}
