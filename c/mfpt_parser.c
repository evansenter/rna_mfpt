#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "params.h"
#include "parser.h"

int count_lines(char* file_path) {
  FILE* file = fopen(file_path, "r");
  int c;
  int line_count = 0;

  if (file == NULL) {
    fprintf(stderr, "File not found.\n");
    fclose(file);
    return 0;
  }

  while ((c = fgetc(file)) != EOF) {
    if (c == '\n') {
      line_count++;
    }
  }

  fclose(file);
  return line_count;
}

void populate_arrays(char* file_path, KLP_MATRIX klp_matrix, MFPT_PARAMS parameters) {
  int i = 0;
  FILE* file = fopen(file_path, "r");
  char* token;
  char line[1024];

  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    klp_matrix.k[i]  = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix.l[i]  = atoi(token);
    token = strtok(NULL, ",");
    klp_matrix.p[i]  = atof(token);

    if (!parameters.energy_based && (klp_matrix.p[i] < 0 || klp_matrix.p[i] > 1)) {
      fprintf(stderr, "Error: line number %d (0-indexed) in the input doesn't satisfy 0 <= probability (%+1.2f) <= 1. Did you forget the -e flag?\n\n", i, klp_matrix.p[i]);
      mfpt_usage();
    }

    i++;
  }

  fclose(file);
}
