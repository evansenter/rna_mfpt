#ifndef MFPT_PARAMS_H
#define MFPT_PARAMS_H

#include "data_structures.h"

MFPT_PARAMS init_mfpt_params();
MFPT_PARAMS parse_mfpt_args(int, char* []);
int mfpt_error_handling(MFPT_PARAMS);
void debug_mfpt_parameters(MFPT_PARAMS);
void mfpt_usage();

#endif
