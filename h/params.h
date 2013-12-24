#ifndef MFPT_PARAMS_H
#define MFPT_PARAMS_H

#include "data_structures.h"

MFPT_PARAMETERS init_mfpt_params();
MFPT_PARAMETERS parse_mfpt_args(int, char*[]);
int mfpt_error_handling(MFPT_PARAMETERS);
void debug_mfpt_parameters(MFPT_PARAMETERS);
void mfpt_usage();

#endif
