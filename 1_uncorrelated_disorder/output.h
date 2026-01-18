#ifndef OUTPUT_H
#define OUTPUT_H

#include "config.h"

char* generate_filename(const SimulationConfig* config);
void write_results(const char* filename, const SimulationConfig* config, 
                   TimeStepStats* stats, int num_time_steps);

#endif // OUTPUT_H
