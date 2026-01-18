#ifndef OUTPUT_H
#define OUTPUT_H

#include "config.h"

char* generate_filename(const SimulationConfig* config);
void write_header(FILE* fp);

#endif  // OUTPUT_H