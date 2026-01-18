#ifndef IO_H
#define IO_H

#include <stdio.h>

typedef struct {
    int L;
    int M;
    double f;
    double p1_min, p1_max;
    double p2_min, p2_max;
    double delta_p;
    unsigned int seed;
} SimulationParams;

// Function declarations
int parse_command_line(int argc, char *argv[], SimulationParams *params);
FILE* create_output_file(const SimulationParams *params);
void write_header(FILE *fp);
void write_data_point(FILE *fp, double p1, double p2, double p_bar, 
                      double gc_mean, double gc_err, 
                      double small_mean, double small_err);
void print_usage(const char *program_name);

#endif
