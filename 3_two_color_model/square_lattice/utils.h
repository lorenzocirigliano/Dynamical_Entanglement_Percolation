#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct {
    double mean;
    double std_error;
} Statistics;

// Function declarations
void init_random(void);
void init_random_seed(unsigned int seed);
double random_uniform(void);
Statistics calculate_statistics(double *data, int n);

#endif
