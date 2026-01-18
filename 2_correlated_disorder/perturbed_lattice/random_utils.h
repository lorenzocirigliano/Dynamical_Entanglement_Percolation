#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include "config.h"

// Better random number generation
typedef struct {
    unsigned long state[4];
} RNG;

void init_rng(RNG* rng, unsigned long seed);
double uniform_random(RNG* rng);
double gaussian_random(RNG* rng, double mu, double sigma);

#endif
