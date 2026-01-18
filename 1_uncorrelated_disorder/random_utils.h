#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include "config.h"

double gaussian_random(double mu, double sigma);
double bernoulli_random(double p, double omega_high);
void init_random(unsigned int seed);
double uniform_random(void);

#endif // RANDOM_UTILS_H
