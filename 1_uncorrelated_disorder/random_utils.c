#include "random_utils.h"

void init_random(unsigned int seed) {
    srand(seed);
}

double uniform_random(void) {
    return (double)rand() / RAND_MAX;
}

// Box-Muller transform for Gaussian random numbers
double gaussian_random(double mu, double sigma) {
    static int has_spare = 0;
    static double spare;
    
    if (has_spare) {
        has_spare = 0;
        return spare * sigma + mu;
    }
    
    has_spare = 1;
    double u = uniform_random() * 2.0 - 1.0;
    double v = uniform_random() * 2.0 - 1.0;
    double s = u * u + v * v;
    
    while (s >= 1.0 || s == 0.0) {
        u = uniform_random() * 2.0 - 1.0;
        v = uniform_random() * 2.0 - 1.0;
        s = u * u + v * v;
    }
    
    double mult = sqrt(-2.0 * log(s) / s);
    spare = v * mult;
    return mu + sigma * u * mult;
}

// Bernoulli random: returns 1.0 with probability p, omega_high with probability (1-p)
double bernoulli_random(double p, double omega_high) {
    double r = uniform_random();
    if (r < p) {
        return 1.0;
    } else {
        return omega_high;
    }
}
