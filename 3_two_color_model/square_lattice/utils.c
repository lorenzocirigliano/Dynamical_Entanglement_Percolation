#include "utils.h"

void init_random(void) {
    srand(time(NULL));
}

void init_random_seed(unsigned int seed) {
    srand(seed);
}

double random_uniform(void) {
    return (double)rand() / RAND_MAX;
}

Statistics calculate_statistics(double *data, int n) {
    Statistics stats;
    double sum = 0.0;
    double sum_sq = 0.0;
    
    for (int i = 0; i < n; i++) {
        sum += data[i];
        sum_sq += data[i] * data[i];
    }
    
    stats.mean = sum / n;
    double variance = (sum_sq / n) - (stats.mean * stats.mean);
    stats.std_error = sqrt(variance / n);
    
    return stats;
}
