#ifndef CONFIG_H
#define CONFIG_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define MAX_FILENAME 256
#define RESULTS_BUFFER_SIZE 100

// Weight function types
typedef enum {
    WEIGHT_POWER_LAW,    // omega(d) = d^alpha
    WEIGHT_EXPONENTIAL,  // omega(d) = exp(-d/lambda)
    WEIGHT_BERNOULLI    // omega(d) = theta(d-lambda) + Omega*theta(lambda-d)
} WeightFunctionType;

typedef struct {
    int L;                     // Lattice size
    double sigma;              // Position perturbation std dev
    double alpha;              // Distance exponent (for power law)
    double lambda;             // Threshold distance (for exponential and Bernoulli)
    double Omega;              // Frequency for d < lambda (Bernoulli only)
    WeightFunctionType weight_type;  // Which weight function to use
    double T_INITIAL;          // Initial time
    double T_FINAL;            // Final time
    double dt;                 // Time step
    int M1;                    // Number of simulations per disorder realization
    int M2;                    // Number of disorder realizations
    bool use_pbc;              // true for PBC, false for OBC
    unsigned long seed;        // RNG seed for reproducibility
} SimulationConfig;

typedef struct {
    double x, y;
} Position;

typedef struct {
    double weight;
    int active;
} Edge;

typedef struct {
    double active_fraction;
    double largest_cluster_fraction;
    double avg_cluster_size;  // Correct average: sum(n_s * s^2) / sum(n_s * s), excluding largest
} PercolationStats;

typedef struct {
    double t;
    double avg_prob;
    PercolationStats correlated;    // Changed from "disorder" for clarity
    PercolationStats uncorrelated;  // Shuffled weights (was "uniform")
} TimeStepResults;

typedef struct {
    TimeStepResults *buffer;
    int capacity;
    int count;
} ResultsBuffer;

// Helper function to compute weight based on type
static inline double compute_weight(double distance, SimulationConfig *config) {
    switch (config->weight_type) {
        case WEIGHT_POWER_LAW:
            return pow(distance, config->alpha);
        case WEIGHT_EXPONENTIAL:
            return config->Omega * exp(-distance / config->lambda);
        case WEIGHT_BERNOULLI:
            // theta(d - lambda) + Omega * theta(lambda - d)
            // If d >= lambda: weight = 1
            // If d < lambda: weight = Omega
            return (distance >= config->lambda) ? 1.0 : config->Omega;
        default:
            return pow(distance, config->alpha);
    }
}

#endif