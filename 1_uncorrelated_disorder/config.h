#ifndef CONFIG_H
#define CONFIG_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define MAX_FILENAME 256

// Distribution types
typedef enum {
    DIST_GAUSSIAN,
    DIST_BERNOULLI
} DistributionType;

// Simulation configuration
typedef struct {
    int L;                      // Lattice size
    DistributionType dist_type; // Distribution type
    double param1;              // mu (Gaussian) or p (Bernoulli)
    double param2;              // sigma (Gaussian) or Omega (Bernoulli)
    int M;                      // Number of disorder realizations
    int D;                      // Number of percolation processes per time
    double T_INITIAL;           // Initial time
    double T_FINAL;             // Final time
    double dt;                  // Time step
    unsigned int seed;          // Random seed
} SimulationConfig;

// Edge structure
typedef struct {
    double weight;
    int active;
} Edge;

// Statistics for a single time step
typedef struct {
    double t;                   // Time
    double avg_prob;            // Average activation probability
    double mean_active;         // Mean active fraction
    double std_active;          // Std dev active fraction
    double mean_largest;        // Mean largest cluster fraction
    double std_largest;         // Std dev largest cluster fraction
    double mean_suscept;        // Mean susceptibility
    double std_suscept;         // Std dev susceptibility
} TimeStepStats;

#endif // CONFIG_H
