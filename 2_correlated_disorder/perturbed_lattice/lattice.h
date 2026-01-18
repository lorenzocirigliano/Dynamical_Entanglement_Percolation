#ifndef LATTICE_H
#define LATTICE_H

#include "config.h"
#include "random_utils.h"

typedef struct {
    Position* positions;       // Flat array: L*L
    Edge* h_edges;            // Flat array: L*L  
    Edge* v_edges;            // Flat array: L*L
    double* h_probs;          // Pre-computed probabilities (correlated)
    double* v_probs;          // Pre-computed probabilities (correlated)
    
    // For uncorrelated (shuffled) version
    double* shuffled_weights;     // Temporary array for shuffling
    double* shuffled_h_probs;    // Pre-computed probabilities (uncorrelated)
    double* shuffled_v_probs;    // Pre-computed probabilities (uncorrelated)
    
    // Edge validity for OBC
    bool* h_edge_valid;       // Which horizontal edges exist
    bool* v_edge_valid;       // Which vertical edges exist
    
    int L;
    bool use_pbc;
    int n_total_edges;        // Total number of valid edges
} Lattice;

// Index mapping for flat arrays
#define GET_INDEX(i, j, L) ((i) * (L) + (j))

Lattice* create_lattice(int L, bool use_pbc);
void destroy_lattice(Lattice* lattice);
void generate_disorder_realization(Lattice* lattice, double sigma, RNG* rng);
void compute_edge_weights(Lattice* lattice, SimulationConfig* config);
void shuffle_weights(Lattice* lattice, RNG* rng);
void precompute_probabilities(Lattice* lattice, double t, bool use_shuffled);
double compute_average_probability(Lattice* lattice, bool use_shuffled);
double compute_distance(Position p1, Position p2, int L, bool use_pbc);

#endif