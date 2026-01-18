#ifndef PERCOLATION_H
#define PERCOLATION_H

#include "config.h"
#include "lattice.h"
#include "union_find.h"
#include "random_utils.h"

// Main simulation function that works with both correlated and uncorrelated weights
PercolationStats simulate_percolation_with_probs(Lattice* lattice, 
                                                 UnionFind* uf, 
                                                 RNG* rng,
                                                 bool use_shuffled);

// Cluster analysis
void get_cluster_statistics(UnionFind* uf, int* largest_size, double* avg_small_size);

#endif