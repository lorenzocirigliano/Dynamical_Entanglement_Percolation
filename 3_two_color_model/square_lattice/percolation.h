#ifndef PERCOLATION_H
#define PERCOLATION_H

#include "lattice.h"

typedef struct {
    double giant_size;
    double avg_small_size;
    int num_small_clusters;
} ClusterStats;

// Function declarations
void activate_edges(Lattice *lattice, double p1, double p2);
ClusterStats calculate_clusters(Lattice *lattice);
int find_root_local(int *parent, int x);
void union_sets_local(int *parent, int *rank, int x, int y);

#endif
