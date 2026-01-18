#ifndef PERCOLATION_H
#define PERCOLATION_H

#include "config.h"
#include "lattice.h"
#include "union_find.h"

void find_clusters(UnionFind* uf, Lattice* lattice);
void get_cluster_statistics(UnionFind* uf, int* largest_size, double* susceptibility_avg);

#endif // PERCOLATION_H
