#include "percolation.h"
#include "utils.h"
#include <stdlib.h>

void activate_edges(Lattice *lattice, double p1, double p2) {
  for (int i = 0; i < lattice->num_edges; i++) {
    double prob = (lattice->edge_types[i] == 1) ? p1 : p2;
    lattice->edge_active[i] = (random_uniform() < prob) ? 1 : 0;
  }
}

int find_root_local(int *parent, int x) {
  if (parent[x] != x) {
    parent[x] = find_root_local(parent, parent[x]);
  }
  return parent[x];
}

void union_sets_local(int *parent, int *rank, int x, int y) {
  int root_x = find_root_local(parent, x);
  int root_y = find_root_local(parent, y);

  if (root_x != root_y) {
    if (rank[root_x] < rank[root_y]) {
      parent[root_x] = root_y;
    } else if (rank[root_x] > rank[root_y]) {
      parent[root_y] = root_x;
    } else {
      parent[root_y] = root_x;
      rank[root_x]++;
    }
  }
}

ClusterStats calculate_clusters(Lattice *lattice) {
  int L = lattice->L;
  int N = L * L;

  // Union-find arrays
  int *parent = malloc(N * sizeof(int));
  int *rank = malloc(N * sizeof(int));

  // Initialize union-find
  for (int i = 0; i < N; i++) {
    parent[i] = i;
    rank[i] = 0;
  }

  // Union nodes connected by active edges
  for (int edge_idx = 0; edge_idx < lattice->num_edges; edge_idx++) {
    if (lattice->edge_active[edge_idx]) {
      int node1, node2;
      get_edge_endpoints_pbc(lattice, edge_idx, &node1, &node2);
      union_sets_local(parent, rank, node1, node2);
    }
  }

  // Count cluster sizes
  int *cluster_sizes = calloc(N, sizeof(int));
  for (int i = 0; i < N; i++) {
    int root = find_root_local(parent, i);
    cluster_sizes[root]++;
  }

  // Find giant component size
  int giant_size = 0;
  for (int i = 0; i < N; i++) {
    if (cluster_sizes[i] > giant_size) {
      giant_size = cluster_sizes[i];
    }
  }

  // Calculate small cluster statistics
  // In percolation theory, the mean cluster size (susceptibility) is:
  // χ = Σ s² n(s) / Σ s n(s)  (excluding giant component)
  // where s is cluster size and n(s) is number of clusters of size s

  long long sum_s_squared_ns = 0; // Σ s² n(s)
  long long sum_s_ns = 0;         // Σ s n(s)
  int num_small_clusters = 0;

  for (int i = 0; i < N; i++) {
    int s = cluster_sizes[i];
    if (s > 0 && s < giant_size) {
      // This is a small cluster (not the giant component)
      sum_s_squared_ns += (long long)s * s; // s² × 1 (one cluster of size s)
      sum_s_ns += s;                        // s × 1
      num_small_clusters++;
    }
  }

  ClusterStats stats;
  stats.giant_size =
      (double)giant_size / N; // Fraction of lattice in giant component

  // Mean cluster size (susceptibility): χ = ⟨s²⟩ / ⟨s⟩
  stats.avg_small_size =
      (sum_s_ns > 0) ? (double)sum_s_squared_ns / sum_s_ns : 0.0;

  stats.num_small_clusters = num_small_clusters;

  free(parent);
  free(rank);
  free(cluster_sizes);
  return stats;
}