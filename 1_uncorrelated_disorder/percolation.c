#include "percolation.h"

// Find connected components using Union-Find
void find_clusters(UnionFind* uf, Lattice* lattice) {
    int L = lattice->L;
    
    // Reset Union-Find
    reset_union_find(uf);
    
    // Process horizontal edges (with periodic boundary conditions)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (lattice->h_edges[i][j].active) {
                int node1 = get_node_index(i, j, L);
                int node2 = get_node_index(i, (j + 1) % L, L); // Periodic in j
                union_sets(uf, node1, node2);
            }
        }
    }
    
    // Process vertical edges (with periodic boundary conditions)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (lattice->v_edges[i][j].active) {
                int node1 = get_node_index(i, j, L);
                int node2 = get_node_index((i + 1) % L, j, L); // Periodic in i
                union_sets(uf, node1, node2);
            }
        }
    }
}

void get_cluster_statistics(UnionFind* uf, int* largest_size, double* susceptibility_avg) {
    int L = uf->L;
    int size = L * L;
    
    int* cluster_sizes = calloc(size, sizeof(int));
    int num_clusters = 0;
    
    // Count cluster sizes
    for (int i = 0; i < size; i++) {
        int root = find(uf, i);
        if (root == i) { // This is a root node
            cluster_sizes[num_clusters++] = uf->size[i];
        }
    }
    
    // Find largest cluster
    *largest_size = 0;
    for (int i = 0; i < num_clusters; i++) {
        if (cluster_sizes[i] > *largest_size) {
            *largest_size = cluster_sizes[i];
        }
    }
    
    // Calculate susceptibility-related quantity: sum(s^2 * n_s) / sum(s * n_s)
    // where n_s is the number of clusters of size s
    // This is equivalent to: sum_{clusters not largest}(s^2) / sum_{clusters not largest}(s)
    long long sum_s_squared = 0; // sum of s^2 for non-largest clusters
    long long sum_s = 0;         // sum of s for non-largest clusters
    
    for (int i = 0; i < num_clusters; i++) {
        if (cluster_sizes[i] != *largest_size) {
            int s = cluster_sizes[i];
            sum_s_squared += (long long)s * s;
            sum_s += s;
        }
    }
    
    *susceptibility_avg = (sum_s > 0) ? (double)sum_s_squared / sum_s : 0.0;
    
    free(cluster_sizes);
}
