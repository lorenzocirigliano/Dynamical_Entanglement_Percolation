#include "percolation.h"

// Calculate proper average cluster size for susceptibility
// S = sum(n_s * s^2) / sum(n_s * s), excluding the largest cluster
// This is the susceptibility-related average cluster size
void get_cluster_statistics(UnionFind* uf, int* largest_size, double* avg_cluster_size) {
    int size = uf->L * uf->L;
    
    *largest_size = 0;
    int largest_root = -1;
    
    // First pass: find the largest cluster
    for (int i = 0; i < size; i++) {
        if (uf->parent[i] == i && uf->size[i] > *largest_size) {
            *largest_size = uf->size[i];
            largest_root = i;
        }
    }
    
    // Second pass: calculate susceptibility-related average
    // S = sum(n_s * s^2) / sum(n_s * s) for all clusters except the largest
    double sum_ns_s2 = 0.0;  // sum of n_s * s^2
    double sum_ns_s = 0.0;   // sum of n_s * s
    
    for (int i = 0; i < size; i++) {
        if (uf->parent[i] == i && i != largest_root) {
            int cluster_size = uf->size[i];
            sum_ns_s2 += (double)cluster_size * cluster_size;  // n_s * s^2 (n_s=1 for each cluster)
            sum_ns_s += (double)cluster_size;                  // n_s * s
        }
    }
    
    // The proper average for susceptibility
    *avg_cluster_size = (sum_ns_s > 0) ? sum_ns_s2 / sum_ns_s : 0.0;
}

// Optimized cluster finding
static inline void find_clusters(Lattice* lattice, UnionFind* uf) {
    int L = lattice->L;
    reset_union_find(uf);
    
    // Process horizontal edges - cache-friendly access pattern
    for (int i = 0; i < L; i++) {
        int row_base = i * L;
        for (int j = 0; j < L; j++) {
            int idx = row_base + j;
            if (lattice->h_edge_valid[idx] && lattice->h_edges[idx].active) {
                int j_next = lattice->use_pbc ? ((j + 1) % L) : (j + 1);
                union_sets(uf, idx, row_base + j_next);
            }
        }
    }
    
    // Process vertical edges
    for (int i = 0; i < L; i++) {
        int row_base = i * L;
        for (int j = 0; j < L; j++) {
            int idx = row_base + j;
            if (lattice->v_edge_valid[idx] && lattice->v_edges[idx].active) {
                int i_next = lattice->use_pbc ? ((i + 1) % L) : (i + 1);
                int next_row_base = i_next * L;
                union_sets(uf, idx, next_row_base + j);
            }
        }
    }
}

PercolationStats simulate_percolation_with_probs(Lattice* lattice, 
                                                  UnionFind* uf, 
                                                  RNG* rng,
                                                  bool use_shuffled) {
    int L = lattice->L;
    int size = L * L;
    int active_edges = 0;
    int total_valid_edges = lattice->n_total_edges;
    
    double *h_probs = use_shuffled ? lattice->shuffled_h_probs : lattice->h_probs;
    double *v_probs = use_shuffled ? lattice->shuffled_v_probs : lattice->v_probs;
    
    // Activate edges based on pre-computed probabilities
    // Process in single loop for better cache usage
    for (int i = 0; i < size; i++) {
        if (lattice->h_edge_valid[i]) {
            lattice->h_edges[i].active = (uniform_random(rng) < h_probs[i]) ? 1 : 0;
            active_edges += lattice->h_edges[i].active;
        } else {
            lattice->h_edges[i].active = 0;
        }
        
        if (lattice->v_edge_valid[i]) {
            lattice->v_edges[i].active = (uniform_random(rng) < v_probs[i]) ? 1 : 0;
            active_edges += lattice->v_edges[i].active;
        } else {
            lattice->v_edges[i].active = 0;
        }
    }
    
    // Find clusters
    find_clusters(lattice, uf);
    
    // Get statistics
    int largest_size;
    double avg_cluster_size;
    get_cluster_statistics(uf, &largest_size, &avg_cluster_size);
    
    PercolationStats stats;
    stats.active_fraction = (total_valid_edges > 0) 
        ? (double)active_edges / total_valid_edges 
        : 0.0;
    stats.largest_cluster_fraction = (double)largest_size / size;
    stats.avg_cluster_size = avg_cluster_size;
    
    return stats;
}
