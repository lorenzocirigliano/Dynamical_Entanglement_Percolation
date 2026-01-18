#include "lattice.h"

double compute_distance(Position p1, Position p2, int L, bool use_pbc) {
    double dx = fabs(p1.x - p2.x);
    double dy = fabs(p1.y - p2.y);
    
    if (use_pbc) {
        if (dx > L / 2.0) dx = L - dx;
        if (dy > L / 2.0) dy = L - dy;
    }
    
    return sqrt(dx * dx + dy * dy);
}

Lattice *create_lattice(int L, bool use_pbc) {
    Lattice *lattice = malloc(sizeof(Lattice));
    if (!lattice) return NULL;
    
    lattice->L = L;
    lattice->use_pbc = use_pbc;
    int size = L * L;
    
    // Determine number of edges based on boundary conditions
    int n_h_edges = use_pbc ? size : L * (L - 1);
    int n_v_edges = use_pbc ? size : (L - 1) * L;
    lattice->n_total_edges = n_h_edges + n_v_edges;
    
    // Allocate flat arrays
    lattice->positions = malloc(size * sizeof(Position));
    lattice->h_edges = malloc(size * sizeof(Edge));
    lattice->v_edges = malloc(size * sizeof(Edge));
    lattice->h_probs = malloc(size * sizeof(double));
    lattice->v_probs = malloc(size * sizeof(double));
    
    // For shuffled weights version
    lattice->shuffled_weights = malloc(lattice->n_total_edges * sizeof(double));
    lattice->shuffled_h_probs = malloc(size * sizeof(double));
    lattice->shuffled_v_probs = malloc(size * sizeof(double));
    
    // Edge validity arrays for OBC
    lattice->h_edge_valid = malloc(size * sizeof(bool));
    lattice->v_edge_valid = malloc(size * sizeof(bool));
    
    if (!lattice->positions || !lattice->h_edges || !lattice->v_edges ||
        !lattice->h_probs || !lattice->v_probs || !lattice->shuffled_weights ||
        !lattice->shuffled_h_probs || !lattice->shuffled_v_probs ||
        !lattice->h_edge_valid || !lattice->v_edge_valid) {
        destroy_lattice(lattice);
        return NULL;
    }
    
    // Initialize edge validity
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = GET_INDEX(i, j, L);
            // Horizontal edges: invalid at right boundary for OBC
            lattice->h_edge_valid[idx] = use_pbc || (j < L - 1);
            // Vertical edges: invalid at bottom boundary for OBC
            lattice->v_edge_valid[idx] = use_pbc || (i < L - 1);
        }
    }
    
    return lattice;
}

void destroy_lattice(Lattice *lattice) {
    if (!lattice) return;
    
    free(lattice->positions);
    free(lattice->h_edges);
    free(lattice->v_edges);
    free(lattice->h_probs);
    free(lattice->v_probs);
    free(lattice->shuffled_weights);
    free(lattice->shuffled_h_probs);
    free(lattice->shuffled_v_probs);
    free(lattice->h_edge_valid);
    free(lattice->v_edge_valid);
    free(lattice);
}

void generate_disorder_realization(Lattice *lattice, double sigma, RNG *rng) {
    int L = lattice->L;
    
    // Generate perturbed positions for this disorder realization
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = GET_INDEX(i, j, L);
            lattice->positions[idx].x = i + gaussian_random(rng, 0.0, sigma);
            lattice->positions[idx].y = j + gaussian_random(rng, 0.0, sigma);
        }
    }
}

void compute_edge_weights(Lattice *lattice, SimulationConfig *config) {
    int L = lattice->L;
    
    // Horizontal edges
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = GET_INDEX(i, j, L);
            
            if (lattice->h_edge_valid[idx]) {
                int j_next = lattice->use_pbc ? ((j + 1) % L) : (j + 1);
                int idx_next = GET_INDEX(i, j_next, L);
                
                double distance = compute_distance(lattice->positions[idx],
                                                  lattice->positions[idx_next], 
                                                  L, lattice->use_pbc);
                lattice->h_edges[idx].weight = compute_weight(distance, config);
            } else {
                lattice->h_edges[idx].weight = 0.0;
            }
        }
    }
    
    // Vertical edges
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = GET_INDEX(i, j, L);
            
            if (lattice->v_edge_valid[idx]) {
                int i_next = lattice->use_pbc ? ((i + 1) % L) : (i + 1);
                int idx_next = GET_INDEX(i_next, j, L);
                
                double distance = compute_distance(lattice->positions[idx],
                                                  lattice->positions[idx_next], 
                                                  L, lattice->use_pbc);
                lattice->v_edges[idx].weight = compute_weight(distance, config);
            } else {
                lattice->v_edges[idx].weight = 0.0;
            }
        }
    }
}

void shuffle_weights(Lattice *lattice, RNG *rng) {
    int L = lattice->L;
    int edge_count = 0;
    
    // Collect all valid edge weights
    for (int idx = 0; idx < L * L; idx++) {
        if (lattice->h_edge_valid[idx]) {
            lattice->shuffled_weights[edge_count++] = lattice->h_edges[idx].weight;
        }
    }
    for (int idx = 0; idx < L * L; idx++) {
        if (lattice->v_edge_valid[idx]) {
            lattice->shuffled_weights[edge_count++] = lattice->v_edges[idx].weight;
        }
    }
    
    // Fisher-Yates shuffle
    for (int i = edge_count - 1; i > 0; i--) {
        int j = (int)(uniform_random(rng) * (i + 1));
        if (j > i) j = i;  // Safety check
        double temp = lattice->shuffled_weights[i];
        lattice->shuffled_weights[i] = lattice->shuffled_weights[j];
        lattice->shuffled_weights[j] = temp;
    }
    
    // Redistribute shuffled weights back to edges
    edge_count = 0;
    for (int idx = 0; idx < L * L; idx++) {
        if (lattice->h_edge_valid[idx]) {
            lattice->h_edges[idx].weight = lattice->shuffled_weights[edge_count++];
        }
    }
    for (int idx = 0; idx < L * L; idx++) {
        if (lattice->v_edge_valid[idx]) {
            lattice->v_edges[idx].weight = lattice->shuffled_weights[edge_count++];
        }
    }
}

void precompute_probabilities(Lattice *lattice, double t, bool use_shuffled) {
    int size = lattice->L * lattice->L;
    double *h_probs = use_shuffled ? lattice->shuffled_h_probs : lattice->h_probs;
    double *v_probs = use_shuffled ? lattice->shuffled_v_probs : lattice->v_probs;
    
    for (int i = 0; i < size; i++) {
        if (lattice->h_edge_valid[i]) {
            double cos_val = cos(lattice->h_edges[i].weight * t);
            h_probs[i] = 1.0 - fabs(cos_val);
        } else {
            h_probs[i] = 0.0;
        }
        
        if (lattice->v_edge_valid[i]) {
            double cos_val = cos(lattice->v_edges[i].weight * t);
            v_probs[i] = 1.0 - fabs(cos_val);
        } else {
            v_probs[i] = 0.0;
        }
    }
}

double compute_average_probability(Lattice *lattice, bool use_shuffled) {
    int size = lattice->L * lattice->L;
    double sum = 0.0;
    int valid_edges = 0;
    
    double *h_probs = use_shuffled ? lattice->shuffled_h_probs : lattice->h_probs;
    double *v_probs = use_shuffled ? lattice->shuffled_v_probs : lattice->v_probs;
    
    for (int i = 0; i < size; i++) {
        if (lattice->h_edge_valid[i]) {
            sum += h_probs[i];
            valid_edges++;
        }
        if (lattice->v_edge_valid[i]) {
            sum += v_probs[i];
            valid_edges++;
        }
    }
    
    return (valid_edges > 0) ? sum / valid_edges : 0.0;
}