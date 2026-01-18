#include "lattice.h"
#include "random_utils.h"

Lattice* create_lattice(int L) {
    Lattice* lattice = malloc(sizeof(Lattice));
    if (!lattice) return NULL;
    
    lattice->L = L;
    
    // Allocate edge arrays
    lattice->h_edges = malloc(L * sizeof(Edge*));
    lattice->v_edges = malloc(L * sizeof(Edge*));
    
    if (!lattice->h_edges || !lattice->v_edges) {
        destroy_lattice(lattice);
        return NULL;
    }
    
    for (int i = 0; i < L; i++) {
        lattice->h_edges[i] = malloc(L * sizeof(Edge));
        lattice->v_edges[i] = malloc(L * sizeof(Edge));
        
        if (!lattice->h_edges[i] || !lattice->v_edges[i]) {
            destroy_lattice(lattice);
            return NULL;
        }
    }
    
    return lattice;
}

void destroy_lattice(Lattice* lattice) {
    if (!lattice) return;
    
    if (lattice->h_edges) {
        for (int i = 0; i < lattice->L; i++) {
            free(lattice->h_edges[i]);
        }
        free(lattice->h_edges);
    }
    
    if (lattice->v_edges) {
        for (int i = 0; i < lattice->L; i++) {
            free(lattice->v_edges[i]);
        }
        free(lattice->v_edges);
    }
    
    free(lattice);
}

void generate_weights(Lattice* lattice, const SimulationConfig* config) {
    int L = lattice->L;
    
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (config->dist_type == DIST_GAUSSIAN) {
                lattice->h_edges[i][j].weight = gaussian_random(config->param1, config->param2);
                lattice->v_edges[i][j].weight = gaussian_random(config->param1, config->param2);
            } else { // DIST_BERNOULLI
                lattice->h_edges[i][j].weight = bernoulli_random(config->param1, config->param2);
                lattice->v_edges[i][j].weight = bernoulli_random(config->param1, config->param2);
            }
        }
    }
}

double activation_probability(double omega, double t) {
    return 1.0 - fabs(cos(omega * t));
}

int get_node_index(int i, int j, int L) {
    return i * L + j;
}

void get_coordinates(int index, int L, int* i, int* j) {
    *i = index / L;
    *j = index % L;
}
