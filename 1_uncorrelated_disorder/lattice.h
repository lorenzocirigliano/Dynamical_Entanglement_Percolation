#ifndef LATTICE_H
#define LATTICE_H

#include "config.h"

typedef struct {
    Edge** h_edges;     // Horizontal edges [L][L]
    Edge** v_edges;     // Vertical edges [L][L]
    int L;
} Lattice;

Lattice* create_lattice(int L);
void destroy_lattice(Lattice* lattice);
void generate_weights(Lattice* lattice, const SimulationConfig* config);
double activation_probability(double omega, double t);
int get_node_index(int i, int j, int L);
void get_coordinates(int index, int L, int* i, int* j);

#endif // LATTICE_H
