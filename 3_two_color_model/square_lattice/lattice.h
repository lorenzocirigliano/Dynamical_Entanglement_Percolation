#ifndef LATTICE_H
#define LATTICE_H

// Disorder mode enumeration
typedef enum {
  MODE_CORRELATED,  // Local constraints enforced
  MODE_ALTERNATING, // Edges alternate along directions
  MODE_UNCORRELATED // Random assignment with probability f
} DisorderMode;

typedef struct {
  int L;                   // Lattice size
  int num_nodes;           // Total number of nodes (L*L)
  int num_edges;           // Total number of edges
  int degree;              // Node degree (4 for square lattice)
  int *edge_types;         // Array of edge types (1 or 2)
  int *edge_active;        // Array indicating if edge is active
  int unconstrained_edges; // Number of edges that couldn't satisfy constraints
} Lattice;

// Core functions
Lattice *create_square_lattice_pbc(int L);
void destroy_lattice(Lattice *lattice);

// Main assignment function
void assign_edge_types(Lattice *lattice, double f, DisorderMode mode);

// Helper functions
int get_horizontal_edge_index(int L, int i, int j);
int get_vertical_edge_index(int L, int i, int j);
void get_edge_endpoints_pbc(const Lattice *lattice, int edge_idx, int *node1,
                            int *node2);
void verify_edge_assignment(const Lattice *lattice, double f,
                            DisorderMode mode);

// Internal assignment functions (called by assign_edge_types)
void assign_edge_types_correlated(Lattice *lattice, double f);
void assign_edge_types_alternating(Lattice *lattice, double f);
void assign_edge_types_uncorrelated(Lattice *lattice, double f);
void assign_edge_types_random_constrained(Lattice *lattice, double f);

#endif