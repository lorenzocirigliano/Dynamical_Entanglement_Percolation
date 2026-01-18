#include "lattice.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Lattice *create_square_lattice_pbc(int L) {
  Lattice *lattice = malloc(sizeof(Lattice));
  if (!lattice)
    return NULL;

  lattice->L = L;
  lattice->num_nodes = L * L;
  lattice->degree = 4;
  lattice->num_edges = 2 * L * L;
  lattice->unconstrained_edges = 0;

  lattice->edge_types = malloc(lattice->num_edges * sizeof(int));
  lattice->edge_active = malloc(lattice->num_edges * sizeof(int));

  if (!lattice->edge_types || !lattice->edge_active) {
    destroy_lattice(lattice);
    return NULL;
  }

  // Initialize all edges as inactive and type 2
  for (int i = 0; i < lattice->num_edges; i++) {
    lattice->edge_types[i] = 2;
    lattice->edge_active[i] = 0;
  }

  return lattice;
}

void destroy_lattice(Lattice *lattice) {
  if (!lattice)
    return;
  free(lattice->edge_types);
  free(lattice->edge_active);
  free(lattice);
}

int get_horizontal_edge_index(int L, int i, int j) { return i * L + j; }

int get_vertical_edge_index(int L, int i, int j) { return L * L + i * L + j; }

void get_edge_endpoints_pbc(const Lattice *lattice, int edge_idx, int *node1,
                            int *node2) {
  int L = lattice->L;

  if (edge_idx < L * L) {
    // Horizontal edge
    int i = edge_idx / L;
    int j = edge_idx % L;
    *node1 = i * L + j;
    *node2 = i * L + ((j + 1) % L);
  } else {
    // Vertical edge
    int vert_idx = edge_idx - L * L;
    int i = vert_idx / L;
    int j = vert_idx % L;
    *node1 = i * L + j;
    *node2 = ((i + 1) % L) * L + j;
  }
}

// ============================================================================
// MODE 1: CORRELATED - Exact local constraints
// ============================================================================
void assign_edge_types_correlated(Lattice *lattice, double f) {
  int L = lattice->L;
  int target_type1_per_node = (int)(f * lattice->degree);

  printf("CORRELATED: f=%.2f (%d type-1 edges per node)\n", f,
         target_type1_per_node);

  // Initialize all edges as type 2
  for (int i = 0; i < lattice->num_edges; i++) {
    lattice->edge_types[i] = 2;
  }

  lattice->unconstrained_edges = 0;

  switch (target_type1_per_node) {
  case 0: // f = 0.0
    printf("  Pattern: All edges type 2\n");
    break;

  case 1: // f = 0.25
    printf("  Pattern: Checkerboard with 1 type-1 edge per node\n");
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        if ((i + j) % 2 == 0) {
          int vert_idx = get_vertical_edge_index(L, i, j);
          lattice->edge_types[vert_idx] = 1;
        }
      }
    }
    break;

  case 2: // f = 0.5
    printf("  Pattern: Horizontal type 1, vertical type 2\n");
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        int horiz_idx = get_horizontal_edge_index(L, i, j);
        lattice->edge_types[horiz_idx] = 1;
      }
    }
    break;

  case 3: // f = 0.75
    printf("  Pattern: 3 type-1 edges per node\n");
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        int horiz_idx = get_horizontal_edge_index(L, i, j);
        lattice->edge_types[horiz_idx] = 1;

        if ((i + j) % 2 == 0) {
          int vert_idx = get_vertical_edge_index(L, i, j);
          lattice->edge_types[vert_idx] = 1;
        }
      }
    }
    break;

  case 4: // f = 1.0
    printf("  Pattern: All edges type 1\n");
    for (int i = 0; i < lattice->num_edges; i++) {
      lattice->edge_types[i] = 1;
    }
    break;

  default:
    printf("  Warning: Non-standard f value, using random constrained\n");
    assign_edge_types_random_constrained(lattice, f);
    break;
  }
}

// ============================================================================
// MODE 2: ALTERNATING - Edges alternate along directions
// ============================================================================
void assign_edge_types_alternating(Lattice *lattice, double f) {
  int L = lattice->L;
  int target_type1_per_node = (int)(f * lattice->degree);

  printf("ALTERNATING: f=%.2f (target %d type-1 per node)\n", f,
         target_type1_per_node);

  // Initialize all edges as type 2
  for (int i = 0; i < lattice->num_edges; i++) {
    lattice->edge_types[i] = 2;
  }

  if (target_type1_per_node == 2) {
    // f = 0.5: Perfect alternating pattern with independent starting colors
    printf("  Pattern: Alternation with independent row/column starts\n");

    // Random starting type for each row (horizontal edges)
    int *row_starts = malloc(L * sizeof(int));
    for (int i = 0; i < L; i++) {
      row_starts[i] = (random_uniform() < 0.5) ? 1 : 2;
    }

    // Random starting type for each column (vertical edges)
    int *col_starts = malloc(L * sizeof(int));
    for (int j = 0; j < L; j++) {
      col_starts[j] = (random_uniform() < 0.5) ? 1 : 2;
    }

    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        // Horizontal edges: alternate along row, starting with row_starts[i]
        int h_idx = get_horizontal_edge_index(L, i, j);
        int h_start = row_starts[i];
        lattice->edge_types[h_idx] = (j % 2 == 0) ? h_start : (3 - h_start);

        // Vertical edges: alternate along column, starting with col_starts[j]
        int v_idx = get_vertical_edge_index(L, i, j);
        int v_start = col_starts[j];
        lattice->edge_types[v_idx] = (i % 2 == 0) ? v_start : (3 - v_start);
      }
    }

    free(row_starts);
    free(col_starts);
  } else if (target_type1_per_node == 0) {
    printf("  Pattern: All edges type 2\n");
    // Already initialized as type 2
  } else if (target_type1_per_node == 4) {
    printf("  Pattern: All edges type 1\n");
    for (int i = 0; i < lattice->num_edges; i++) {
      lattice->edge_types[i] = 1;
    }
  } else {
    // For f != 0.5, alternation is imperfect
    printf("  Warning: f=%.2f cannot achieve perfect alternation\n", f);
    printf("  Using best approximation with alternating pattern\n");

    // Default to alternating pattern (will give ~50% type 1)
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        int h_idx = get_horizontal_edge_index(L, i, j);
        int v_idx = get_vertical_edge_index(L, i, j);
        lattice->edge_types[h_idx] = (j % 2 == 0) ? 1 : 2;
        lattice->edge_types[v_idx] = (i % 2 == 0) ? 1 : 2;
      }
    }
  }

  if (L % 2 == 1) {
    printf("  Note: Odd L=%d may cause PBC violations in alternation\n", L);
  }
}

// ============================================================================
// MODE 3: UNCORRELATED - Random assignment with probability f
// ============================================================================
void assign_edge_types_uncorrelated(Lattice *lattice, double f) {
  printf("UNCORRELATED: Each edge independently type-1 with probability %.3f\n",
         f);

  // Each edge independently becomes type 1 with probability f
  for (int i = 0; i < lattice->num_edges; i++) {
    lattice->edge_types[i] = (random_uniform() < f) ? 1 : 2;
  }

  lattice->unconstrained_edges = 0;
}

// ============================================================================
// MAIN ASSIGNMENT FUNCTION
// ============================================================================
void assign_edge_types(Lattice *lattice, double f, DisorderMode mode) {
  printf("\n=== EDGE TYPE ASSIGNMENT ===\n");

  switch (mode) {
  case MODE_CORRELATED:
    assign_edge_types_correlated(lattice, f);
    break;
  case MODE_ALTERNATING:
    assign_edge_types_alternating(lattice, f);
    break;
  case MODE_UNCORRELATED:
    assign_edge_types_uncorrelated(lattice, f);
    break;
  default:
    printf("Error: Unknown disorder mode\n");
    assign_edge_types_correlated(lattice, f);
  }

  verify_edge_assignment(lattice, f, mode);
}

// ============================================================================
// HELPER: Random constrained assignment
// ============================================================================
void assign_edge_types_random_constrained(Lattice *lattice, double f) {
  int L = lattice->L;
  int target_type1_per_node = (int)(f * lattice->degree);

  for (int i = 0; i < lattice->num_edges; i++) {
    lattice->edge_types[i] = 2;
  }

  int *node_type1_count = calloc(L * L, sizeof(int));
  lattice->unconstrained_edges = 0;

  int *edge_order = malloc(lattice->num_edges * sizeof(int));
  for (int i = 0; i < lattice->num_edges; i++) {
    edge_order[i] = i;
  }

  // Fisher-Yates shuffle
  for (int i = lattice->num_edges - 1; i > 0; i--) {
    int j = (int)(random_uniform() * (i + 1));
    int temp = edge_order[i];
    edge_order[i] = edge_order[j];
    edge_order[j] = temp;
  }

  for (int idx = 0; idx < lattice->num_edges; idx++) {
    int edge_idx = edge_order[idx];
    int node1, node2;
    get_edge_endpoints_pbc(lattice, edge_idx, &node1, &node2);

    if (node_type1_count[node1] < target_type1_per_node &&
        node_type1_count[node2] < target_type1_per_node) {

      lattice->edge_types[edge_idx] = 1;
      node_type1_count[node1]++;
      node_type1_count[node2]++;
    }
  }

  for (int node = 0; node < L * L; node++) {
    if (node_type1_count[node] != target_type1_per_node) {
      lattice->unconstrained_edges++;
    }
  }

  free(node_type1_count);
  free(edge_order);
}

// ============================================================================
// VERIFICATION
// ============================================================================
void verify_edge_assignment(const Lattice *lattice, double f,
                            DisorderMode mode) {
  int L = lattice->L;
  int target_type1_per_node = (int)(f * lattice->degree);
  int *node_type1_count = calloc(L * L, sizeof(int));

  int total_type1_edges = 0;
  for (int edge_idx = 0; edge_idx < lattice->num_edges; edge_idx++) {
    if (lattice->edge_types[edge_idx] == 1) {
      total_type1_edges++;
      int node1, node2;
      get_edge_endpoints_pbc(lattice, edge_idx, &node1, &node2);
      node_type1_count[node1]++;
      node_type1_count[node2]++;
    }
  }

  // Calculate statistics
  double mean = 0, variance = 0;
  int satisfied_nodes = 0;
  int min_degree = 4, max_degree = 0;

  for (int i = 0; i < L * L; i++) {
    mean += node_type1_count[i];
    if (node_type1_count[i] == target_type1_per_node)
      satisfied_nodes++;
    if (node_type1_count[i] < min_degree)
      min_degree = node_type1_count[i];
    if (node_type1_count[i] > max_degree)
      max_degree = node_type1_count[i];
  }
  mean /= (L * L);

  for (int i = 0; i < L * L; i++) {
    double diff = node_type1_count[i] - mean;
    variance += diff * diff;
  }
  variance /= (L * L);

  printf("\nVerification Results:\n");
  printf("  Total type-1 edges: %d/%d (%.3f, target %.3f)\n", total_type1_edges,
         lattice->num_edges, (double)total_type1_edges / lattice->num_edges, f);
  printf("  Type-1 edges per node: mean=%.3f, std=%.3f\n", mean,
         sqrt(variance));

  if (mode == MODE_CORRELATED) {
    printf("  Constraint satisfaction: %d/%d nodes (%.1f%%)\n", satisfied_nodes,
           L * L, 100.0 * satisfied_nodes / (L * L));
  } else if (mode == MODE_ALTERNATING) {
    // Check alternation violations
    int h_violations = 0, v_violations = 0;

    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
        int edge1 = get_horizontal_edge_index(L, i, j);
        int edge2 = get_horizontal_edge_index(L, i, (j + 1) % L);
        if (lattice->edge_types[edge1] == lattice->edge_types[edge2]) {
          h_violations++;
        }
      }
    }

    for (int j = 0; j < L; j++) {
      for (int i = 0; i < L; i++) {
        int edge1 = get_vertical_edge_index(L, i, j);
        int edge2 = get_vertical_edge_index(L, (i + 1) % L, j);
        if (lattice->edge_types[edge1] == lattice->edge_types[edge2]) {
          v_violations++;
        }
      }
    }

    printf("  Alternation violations: H=%d, V=%d (%.1f%%)\n", h_violations,
           v_violations, 100.0 * (h_violations + v_violations) / (2 * L * L));
  } else if (mode == MODE_UNCORRELATED) {
    printf("  Per-node range: [%d, %d] type-1 edges\n", min_degree, max_degree);
    printf("  Expected binomial std: %.3f\n", sqrt(4 * f * (1 - f)));
  }

  free(node_type1_count);
}