#include "io.h"
#include "lattice.h"
#include "percolation.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Extended parameters structure
typedef struct {
  SimulationParams base;
  DisorderMode mode;
  char mode_name[32];
} ExtendedParams;

void print_extended_usage(const char *program_name) {
  printf("Usage: %s mode L M f p1_min p1_max p2_min p2_max delta_p seed\n",
         program_name);
  printf("\n");
  printf("  mode    : Disorder mode - 'correlated', 'alternating', or "
         "'uncorrelated'\n");
  printf("  L       : Lattice size\n");
  printf("  M       : Number of Monte Carlo realizations\n");
  printf("  f       : Fraction parameter (0-1, mapped to 0, 0.25, 0.5, 0.75, "
         "or 1.0)\n");
  printf("  p1_min  : Minimum probability for type 1 edges\n");
  printf("  p1_max  : Maximum probability for type 1 edges\n");
  printf("  p2_min  : Minimum probability for type 2 edges\n");
  printf("  p2_max  : Maximum probability for type 2 edges\n");
  printf("  delta_p : Step size for probability sweep\n");
  printf("  seed    : Random number generator seed (unsigned integer)\n");
  printf("\n");
  printf("Disorder modes:\n");
  printf("  correlated   : Each node has exactly f*4 type-1 edges\n");
  printf("  alternating  : Edges alternate type along each direction\n");
  printf(
      "  uncorrelated : Each edge independently type-1 with probability f\n");
  printf("\n");
  printf("Example:\n");
  printf("  %s alternating 50 1000 0.5 0.0 1.0 0.0 1.0 0.05 12345\n",
         program_name);
}

int parse_extended_command_line(int argc, char *argv[],
                                ExtendedParams *params) {
  if (argc != 11) {
    print_extended_usage(argv[0]);
    return 0;
  }

  // Parse disorder mode
  if (strcmp(argv[1], "correlated") == 0) {
    params->mode = MODE_CORRELATED;
    strcpy(params->mode_name, "correlated");
  } else if (strcmp(argv[1], "alternating") == 0) {
    params->mode = MODE_ALTERNATING;
    strcpy(params->mode_name, "alternating");
  } else if (strcmp(argv[1], "uncorrelated") == 0) {
    params->mode = MODE_UNCORRELATED;
    strcpy(params->mode_name, "uncorrelated");
  } else {
    fprintf(stderr, "Error: Invalid mode '%s'\n", argv[1]);
    fprintf(stderr,
            "Mode must be 'correlated', 'alternating', or 'uncorrelated'\n");
    return 0;
  }

  // Parse remaining parameters
  params->base.L = atoi(argv[2]);
  params->base.M = atoi(argv[3]);
  params->base.f = atof(argv[4]);
  params->base.p1_min = atof(argv[5]);
  params->base.p1_max = atof(argv[6]);
  params->base.p2_min = atof(argv[7]);
  params->base.p2_max = atof(argv[8]);
  params->base.delta_p = atof(argv[9]);
  params->base.seed = (unsigned int)atoi(argv[10]);

  // Map f to nearest standard value
  printf("Input f value: %.6f\n", params->base.f);
  if (params->base.f <= 0.125) {
    params->base.f = 0.0;
  } else if (params->base.f <= 0.375) {
    params->base.f = 0.25;
  } else if (params->base.f <= 0.625) {
    params->base.f = 0.5;
  } else if (params->base.f <= 0.875) {
    params->base.f = 0.75;
  } else {
    params->base.f = 1.0;
  }
  printf("Mapped f value: %.2f\n", params->base.f);

  // Validate parameters
  if (params->base.L <= 0 || params->base.M <= 0 || params->base.delta_p <= 0) {
    fprintf(stderr, "Error: L, M, and delta_p must be positive\n");
    return 0;
  }

  return 1;
}

FILE *create_extended_output_file(const ExtendedParams *params) {
  char filename[256];
  snprintf(filename, sizeof(filename), "phase_diagram_%s_L%d_M%d_f%.2lf.dat",
           params->mode_name, params->base.L, params->base.M, params->base.f);

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Error: Cannot create output file %s\n", filename);
    return NULL;
  }

  printf("Results will be written to %s\n", filename);
  return fp;
}

int main(int argc, char *argv[]) {
  ExtendedParams params;

  if (!parse_extended_command_line(argc, argv, &params)) {
    return 1;
  }

  // Initialize RNG with provided seed for reproducibility
  init_random_seed(params.base.seed);

  FILE *fp = create_extended_output_file(&params);
  if (!fp)
    return 1;

  // Write header with mode information
  fprintf(fp, "# Percolation simulation - %s disorder\n", params.mode_name);
  fprintf(fp, "# L=%d, M=%d, f=%.2f, seed=%u\n", params.base.L, params.base.M,
          params.base.f, params.base.seed);
  fprintf(fp, "# p1 p2 p_bar size_of_GC err_on_size_of_GC "
              "average_small_cluster_size err_on_average_small_cluster\n");

  printf("\n========================================\n");
  printf("PERCOLATION SIMULATION - %s DISORDER\n", params.mode_name);
  printf("========================================\n");
  printf("Parameters: L=%d, M=%d, f=%.2f, seed=%u\n", params.base.L,
         params.base.M, params.base.f, params.base.seed);
  printf("Probability ranges: p1=[%.2f,%.2f], p2=[%.2f,%.2f], step=%.3f\n",
         params.base.p1_min, params.base.p1_max, params.base.p2_min,
         params.base.p2_max, params.base.delta_p);

  // Create lattice and assign edge types ONCE
  Lattice *lattice = create_square_lattice_pbc(params.base.L);
  assign_edge_types(lattice, params.base.f, params.mode);

  printf("\nStarting parameter sweep...\n");

  // Count total points
  int total_points = 0;
  int completed_points = 0;

  for (double p1 = params.base.p1_min; p1 <= params.base.p1_max + 1e-9;
       p1 += params.base.delta_p) {
    for (double p2 = params.base.p2_min; p2 <= params.base.p2_max + 1e-9;
         p2 += params.base.delta_p) {
      total_points++;
    }
  }

  printf("Total parameter points: %d\n\n", total_points);

  // Main simulation loop
  for (double p1 = params.base.p1_min; p1 <= params.base.p1_max + 1e-9;
       p1 += params.base.delta_p) {
    for (double p2 = params.base.p2_min; p2 <= params.base.p2_max + 1e-9;
         p2 += params.base.delta_p) {

      double p_bar = params.base.f * p1 + (1.0 - params.base.f) * p2;

      double *gc_sizes = malloc(params.base.M * sizeof(double));
      double *small_sizes = malloc(params.base.M * sizeof(double));

      // Monte Carlo realizations
      for (int m = 0; m < params.base.M; m++) {
        activate_edges(lattice, p1, p2);
        ClusterStats stats = calculate_clusters(lattice);
        gc_sizes[m] = stats.giant_size;
        small_sizes[m] = stats.avg_small_size;
      }

      // Calculate statistics
      Statistics gc_stats = calculate_statistics(gc_sizes, params.base.M);
      Statistics small_stats = calculate_statistics(small_sizes, params.base.M);

      // Write results
      fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", p1, p2, p_bar,
              gc_stats.mean, gc_stats.std_error, small_stats.mean,
              small_stats.std_error);
      fflush(fp);

      free(gc_sizes);
      free(small_sizes);

      completed_points++;
      if (completed_points % 10 == 0 || completed_points == total_points) {
        printf("Progress: %d/%d (%.1f%%) - p1=%.3f, p2=%.3f\r",
               completed_points, total_points,
               100.0 * completed_points / total_points, p1, p2);
        fflush(stdout);
      }
    }
    fprintf(fp, "\n");
    fflush(fp);
  }

  printf("\n\nSimulation completed successfully!\n");

  // Clean up
  destroy_lattice(lattice);
  fclose(fp);

  return 0;
}