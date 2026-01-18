#include "io.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int parse_command_line(int argc, char *argv[], SimulationParams *params) {
  if (argc != 9) {
    print_usage(argv[0]);
    return 0;
  }

  params->L = atoi(argv[1]);
  params->M = atoi(argv[2]);
  params->f = atof(argv[3]);
  params->p1_min = atof(argv[4]);
  params->p1_max = atof(argv[5]);
  params->p2_min = atof(argv[6]);
  params->p2_max = atof(argv[7]);
  params->delta_p = atof(argv[8]);

  // Debug: Print the input f value
  printf("Input f value: %.6f\n", params->f);

  // Map f to the nearest standard value (0, 0.25, 0.5, 0.75, 1.0)
  if (params->f <= 0.125) {
    params->f = 0.0;
  } else if (params->f <= 0.375) {
    params->f = 0.25;
  } else if (params->f <= 0.625) {
    params->f = 0.5;
  } else if (params->f <= 0.875) {
    params->f = 0.75;
  } else {
    params->f = 1.0;
  }

  printf("Mapped f value: %.2f\n", params->f);

  // Validate parameters
  if (params->L <= 0 || params->M <= 0 || params->delta_p <= 0) {
    fprintf(stderr, "Error: L, M, and delta_p must be positive\n");
    return 0;
  }

  return 1;
}

FILE *create_output_file(const SimulationParams *params) {
  char filename[256];
  snprintf(filename, sizeof(filename), "phase_diagram_squareL%d_M%d_f%.2lf.dat",
           params->L, params->M, params->f);

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Error: Cannot create output file %s\n", filename);
    return NULL;
  }

  printf("Results will be written to %s\n", filename);
  return fp;
}

void write_header(FILE *fp) {
  fprintf(fp, "# p1 p2 p_bar size_of_GC err_on_size_of_GC "
              "average_small_cluster_size err_on_average_small_cluster\n");
}

void write_data_point(FILE *fp, double p1, double p2, double p_bar,
                      double gc_mean, double gc_err, double small_mean,
                      double small_err) {
  fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", p1, p2, p_bar, gc_mean,
          gc_err, small_mean, small_err);
  fflush(fp);
}

void print_usage(const char *program_name) {
  printf("Usage: %s L M f p1_min p1_max p2_min p2_max delta_p\n", program_name);
  printf("  L       : Lattice size\n");
  printf("  M       : Number of Monte Carlo realizations\n");
  printf("  f       : Fraction parameter (0-1, will be mapped to 0, 0.25, 0.5, "
         "0.75, or 1.0)\n");
  printf("  p1_min  : Minimum probability for type 1 edges\n");
  printf("  p1_max  : Maximum probability for type 1 edges\n");
  printf("  p2_min  : Minimum probability for type 2 edges\n");
  printf("  p2_max  : Maximum probability for type 2 edges\n");
  printf("  delta_p : Step size for probability sweep\n");
  printf("\n");
  printf("Standard f values and their meanings:\n");
  printf("  f = 0.00: All edges are type 2\n");
  printf("  f = 0.25: 1 type-1 edge per node, 3 type-2 edges per node\n");
  printf("  f = 0.50: 2 type-1 edges per node, 2 type-2 edges per node\n");
  printf("  f = 0.75: 3 type-1 edges per node, 1 type-2 edge per node\n");
  printf("  f = 1.00: All edges are type 1\n");
}
