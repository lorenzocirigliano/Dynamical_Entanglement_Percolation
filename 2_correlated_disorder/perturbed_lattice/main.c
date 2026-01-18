#include "config.h"
#include "lattice.h"
#include "output.h"
#include "percolation.h"
#include "random_utils.h"
#include "union_find.h"
#include <getopt.h>
#include <time.h>

void print_usage(char *program_name) {
  printf("Usage: %s -L <size> -s <sigma> -w <weight_type> [options]\n",
         program_name);
  printf("\nRequired parameters:\n");
  printf("  -L <int>    : Lattice size\n");
  printf("  -s <float>  : Position perturbation std dev\n");
  printf("  -w <type>   : Weight function type (power|exp|bernoulli)\n");
  printf("\nWeight function parameters:\n");
  printf("  For 'power':     -a <float> : Power law exponent (default: 1.0)\n");
  printf("  For 'exp':       -l <float> : Exponential decay length\n");
  printf("                   -O <float> : Frequency Omega (prefactor)\n");
  printf("  For 'bernoulli': -l <float> : Threshold distance lambda\n");
  printf("                   -O <float> : Frequency Omega for d < lambda\n");
  printf("\nOptional parameters:\n");
  printf("  -I <float>  : Initial time (default: 0.0)\n");
  printf("  -F <float>  : Final time (default: 10.0)\n");
  printf("  -d <float>  : Time step dt (default: 0.01)\n");
  printf("  -M <int>    : M1 - activation realizations (default: 100)\n");
  printf("  -D <int>    : M2 - disorder realizations (default: 100)\n");
  printf("  -b          : Use open boundary conditions (default: PBC)\n");
  printf(
      "  -r <seed>   : RNG seed for reproducibility (default: time-based)\n");
  printf("  -h          : Show this help message\n");
  printf("\nExamples:\n");
  printf("  Power law:    %s -L 100 -s 0.1 -w power -a 2.0\n", program_name);
  printf("  Exponential:  %s -L 100 -s 0.1 -w exp -l 1.5 -O 1.0\n",
         program_name);
  printf("  Bernoulli:    %s -L 100 -s 0.1 -w bernoulli -l 1.0 -O 0.5\n",
         program_name);
}

int main(int argc, char *argv[]) {
  // Default configuration
  SimulationConfig config = {
      .L = 0,        // Must be specified
      .sigma = -1.0, // Must be specified
      .alpha = 1.0,  // Default power law exponent
      .lambda = 0.0, // Must be specified for exp/bernoulli
      .Omega = 0.0,  // Must be specified for bernoulli
      .weight_type = WEIGHT_POWER_LAW,
      .T_INITIAL = 0.0,
      .T_FINAL = 10.0,
      .dt = 0.01,
      .M1 = 100,
      .M2 = 100,
      .use_pbc = true,
      .seed = 0 // 0 means use time-based seed
  };

  char *weight_type_str = NULL;
  int opt;

  while ((opt = getopt(argc, argv, "L:s:w:a:l:O:I:F:d:M:D:br:h")) != -1) {
    switch (opt) {
    case 'L':
      config.L = atoi(optarg);
      break;
    case 's':
      config.sigma = atof(optarg);
      break;
    case 'w':
      weight_type_str = optarg;
      break;
    case 'a':
      config.alpha = atof(optarg);
      break;
    case 'l':
      config.lambda = atof(optarg);
      break;
    case 'O':
      config.Omega = atof(optarg);
      break;
    case 'I':
      config.T_INITIAL = atof(optarg);
      break;
    case 'F':
      config.T_FINAL = atof(optarg);
      break;
    case 'd':
      config.dt = atof(optarg);
      break;
    case 'M':
      config.M1 = atoi(optarg);
      break;
    case 'D':
      config.M2 = atoi(optarg);
      break;
    case 'b':
      config.use_pbc = false;
      break;
    case 'r':
      config.seed = strtoul(optarg, NULL, 10);
      break;
    case 'h':
      print_usage(argv[0]);
      return 0;
    default:
      print_usage(argv[0]);
      return 1;
    }
  }

  // Validate and set weight type
  if (!weight_type_str) {
    printf("Error: Weight function type (-w) must be specified\n\n");
    print_usage(argv[0]);
    return 1;
  }

  if (strcmp(weight_type_str, "power") == 0) {
    config.weight_type = WEIGHT_POWER_LAW;
  } else if (strcmp(weight_type_str, "exp") == 0) {
    config.weight_type = WEIGHT_EXPONENTIAL;
    if (config.lambda <= 0) {
      printf("Error: Exponential weight requires lambda (-l) > 0\n");
      return 1;
    }
    if (config.Omega <= 0) {
      printf("Error: Exponential weight requires Omega (-O) > 0\n");
      return 1;
    }
  } else if (strcmp(weight_type_str, "bernoulli") == 0) {
    config.weight_type = WEIGHT_BERNOULLI;
    if (config.lambda <= 0) {
      printf("Error: Bernoulli weight requires lambda (-l) > 0\n");
      return 1;
    }
    if (config.Omega <= 0) {
      printf("Error: Bernoulli weight requires Omega (-O) > 0\n");
      return 1;
    }
  } else {
    printf(
        "Error: Unknown weight type '%s'. Use 'power', 'exp', or 'bernoulli'\n",
        weight_type_str);
    return 1;
  }

  // Validate required parameters
  if (config.L <= 0 || config.sigma < 0) {
    printf("Error: L and sigma must be specified and valid\n\n");
    print_usage(argv[0]);
    return 1;
  }

  // Validate time parameters
  if (config.T_FINAL <= config.T_INITIAL) {
    printf("Error: T_FINAL must be greater than T_INITIAL\n");
    return 1;
  }

  // Set seed if not provided
  if (config.seed == 0) {
    config.seed = (unsigned long)time(NULL);
  }

  // Print configuration
  printf("Configuration:\n");
  printf("  Lattice: L=%d, %s boundary conditions\n", config.L,
         config.use_pbc ? "periodic" : "open");
  printf("  Disorder: sigma=%.3f\n", config.sigma);

  if (config.weight_type == WEIGHT_POWER_LAW) {
    printf("  Weight function: d^%.2f\n", config.alpha);
  } else if (config.weight_type == WEIGHT_EXPONENTIAL) {
    printf("  Weight function: %.2f*exp(-d/%.2f)\n", config.Omega,
           config.lambda);
  } else {
    printf("  Weight function: Bernoulli with lambda=%.2f, Omega=%.2f\n",
           config.lambda, config.Omega);
  }

  printf("  Time: T_INITIAL=%.2f, T_FINAL=%.2f, dt=%.3f\n", config.T_INITIAL,
         config.T_FINAL, config.dt);
  printf("  Statistics: M1=%d, M2=%d\n", config.M1, config.M2);
  printf("  RNG seed: %lu\n", config.seed);

  // Calculate number of time steps
  int num_time_steps =
      (int)((config.T_FINAL - config.T_INITIAL + config.dt / 2) / config.dt) +
      1;
  printf("  Time steps: %d\n", num_time_steps);

  // Initialize timing
  clock_t start_time = clock();

  // Initialize RNG with provided seed
  RNG rng;
  init_rng(&rng, config.seed);

  // Create lattice and structures
  Lattice *lattice = create_lattice(config.L, config.use_pbc);
  if (!lattice) {
    printf("Error creating lattice\n");
    return 1;
  }

  UnionFind *uf = create_union_find(config.L);
  if (!uf) {
    printf("Error creating union-find structure\n");
    destroy_lattice(lattice);
    return 1;
  }

  // Pre-allocate arrays for weight storage (avoid repeated allocations)
  int size = config.L * config.L;
  double *original_h = malloc(size * sizeof(double));
  double *original_v = malloc(size * sizeof(double));
  double *shuffled_h = malloc(size * sizeof(double));
  double *shuffled_v = malloc(size * sizeof(double));

  if (!original_h || !original_v || !shuffled_h || !shuffled_v) {
    printf("Error allocating weight storage arrays\n");
    free(original_h);
    free(original_v);
    free(shuffled_h);
    free(shuffled_v);
    destroy_lattice(lattice);
    destroy_union_find(uf);
    return 1;
  }

  // Allocate arrays to store results
  TimeStepResults *averaged_results =
      calloc(num_time_steps, sizeof(TimeStepResults));
  TimeStepResults *disorder_results =
      calloc(num_time_steps, sizeof(TimeStepResults));

  if (!averaged_results || !disorder_results) {
    printf("Error allocating result arrays\n");
    free(original_h);
    free(original_v);
    free(shuffled_h);
    free(shuffled_v);
    destroy_lattice(lattice);
    destroy_union_find(uf);
    return 1;
  }

  // Initialize time array
  for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
    averaged_results[t_idx].t = config.T_INITIAL + t_idx * config.dt;
  }

  printf("\nStarting simulation with %d disorder realizations...\n", config.M2);

  // Main loop: M2 disorder realizations
  for (int m2 = 0; m2 < config.M2; m2++) {
    if ((m2 + 1) % 10 == 0 || m2 == 0) {
      double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
      double estimated_total = elapsed * config.M2 / (m2 + 1);
      printf("  Disorder %d/%d (%.1f%%, est. %.1f s remaining)\n", m2 + 1,
             config.M2, 100.0 * (m2 + 1) / config.M2,
             estimated_total - elapsed);
    }

    // Generate disorder realization and compute original weights
    generate_disorder_realization(lattice, config.sigma, &rng);
    compute_edge_weights(lattice, &config);

    // Store original weights
    for (int i = 0; i < size; i++) {
      original_h[i] = lattice->h_edges[i].weight;
      original_v[i] = lattice->v_edges[i].weight;
    }

    // Create shuffled weights ONCE per disorder realization
    shuffle_weights(lattice, &rng);
    for (int i = 0; i < size; i++) {
      shuffled_h[i] = lattice->h_edges[i].weight;
      shuffled_v[i] = lattice->v_edges[i].weight;
    }

    // Clear disorder results for this realization
    for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
      disorder_results[t_idx].t = config.T_INITIAL + t_idx * config.dt;
      disorder_results[t_idx].avg_prob = 0;
      disorder_results[t_idx].correlated = (PercolationStats){0};
      disorder_results[t_idx].uncorrelated = (PercolationStats){0};
    }

    // Time evolution for this disorder realization
    for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
      double t = config.T_INITIAL + t_idx * config.dt;

      // Restore and precompute probabilities for CORRELATED version
      for (int i = 0; i < size; i++) {
        lattice->h_edges[i].weight = original_h[i];
        lattice->v_edges[i].weight = original_v[i];
      }
      precompute_probabilities(lattice, t, false);
      double avg_prob_corr = compute_average_probability(lattice, false);

      // Restore and precompute probabilities for UNCORRELATED version
      for (int i = 0; i < size; i++) {
        lattice->h_edges[i].weight = shuffled_h[i];
        lattice->v_edges[i].weight = shuffled_v[i];
      }
      precompute_probabilities(lattice, t, true);
      double avg_prob_uncorr = compute_average_probability(lattice, true);

      // Store average probability
      disorder_results[t_idx].avg_prob =
          (avg_prob_corr + avg_prob_uncorr) / 2.0;

      // M1 activation realizations
      for (int m1 = 0; m1 < config.M1; m1++) {
        // Correlated version
        PercolationStats corr_stats =
            simulate_percolation_with_probs(lattice, uf, &rng, false);
        disorder_results[t_idx].correlated.active_fraction +=
            corr_stats.active_fraction;
        disorder_results[t_idx].correlated.largest_cluster_fraction +=
            corr_stats.largest_cluster_fraction;
        disorder_results[t_idx].correlated.avg_cluster_size +=
            corr_stats.avg_cluster_size;

        // Uncorrelated version
        PercolationStats uncorr_stats =
            simulate_percolation_with_probs(lattice, uf, &rng, true);
        disorder_results[t_idx].uncorrelated.active_fraction +=
            uncorr_stats.active_fraction;
        disorder_results[t_idx].uncorrelated.largest_cluster_fraction +=
            uncorr_stats.largest_cluster_fraction;
        disorder_results[t_idx].uncorrelated.avg_cluster_size +=
            uncorr_stats.avg_cluster_size;
      }

      // Average over M1
      disorder_results[t_idx].correlated.active_fraction /= config.M1;
      disorder_results[t_idx].correlated.largest_cluster_fraction /= config.M1;
      disorder_results[t_idx].correlated.avg_cluster_size /= config.M1;
      disorder_results[t_idx].uncorrelated.active_fraction /= config.M1;
      disorder_results[t_idx].uncorrelated.largest_cluster_fraction /=
          config.M1;
      disorder_results[t_idx].uncorrelated.avg_cluster_size /= config.M1;
    }

    // Accumulate results from this disorder realization
    for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
      averaged_results[t_idx].avg_prob += disorder_results[t_idx].avg_prob;
      averaged_results[t_idx].correlated.active_fraction +=
          disorder_results[t_idx].correlated.active_fraction;
      averaged_results[t_idx].correlated.largest_cluster_fraction +=
          disorder_results[t_idx].correlated.largest_cluster_fraction;
      averaged_results[t_idx].correlated.avg_cluster_size +=
          disorder_results[t_idx].correlated.avg_cluster_size;
      averaged_results[t_idx].uncorrelated.active_fraction +=
          disorder_results[t_idx].uncorrelated.active_fraction;
      averaged_results[t_idx].uncorrelated.largest_cluster_fraction +=
          disorder_results[t_idx].uncorrelated.largest_cluster_fraction;
      averaged_results[t_idx].uncorrelated.avg_cluster_size +=
          disorder_results[t_idx].uncorrelated.avg_cluster_size;
    }
  }

  // Average over M2 disorder realizations
  printf("Averaging results...\n");
  for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
    averaged_results[t_idx].avg_prob /= config.M2;
    averaged_results[t_idx].correlated.active_fraction /= config.M2;
    averaged_results[t_idx].correlated.largest_cluster_fraction /= config.M2;
    averaged_results[t_idx].correlated.avg_cluster_size /= config.M2;
    averaged_results[t_idx].uncorrelated.active_fraction /= config.M2;
    averaged_results[t_idx].uncorrelated.largest_cluster_fraction /= config.M2;
    averaged_results[t_idx].uncorrelated.avg_cluster_size /= config.M2;
  }

  // Write results to file
  char *filename = generate_filename(&config);
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    printf("Error opening output file: %s\n", filename);
    free(filename);
    free(original_h);
    free(original_v);
    free(shuffled_h);
    free(shuffled_v);
    free(averaged_results);
    free(disorder_results);
    destroy_lattice(lattice);
    destroy_union_find(uf);
    return 1;
  }

  write_header(fp);

  for (int t_idx = 0; t_idx < num_time_steps; t_idx++) {
    const TimeStepResults *r = &averaged_results[t_idx];
    fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", r->t, r->avg_prob,
            r->correlated.active_fraction,
            r->correlated.largest_cluster_fraction,
            r->correlated.avg_cluster_size, r->uncorrelated.active_fraction,
            r->uncorrelated.largest_cluster_fraction,
            r->uncorrelated.avg_cluster_size);
  }

  fclose(fp);

  // Print summary
  double total_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
  printf("\nSimulation completed successfully in %.2f seconds!\n", total_time);
  printf("Results saved to: %s\n", filename);
  printf("Summary:\n");
  printf("  - %d disorder realizations\n", config.M2);
  printf("  - %d activation realizations per time step\n", config.M1);
  printf("  - %d time steps from %.2f to %.2f\n", num_time_steps,
         config.T_INITIAL, config.T_FINAL);
  int total_sims = 2 * num_time_steps * config.M2 * config.M1;
  printf("  - Total simulations: %d (%.0f sims/sec)\n", total_sims,
         total_sims / total_time);

  // Cleanup
  free(filename);
  free(original_h);
  free(original_v);
  free(shuffled_h);
  free(shuffled_v);
  free(averaged_results);
  free(disorder_results);
  destroy_lattice(lattice);
  destroy_union_find(uf);

  return 0;
}