#include "config.h"
#include "lattice.h"
#include "union_find.h"
#include "percolation.h"
#include "random_utils.h"
#include "output.h"

void print_usage(char* program_name) {
    printf("Usage:\n");
    printf("  Gaussian:  %s gaussian L mu sigma M D T_INITIAL T_FINAL dt seed\n", program_name);
    printf("  Bernoulli: %s bernoulli L p Omega M D T_INITIAL T_FINAL dt seed\n", program_name);
    printf("\nParameters:\n");
    printf("  L: lattice size\n");
    printf("  Gaussian mode:\n");
    printf("    mu: mean frequency\n");
    printf("    sigma: frequency standard deviation\n");
    printf("  Bernoulli mode:\n");
    printf("    p: probability of frequency = 1\n");
    printf("    Omega: high frequency value (with probability 1-p)\n");
    printf("  M: number of disorder realizations\n");
    printf("  D: number of percolation processes per time per disorder\n");
    printf("  T_INITIAL: initial time\n");
    printf("  T_FINAL: final time\n");
    printf("  dt: time step\n");
    printf("  seed: random seed for reproducibility\n");
}

int main(int argc, char* argv[]) {
    // Check for correct number of arguments
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse configuration
    SimulationConfig config;
    
    // Read distribution type
    char* dist_type = argv[1];
    if (strcmp(dist_type, "gaussian") == 0) {
        config.dist_type = DIST_GAUSSIAN;
    } else if (strcmp(dist_type, "bernoulli") == 0) {
        config.dist_type = DIST_BERNOULLI;
    } else {
        printf("Error: Distribution type must be 'gaussian' or 'bernoulli'\n");
        return 1;
    }
    
    // Check argument count
    if (argc != 11) {
        printf("Error: Incorrect number of arguments\n\n");
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse parameters
    config.L = atoi(argv[2]);
    config.param1 = atof(argv[3]);
    config.param2 = atof(argv[4]);
    config.M = atoi(argv[5]);
    config.D = atoi(argv[6]);
    config.T_INITIAL = atof(argv[7]);
    config.T_FINAL = atof(argv[8]);
    config.dt = atof(argv[9]);
    config.seed = (unsigned int)atoi(argv[10]);
    
    // Validate time parameters
    if (config.T_FINAL <= config.T_INITIAL) {
        printf("Error: T_FINAL must be greater than T_INITIAL\n");
        return 1;
    }
    
    // Initialize random seed
    init_random(config.seed);
    
    int num_time_steps = (int)((config.T_FINAL - config.T_INITIAL + config.dt / 2) / config.dt) + 1;
    
    // Print simulation parameters
    printf("Starting simulation with:\n");
    printf("  Random seed: %u\n", config.seed);
    printf("  Lattice size: %d x %d\n", config.L, config.L);
    
    if (config.dist_type == DIST_GAUSSIAN) {
        printf("  Distribution: Gaussian\n");
        printf("  Frequency distribution: Gaussian(mu=%.2f, sigma=%.2f)\n",
               config.param1, config.param2);
    } else {
        printf("  Distribution: Bernoulli\n");
        printf("  Frequency distribution: 1 with prob=%.4f, %.2f with prob=%.4f\n",
               config.param1, config.param2, 1.0 - config.param1);
    }
    
    printf("  Activation: p(omega,t) = 1 - |cos(omega*t)|\n");
    printf("  Disorder realizations (M): %d\n", config.M);
    printf("  Percolation processes per time step (D): %d\n", config.D);
    printf("  Time range: [%.2f, %.2f] with dt=%.4f (%d time steps)\n",
           config.T_INITIAL, config.T_FINAL, config.dt, num_time_steps);
    printf("  Total percolation processes: %d\n", config.M * config.D * num_time_steps);
    
    // Allocate arrays for accumulating results
    int total_edges = 2 * config.L * config.L;
    int total_runs = config.M * config.D;
    
    double* avg_prob_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_active_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_active_sq_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_largest_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_largest_sq_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_suscept_vs_time = calloc(num_time_steps, sizeof(double));
    double* sum_suscept_sq_vs_time = calloc(num_time_steps, sizeof(double));
    
    // Create lattice and union-find structures
    Lattice* lattice = create_lattice(config.L);
    UnionFind* uf = create_union_find(config.L);
    
    if (!lattice || !uf) {
        printf("Error creating lattice or union-find structure\n");
        return 1;
    }
    
    // OUTER LOOP: Over disorder realizations (M)
    for (int m = 0; m < config.M; m++) {
        printf("\n=== Disorder realization %d/%d ===\n", m + 1, config.M);
        
        // Generate NEW random weights for this disorder realization
        generate_weights(lattice, &config);
        
        // Time evolution loop for THIS disorder realization
        int time_idx = 0;
        for (double t = config.T_INITIAL; t <= config.T_FINAL; t += config.dt, time_idx++) {
            
            // Calculate average probability for THIS disorder at THIS time
            double avg_prob = 0.0;
            for (int i = 0; i < config.L; i++) {
                for (int j = 0; j < config.L; j++) {
                    avg_prob += activation_probability(lattice->h_edges[i][j].weight, t);
                    avg_prob += activation_probability(lattice->v_edges[i][j].weight, t);
                }
            }
            avg_prob /= total_edges;
            avg_prob_vs_time[time_idx] += avg_prob;
            
            // INNER LOOP: Over percolation processes (D) for this time
            for (int d = 0; d < config.D; d++) {
                
                // Activate edges stochastically
                int active_edges = 0;
                for (int i = 0; i < config.L; i++) {
                    for (int j = 0; j < config.L; j++) {
                        // Horizontal edges
                        double prob_h = activation_probability(lattice->h_edges[i][j].weight, t);
                        lattice->h_edges[i][j].active = uniform_random() < prob_h;
                        if (lattice->h_edges[i][j].active)
                            active_edges++;
                        
                        // Vertical edges
                        double prob_v = activation_probability(lattice->v_edges[i][j].weight, t);
                        lattice->v_edges[i][j].active = uniform_random() < prob_v;
                        if (lattice->v_edges[i][j].active)
                            active_edges++;
                    }
                }
                
                // Find clusters
                find_clusters(uf, lattice);
                int largest_size;
                double suscept_avg;
                get_cluster_statistics(uf, &largest_size, &suscept_avg);
                
                // Compute observables
                double active_frac = (double)active_edges / total_edges;
                double largest_frac = (double)largest_size / (config.L * config.L);
                
                // Accumulate sums and sums of squares for statistics
                sum_active_vs_time[time_idx] += active_frac;
                sum_active_sq_vs_time[time_idx] += active_frac * active_frac;
                
                sum_largest_vs_time[time_idx] += largest_frac;
                sum_largest_sq_vs_time[time_idx] += largest_frac * largest_frac;
                
                sum_suscept_vs_time[time_idx] += suscept_avg;
                sum_suscept_sq_vs_time[time_idx] += suscept_avg * suscept_avg;
            }
            
            // Progress indicator every 10 time steps
            if ((time_idx + 1) % 10 == 0 || time_idx == num_time_steps - 1) {
                printf("\r  Time step %d/%d (t=%.4f)", time_idx + 1, num_time_steps, t);
                fflush(stdout);
            }
        }
        printf("\n");
    }
    
    // Compute statistics and prepare output
    TimeStepStats* stats = malloc(num_time_steps * sizeof(TimeStepStats));
    
    int time_idx = 0;
    for (double t = config.T_INITIAL; t <= config.T_FINAL; t += config.dt, time_idx++) {
        stats[time_idx].t = t;
        stats[time_idx].avg_prob = avg_prob_vs_time[time_idx] / config.M;
        
        // Compute means
        double mean_active = sum_active_vs_time[time_idx] / total_runs;
        double mean_largest = sum_largest_vs_time[time_idx] / total_runs;
        double mean_suscept = sum_suscept_vs_time[time_idx] / total_runs;
        
        // Compute standard deviations using: std = sqrt(E[X^2] - E[X]^2)
        double var_active = (sum_active_sq_vs_time[time_idx] / total_runs) -
                           (mean_active * mean_active);
        double std_active = (var_active > 0) ? sqrt(var_active) : 0.0;
        
        double var_largest = (sum_largest_sq_vs_time[time_idx] / total_runs) -
                            (mean_largest * mean_largest);
        double std_largest = (var_largest > 0) ? sqrt(var_largest) : 0.0;
        
        double var_suscept = (sum_suscept_sq_vs_time[time_idx] / total_runs) -
                            (mean_suscept * mean_suscept);
        double std_suscept = (var_suscept > 0) ? sqrt(var_suscept) : 0.0;
        
        stats[time_idx].mean_active = mean_active;
        stats[time_idx].std_active = std_active;
        stats[time_idx].mean_largest = mean_largest;
        stats[time_idx].std_largest = std_largest;
        stats[time_idx].mean_suscept = mean_suscept;
        stats[time_idx].std_suscept = std_suscept;
    }
    
    // Write results to file
    char* filename = generate_filename(&config);
    write_results(filename, &config, stats, num_time_steps);
    
    // Print summary
    printf("\n=== Simulation completed successfully! ===\n");
    printf("Results saved to %s\n", filename);
    printf("\nColumns:\n");
    printf("  1. t - time\n");
    printf("  2. avg_prob - average activation probability (averaged over M=%d disorder realizations)\n",
           config.M);
    printf("  3. active_frac - fraction of active edges (mean)\n");
    printf("  4. active_std - standard deviation of active_frac\n");
    printf("  5. largest_frac - largest cluster fraction (mean)\n");
    printf("  6. largest_std - standard deviation of largest_frac\n");
    printf("  7. suscept_avg - susceptibility-related quantity: sum(s^2*n_s)/sum(s*n_s) (mean)\n");
    printf("  8. suscept_std - standard deviation of suscept_avg\n");
    printf("\nAll statistics averaged over %d disorder realizations × %d percolation processes = %d total runs\n",
           config.M, config.D, config.M * config.D);
    
    // Cleanup
    free(filename);
    free(stats);
    free(avg_prob_vs_time);
    free(sum_active_vs_time);
    free(sum_active_sq_vs_time);
    free(sum_largest_vs_time);
    free(sum_largest_sq_vs_time);
    free(sum_suscept_vs_time);
    free(sum_suscept_sq_vs_time);
    destroy_lattice(lattice);
    destroy_union_find(uf);
    
    return 0;
}
