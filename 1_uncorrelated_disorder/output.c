#include "output.h"

char* generate_filename(const SimulationConfig* config) {
    char* filename = malloc(MAX_FILENAME);
    if (!filename) return NULL;
    
    if (config->dist_type == DIST_GAUSSIAN) {
        snprintf(filename, MAX_FILENAME,
                "square_uncorrelated_gaussian_L%d_mu%.2lf_sigma%.2lf_Ti%.1f_Tf%.1f_M%d_D%d_seed%u.dat",
                config->L, config->param1, config->param2, config->T_INITIAL, config->T_FINAL, 
                config->M, config->D, config->seed);
    } else { // DIST_BERNOULLI
        snprintf(filename, MAX_FILENAME,
                "square_uncorrelated_bernoulli_L%d_p%.4lf_Omega%.2lf_Ti%.1f_Tf%.1f_M%d_D%d_seed%u.dat",
                config->L, config->param1, config->param2, config->T_INITIAL, config->T_FINAL,
                config->M, config->D, config->seed);
    }
    
    return filename;
}

void write_results(const char* filename, const SimulationConfig* config,
                   TimeStepStats* stats, int num_time_steps) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        printf("Error opening output file\n");
        return;
    }
    
    // Write header
    fprintf(fp, "# Random seed: %u\n", config->seed);
    
    if (config->dist_type == DIST_GAUSSIAN) {
        fprintf(fp, "# Distribution: Gaussian(mu=%.2f, sigma=%.2f)\n", 
                config->param1, config->param2);
    } else {
        fprintf(fp, "# Distribution: Bernoulli(p=%.4f: freq=1, prob=%.4f: freq=%.2f)\n",
                config->param1, 1.0 - config->param1, config->param2);
    }
    
    fprintf(fp, "# M=%d disorder realizations, D=%d percolation processes per time\n",
            config->M, config->D);
    fprintf(fp, "# Structure: Generate disorder -> Complete time evolution with D runs per time -> Repeat M times\n");
    fprintf(fp, "# Activation: p(omega,t) = 1 - |cos(omega*t)|\n");
    fprintf(fp, "# L=%d, T_INITIAL=%.2f, T_FINAL=%.2f, dt=%.4f\n", 
            config->L, config->T_INITIAL, config->T_FINAL, config->dt);
    fprintf(fp, "# Susceptibility formula: sum(s^2*n_s) / sum(s*n_s) excluding largest cluster\n");
    fprintf(fp, "# t avg_prob active_frac active_std largest_frac largest_std suscept_avg suscept_std\n");
    
    // Write data
    for (int i = 0; i < num_time_steps; i++) {
        fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
                stats[i].t, stats[i].avg_prob,
                stats[i].mean_active, stats[i].std_active,
                stats[i].mean_largest, stats[i].std_largest,
                stats[i].mean_suscept, stats[i].std_suscept);
    }
    
    fclose(fp);
}
