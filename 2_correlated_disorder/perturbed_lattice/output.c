#include "output.h"

char *generate_filename(const SimulationConfig *config) {
    char *filename = malloc(MAX_FILENAME);
    if (!filename) return NULL;
    
    const char *bc_str = config->use_pbc ? "PBC" : "OBC";
    
    // Create different filename formats based on weight type
    switch (config->weight_type) {
        case WEIGHT_POWER_LAW:
            snprintf(filename, MAX_FILENAME,
                    "square_correlated_power_L%d_%s_sigma%.3f_alpha%.2f_Ti%.1f_Tf%.1f_dt%.3f_M1%d_M2%d_seed%lu.dat",
                    config->L, bc_str, config->sigma, config->alpha, 
                    config->T_INITIAL, config->T_FINAL, config->dt, config->M1, config->M2, config->seed);
            break;
            
        case WEIGHT_EXPONENTIAL:
            snprintf(filename, MAX_FILENAME,
                    "square_correlated_exp_L%d_%s_sigma%.3f_lambda%.2f_Omega%.2f_Ti%.1f_Tf%.1f_dt%.3f_M1%d_M2%d_seed%lu.dat",
                    config->L, bc_str, config->sigma, config->lambda, config->Omega,
                    config->T_INITIAL, config->T_FINAL, config->dt, config->M1, config->M2, config->seed);
            break;
            
        case WEIGHT_BERNOULLI:
            snprintf(filename, MAX_FILENAME,
                    "square_correlated_bernoulli_L%d_%s_sigma%.3f_lambda%.2f_Omega%.2f_Ti%.1f_Tf%.1f_dt%.3f_M1%d_M2%d_seed%lu.dat",
                    config->L, bc_str, config->sigma, config->lambda, config->Omega,
                    config->T_INITIAL, config->T_FINAL, config->dt, config->M1, config->M2, config->seed);
            break;
            
        default:
            snprintf(filename, MAX_FILENAME,
                    "square_correlated_unknown_L%d_%s_sigma%.3f_Ti%.1f_Tf%.1f_dt%.3f_M1%d_M2%d_seed%lu.dat",
                    config->L, bc_str, config->sigma, 
                    config->T_INITIAL, config->T_FINAL, config->dt, config->M1, config->M2, config->seed);
    }
    
    return filename;
}

void write_header(FILE *fp) {
    fprintf(fp, "# Percolation simulation results\n");
    fprintf(fp, "# Column 1: time t\n");
    fprintf(fp, "# Column 2: average activation probability\n");
    fprintf(fp, "# Columns 3-5: CORRELATED - p(t), P(t), S(t)\n");
    fprintf(fp, "#   3: p(t) = fraction of active edges\n");
    fprintf(fp, "#   4: P(t) = fraction of nodes in largest cluster\n");
    fprintf(fp, "#   5: S(t) = average cluster size (susceptibility-related: sum(n_s*s^2)/sum(n_s*s), excl. largest)\n");
    fprintf(fp, "# Columns 6-8: UNCORRELATED (shuffled) - p(t), P(t), S(t)\n");
    fprintf(fp, "#   6: p(t) = fraction of active edges\n");
    fprintf(fp, "#   7: P(t) = fraction of nodes in largest cluster\n");
    fprintf(fp, "#   8: S(t) = average cluster size (susceptibility-related, excl. largest)\n");
    fprintf(fp, "# Structure: M2 disorder realizations, M1 activation realizations per time\n");
    fprintf(fp, "# t avg_prob corr_p corr_P corr_S uncorr_p uncorr_P uncorr_S\n");
}