#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Configuration structure */
typedef struct {
  char model[32]; /* "powerlaw" or "heaviside" */
  double sigma_min;
  double sigma_max;
  int sigma_steps;
  double alpha_min; /* For powerlaw */
  double alpha_max;
  int alpha_steps;
  double lambda_min; /* For heaviside */
  double lambda_max;
  int lambda_steps;
  double omega_min;
  double omega_max;
  int omega_steps;
  long n_samples;
} Config;

/* Statistics structure */
typedef struct {
  double mean_omega;
  double var_omega;
  double cov_collinear;
  double rho_collinear;
  double cov_perp;
  double rho_perp;
} Statistics;

/* Function to apply power-law transformation */
static inline double omega_powerlaw(double d, double alpha) {
  return pow(d, alpha);
}

/* Function to apply Heaviside transformation */
static inline double omega_heaviside(double d, double lambda, double Omega) {
  return (d > lambda) ? 1.0 : Omega;
}

/* Compute distance between two points */
static inline double distance(double x1, double y1, double x2, double y2) {
  double dx = x2 - x1;
  double dy = y2 - y1;
  return sqrt(dx * dx + dy * dy);
}

/* Compute statistics using Welford's online algorithm for O(1) memory and
 * numerical stability */
void compute_statistics_powerlaw(double sigma, double alpha, long n_samples,
                                 gsl_rng *rng, Statistics *stats) {
  /* Welford's algorithm for pooled omega statistics */
  long n_omega = 0;
  double mean_omega = 0.0;
  double m2_omega = 0.0;

  /* Welford's algorithm for individual edge statistics */
  long n = 0;
  double mean_ab = 0.0, mean_bc = 0.0, mean_bd = 0.0;
  double m2_ab = 0.0, m2_bc = 0.0, m2_bd = 0.0;
  double c_ab_bc = 0.0, c_ab_bd = 0.0;

  for (long i = 0; i < n_samples; i++) {
    /* Generate four points: A=(0,0), B=(1,0), C=(2,0), D=(1,1) with Gaussian
     * noise */
    double ax = gsl_ran_gaussian(rng, sigma);
    double ay = gsl_ran_gaussian(rng, sigma);
    double bx = 1.0 + gsl_ran_gaussian(rng, sigma);
    double by = gsl_ran_gaussian(rng, sigma);
    double cx = 2.0 + gsl_ran_gaussian(rng, sigma);
    double cy = gsl_ran_gaussian(rng, sigma);
    double dx = 1.0 + gsl_ran_gaussian(rng, sigma);
    double dy = 1.0 + gsl_ran_gaussian(rng, sigma);

    /* Compute distances */
    double d_ab = distance(ax, ay, bx, by);
    double d_bc = distance(bx, by, cx, cy);
    double d_bd = distance(bx, by, dx, dy);

    /* Apply transformation */
    double w_ab = omega_powerlaw(d_ab, alpha);
    double w_bc = omega_powerlaw(d_bc, alpha);
    double w_bd = omega_powerlaw(d_bd, alpha);

    /* Update pooled statistics for all omega values */
    n_omega++;
    double delta = w_ab - mean_omega;
    mean_omega += delta / n_omega;
    m2_omega += delta * (w_ab - mean_omega);

    n_omega++;
    delta = w_bc - mean_omega;
    mean_omega += delta / n_omega;
    m2_omega += delta * (w_bc - mean_omega);

    n_omega++;
    delta = w_bd - mean_omega;
    mean_omega += delta / n_omega;
    m2_omega += delta * (w_bd - mean_omega);

    /* Update statistics for correlations (standard Welford) */
    n++;

    double delta_ab = w_ab - mean_ab;
    double delta_bc = w_bc - mean_bc;
    double delta_bd = w_bd - mean_bd;

    mean_ab += delta_ab / n;
    mean_bc += delta_bc / n;
    mean_bd += delta_bd / n;

    double delta2_ab = w_ab - mean_ab;
    double delta2_bc = w_bc - mean_bc;
    double delta2_bd = w_bd - mean_bd;

    m2_ab += delta_ab * delta2_ab;
    m2_bc += delta_bc * delta2_bc;
    m2_bd += delta_bd * delta2_bd;

    c_ab_bc += delta_ab * delta2_bc;
    c_ab_bd += delta_ab * delta2_bd;
  }

  /* Finalize statistics */
  stats->mean_omega = mean_omega;
  stats->var_omega = m2_omega / n_omega;

  double var_ab = m2_ab / n;
  double var_bc = m2_bc / n;
  double var_bd = m2_bd / n;

  stats->cov_collinear = c_ab_bc / n;
  stats->cov_perp = c_ab_bd / n;

  /* Compute correlation coefficients with safety checks */
  double denom_collinear = sqrt(var_ab * var_bc);
  double denom_perp = sqrt(var_ab * var_bd);

  stats->rho_collinear = (denom_collinear > 1e-15)
                             ? (stats->cov_collinear / denom_collinear)
                             : 0.0;
  stats->rho_perp = (denom_perp > 1e-15) ? (stats->cov_perp / denom_perp) : 0.0;
}

void compute_statistics_heaviside(double sigma, double lambda, double Omega,
                                  long n_samples, gsl_rng *rng,
                                  Statistics *stats) {
  /* Allocate arrays for storing samples */
  double *omega_all = malloc(3 * n_samples * sizeof(double));
  double *omega_ab = malloc(n_samples * sizeof(double));
  double *omega_bc = malloc(n_samples * sizeof(double));
  double *omega_bd = malloc(n_samples * sizeof(double));

  if (!omega_all || !omega_ab || !omega_bc || !omega_bd) {
    fprintf(stderr, "Error: Memory allocation failed\n");
    exit(1);
  }

  long total_omega_samples = 0;

  for (long i = 0; i < n_samples; i++) {
    /* Generate four points with Gaussian noise */
    double ax = gsl_ran_gaussian(rng, sigma);
    double ay = gsl_ran_gaussian(rng, sigma);
    double bx = 1.0 + gsl_ran_gaussian(rng, sigma);
    double by = gsl_ran_gaussian(rng, sigma);
    double cx = 2.0 + gsl_ran_gaussian(rng, sigma);
    double cy = gsl_ran_gaussian(rng, sigma);
    double dx = 1.0 + gsl_ran_gaussian(rng, sigma);
    double dy = 1.0 + gsl_ran_gaussian(rng, sigma);

    /* Compute distances */
    double d_ab = distance(ax, ay, bx, by);
    double d_bc = distance(bx, by, cx, cy);
    double d_bd = distance(bx, by, dx, dy);

    /* Apply transformation */
    double w_ab = omega_heaviside(d_ab, lambda, Omega);
    double w_bc = omega_heaviside(d_bc, lambda, Omega);
    double w_bd = omega_heaviside(d_bd, lambda, Omega);

    /* Store samples */
    omega_all[total_omega_samples++] = w_ab;
    omega_all[total_omega_samples++] = w_bc;
    omega_all[total_omega_samples++] = w_bd;

    omega_ab[i] = w_ab;
    omega_bc[i] = w_bc;
    omega_bd[i] = w_bd;
  }

  /* Two-pass algorithm for mean */
  double sum_omega = 0.0;
  for (long i = 0; i < total_omega_samples; i++) {
    sum_omega += omega_all[i];
  }
  stats->mean_omega = sum_omega / total_omega_samples;

  /* Second pass for variance */
  double sum_sq_dev = 0.0;
  for (long i = 0; i < total_omega_samples; i++) {
    double dev = omega_all[i] - stats->mean_omega;
    sum_sq_dev += dev * dev;
  }
  stats->var_omega = sum_sq_dev / total_omega_samples;

  /* Compute means for correlation calculation */
  double mean_ab = 0.0, mean_bc = 0.0, mean_bd = 0.0;
  for (long i = 0; i < n_samples; i++) {
    mean_ab += omega_ab[i];
    mean_bc += omega_bc[i];
    mean_bd += omega_bd[i];
  }
  mean_ab /= n_samples;
  mean_bc /= n_samples;
  mean_bd /= n_samples;

  /* Compute variances and covariances */
  double var_ab = 0.0, var_bc = 0.0, var_bd = 0.0;
  double cov_ab_bc = 0.0, cov_ab_bd = 0.0;

  for (long i = 0; i < n_samples; i++) {
    double dev_ab = omega_ab[i] - mean_ab;
    double dev_bc = omega_bc[i] - mean_bc;
    double dev_bd = omega_bd[i] - mean_bd;

    var_ab += dev_ab * dev_ab;
    var_bc += dev_bc * dev_bc;
    var_bd += dev_bd * dev_bd;

    cov_ab_bc += dev_ab * dev_bc;
    cov_ab_bd += dev_ab * dev_bd;
  }

  var_ab /= n_samples;
  var_bc /= n_samples;
  var_bd /= n_samples;
  cov_ab_bc /= n_samples;
  cov_ab_bd /= n_samples;

  stats->cov_collinear = cov_ab_bc;
  stats->cov_perp = cov_ab_bd;

  /* Compute correlation coefficients with safety checks */
  double denom_collinear = sqrt(var_ab * var_bc);
  double denom_perp = sqrt(var_ab * var_bd);

  stats->rho_collinear =
      (denom_collinear > 1e-15) ? (cov_ab_bc / denom_collinear) : 0.0;
  stats->rho_perp = (denom_perp > 1e-15) ? (cov_ab_bd / denom_perp) : 0.0;

  /* Free allocated memory */
  free(omega_all);
  free(omega_ab);
  free(omega_bc);
  free(omega_bd);
}

/* Generate output filename */
void generate_filename(const Config *cfg, char *filename, size_t maxlen) {
  if (strcmp(cfg->model, "powerlaw") == 0) {
    snprintf(filename, maxlen,
             "results_powerlaw_sigma_%.3f_%.3f_alpha_%.2f_%.2f_N_%ld.dat",
             cfg->sigma_min, cfg->sigma_max, cfg->alpha_min, cfg->alpha_max,
             cfg->n_samples);
  } else {
    snprintf(filename, maxlen,
             "results_heaviside_sigma_%.3f_%.3f_lambda_%.3f_%.3f_omega_%.3f_%."
             "3f_N_%ld.dat",
             cfg->sigma_min, cfg->sigma_max, cfg->lambda_min, cfg->lambda_max,
             cfg->omega_min, cfg->omega_max, cfg->n_samples);
  }
}

/* Print usage */
void print_usage(const char *progname) {
  fprintf(stderr, "Usage: %s [OPTIONS]\n\n", progname);
  fprintf(stderr, "Required options:\n");
  fprintf(stderr,
          "  --model MODEL          Model type: 'powerlaw' or 'heaviside'\n");
  fprintf(stderr, "  --sigma-min VALUE      Minimum sigma value\n");
  fprintf(stderr, "  --sigma-max VALUE      Maximum sigma value\n");
  fprintf(stderr, "  --sigma-steps N        Number of sigma steps\n");
  fprintf(stderr, "  --n-samples N          Number of Monte Carlo samples per "
                  "parameter set\n\n");
  fprintf(stderr, "Power-law model options:\n");
  fprintf(stderr, "  --alpha-min VALUE      Minimum alpha value\n");
  fprintf(stderr, "  --alpha-max VALUE      Maximum alpha value\n");
  fprintf(stderr, "  --alpha-steps N        Number of alpha steps\n\n");
  fprintf(stderr, "Heaviside model options:\n");
  fprintf(stderr, "  --lambda-min VALUE     Minimum lambda value\n");
  fprintf(stderr, "  --lambda-max VALUE     Maximum lambda value\n");
  fprintf(stderr, "  --lambda-steps N       Number of lambda steps\n");
  fprintf(stderr, "  --omega-min VALUE      Minimum Omega value\n");
  fprintf(stderr, "  --omega-max VALUE      Maximum Omega value\n");
  fprintf(stderr, "  --omega-steps N        Number of Omega steps\n");
}

int main(int argc, char *argv[]) {
  Config cfg = {0};

  /* Default values */
  cfg.sigma_steps = 10;
  cfg.alpha_steps = 10;
  cfg.lambda_steps = 10;
  cfg.omega_steps = 10;
  cfg.n_samples = 1000000;

  /* Long options for getopt_long */
  static struct option long_options[] = {
      {"model", required_argument, 0, 'm'},
      {"sigma-min", required_argument, 0, 's'},
      {"sigma-max", required_argument, 0, 'S'},
      {"sigma-steps", required_argument, 0, 'n'},
      {"alpha-min", required_argument, 0, 'a'},
      {"alpha-max", required_argument, 0, 'A'},
      {"alpha-steps", required_argument, 0, 'p'},
      {"lambda-min", required_argument, 0, 'l'},
      {"lambda-max", required_argument, 0, 'L'},
      {"lambda-steps", required_argument, 0, 'q'},
      {"omega-min", required_argument, 0, 'o'},
      {"omega-max", required_argument, 0, 'O'},
      {"omega-steps", required_argument, 0, 'r'},
      {"n-samples", required_argument, 0, 'N'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};

  int opt;
  int option_index = 0;

  while ((opt = getopt_long(argc, argv, "hm:s:S:n:a:A:p:l:L:q:o:O:r:N:",
                            long_options, &option_index)) != -1) {
    switch (opt) {
    case 'm':
      strncpy(cfg.model, optarg, sizeof(cfg.model) - 1);
      break;
    case 's':
      cfg.sigma_min = atof(optarg);
      break;
    case 'S':
      cfg.sigma_max = atof(optarg);
      break;
    case 'n':
      cfg.sigma_steps = atoi(optarg);
      break;
    case 'a':
      cfg.alpha_min = atof(optarg);
      break;
    case 'A':
      cfg.alpha_max = atof(optarg);
      break;
    case 'p':
      cfg.alpha_steps = atoi(optarg);
      break;
    case 'l':
      cfg.lambda_min = atof(optarg);
      break;
    case 'L':
      cfg.lambda_max = atof(optarg);
      break;
    case 'q':
      cfg.lambda_steps = atoi(optarg);
      break;
    case 'o':
      cfg.omega_min = atof(optarg);
      break;
    case 'O':
      cfg.omega_max = atof(optarg);
      break;
    case 'r':
      cfg.omega_steps = atoi(optarg);
      break;
    case 'N':
      cfg.n_samples = atol(optarg);
      break;
    case 'h':
    default:
      print_usage(argv[0]);
      return 1;
    }
  }

  /* Validate inputs */
  if (strlen(cfg.model) == 0) {
    fprintf(stderr, "Error: --model is required\n");
    print_usage(argv[0]);
    return 1;
  }

  if (strcmp(cfg.model, "powerlaw") != 0 &&
      strcmp(cfg.model, "heaviside") != 0) {
    fprintf(stderr, "Error: model must be 'powerlaw' or 'heaviside'\n");
    return 1;
  }

  /* Initialize GSL RNG with high-quality generator */
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, time(NULL));

  /* Generate output filename */
  char filename[512];
  generate_filename(&cfg, filename, sizeof(filename));

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Error: Cannot open file %s for writing\n", filename);
    gsl_rng_free(rng);
    return 1;
  }

  /* Write header */
  fprintf(fp, "# Lattice frequency correlation analysis\n");
  fprintf(fp, "# Model: %s\n", cfg.model);
  fprintf(fp, "# Sigma range: [%.6f, %.6f], steps: %d\n", cfg.sigma_min,
          cfg.sigma_max, cfg.sigma_steps);
  fprintf(fp, "# N_samples: %ld\n", cfg.n_samples);

  if (strcmp(cfg.model, "powerlaw") == 0) {
    fprintf(fp, "# Alpha range: [%.6f, %.6f], steps: %d\n", cfg.alpha_min,
            cfg.alpha_max, cfg.alpha_steps);
    fprintf(fp, "# Columns: sigma alpha mean_omega var_omega cov_collinear "
                "rho_collinear cov_perp rho_perp\n");

    /* Parameter sweep */
    for (int i = 0; i < cfg.sigma_steps; i++) {
      double sigma;
      if (cfg.sigma_steps == 1) {
        sigma = cfg.sigma_min;
      } else {
        sigma = cfg.sigma_min +
                i * (cfg.sigma_max - cfg.sigma_min) / (cfg.sigma_steps - 1);
      }

      for (int j = 0; j < cfg.alpha_steps; j++) {
        double alpha;
        if (cfg.alpha_steps == 1) {
          alpha = cfg.alpha_min;
        } else {
          alpha = cfg.alpha_min +
                  j * (cfg.alpha_max - cfg.alpha_min) / (cfg.alpha_steps - 1);
        }

        Statistics stats;
        compute_statistics_powerlaw(sigma, alpha, cfg.n_samples, rng, &stats);

        fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", sigma, alpha,
                stats.mean_omega, stats.var_omega, stats.cov_collinear,
                stats.rho_collinear, stats.cov_perp, stats.rho_perp);

        printf("Progress: sigma=%.4f, alpha=%.4f\n", sigma, alpha);
      }
      fprintf(fp, "\n"); /* Blank line for gnuplot pm3d */
    }
  } else { /* heaviside */
    fprintf(fp, "# Lambda range: [%.6f, %.6f], steps: %d\n", cfg.lambda_min,
            cfg.lambda_max, cfg.lambda_steps);
    fprintf(fp, "# Omega range: [%.6f, %.6f], steps: %d\n", cfg.omega_min,
            cfg.omega_max, cfg.omega_steps);
    fprintf(fp, "# Columns: sigma lambda Omega mean_omega var_omega "
                "cov_collinear rho_collinear cov_perp rho_perp\n");

    /* Parameter sweep */
    for (int i = 0; i < cfg.sigma_steps; i++) {
      double sigma;
      if (cfg.sigma_steps == 1) {
        sigma = cfg.sigma_min;
      } else {
        sigma = cfg.sigma_min +
                i * (cfg.sigma_max - cfg.sigma_min) / (cfg.sigma_steps - 1);
      }

      for (int j = 0; j < cfg.lambda_steps; j++) {
        double lambda;
        if (cfg.lambda_steps == 1) {
          lambda = cfg.lambda_min;
        } else {
          lambda = cfg.lambda_min + j * (cfg.lambda_max - cfg.lambda_min) /
                                        (cfg.lambda_steps - 1);
        }

        for (int k = 0; k < cfg.omega_steps; k++) {
          double Omega;
          if (cfg.omega_steps == 1) {
            Omega = cfg.omega_min;
          } else {
            Omega = cfg.omega_min +
                    k * (cfg.omega_max - cfg.omega_min) / (cfg.omega_steps - 1);
          }

          Statistics stats;
          compute_statistics_heaviside(sigma, lambda, Omega, cfg.n_samples, rng,
                                       &stats);

          fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", sigma,
                  lambda, Omega, stats.mean_omega, stats.var_omega,
                  stats.cov_collinear, stats.rho_collinear, stats.cov_perp,
                  stats.rho_perp);

          printf("Progress: sigma=%.4f, lambda=%.4f, Omega=%.4f\n", sigma,
                 lambda, Omega);
        }
      }
    }
  }

  fclose(fp);
  gsl_rng_free(rng);

  printf("\nResults written to: %s\n", filename);

  return 0;
}