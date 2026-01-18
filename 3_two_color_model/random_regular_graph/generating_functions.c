#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TMAX 6.29

// Global variables for parameters
double alpha;
int k;

// Function prototypes
double p1(double t) { return 1 - fabs(cos(t)); }
double p2(double t) { return 1 - fabs(cos(alpha * t)); }

// Solve the coupled system m1-m2 using fixed point iteration
void solve_coupled_m1_m2(double p1_val, double p2_val, int k1, int k2,
                         double *m1, double *m2) {
  double m1_old = 0.1, m2_old = 0.1;
  double m1_new, m2_new;
  double tol = 1e-25;
  unsigned int max_iter = 1000000000;
  int iter = 0;

  do {
    // Updated recursive equations for general k1, k2
    if (k1 > 0 && k2 > 0) {
      m1_new = 1 - pow((1 - p1_val * m1_old), k1 - 1) *
                       pow((1 - p2_val * m2_old), k2);
      m2_new = 1 - pow((1 - p1_val * m1_old), k1) *
                       pow((1 - p2_val * m2_old), k2 - 1);
    } else if (k1 == 0) {
      // Pure alpha case: only weight-2 edges
      m1_new = 0.0; // No weight-1 edges
      m2_new = 1 - pow((1 - p2_val * m2_old), k - 1);
    } else if (k2 == 0) {
      // Pure weight-1 case: only weight-1 edges
      m1_new = 1 - pow((1 - p1_val * m1_old), k - 1);
      m2_new = 0.0; // No weight-2 edges
    } else {
      // This shouldn't happen for valid k1, k2
      m1_new = m2_new = 0.0;
    }

    // Check convergence for both variables
    if (fabs(m1_new - m1_old) < tol && fabs(m2_new - m2_old) < tol) {
      break;
    }

    m1_old = m1_new;
    m2_old = m2_new;
    iter++;
  } while (iter < max_iter);

  *m1 = m1_new;
  *m2 = m2_new;
}

// Solve for m_bar(t) - this one is independent
double solve_m_bar(double p_bar_val) {
  double m_old = 0.1, m_new;
  double tol = 1e-25;
  int max_iter = 1000000000;
  int iter = 0;

  do {
    m_new = 1 - pow(1 - p_bar_val * m_old, k - 1);
    if (fabs(m_new - m_old) < tol)
      break;
    m_old = m_new;
    iter++;
  } while (iter < max_iter);

  return m_new;
}

// Calculate S_corr(t) with updated formula
double S_corr(double p1_val, double m1_val, double p2_val, double m2_val,
              int k1, int k2) {
  if (k1 == 0) {
    // Pure alpha case
    return 1 - pow((1 - p2_val * m2_val), k);
  } else if (k2 == 0) {
    // Pure weight-1 case
    return 1 - pow((1 - p1_val * m1_val), k);
  } else {
    // General case
    return 1 - pow((1 - p1_val * m1_val), k1) * pow((1 - p2_val * m2_val), k2);
  }
}

// Calculate S_uncorr(t)
double S_uncorr(double p_bar_val, double m_bar_val) {
  return 1 - pow((1 - p_bar_val * m_bar_val), k);
}

void print_usage(const char *program_name) {
  printf("Usage: %s k alpha\n", program_name);
  printf("\nParameters:\n");
  printf("  k       - Degree of regular graph (integer > 0)\n");
  printf("  alpha   - Weight parameter for second edge type (double > 0)\n");
  printf("\nThe program will solve for all allowed f values: f = i/k for i = "
         "0, 1, ..., k\n");
  printf("Output files: GF_RRG_k{k}_alpha{alpha:.2f}_f{f:.2f}.dat\n");
  printf("\nExample:\n");
  printf("  %s 4 2.5\n", program_name);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Error: Expected 2 arguments, got %d\n", argc - 1);
    print_usage(argv[0]);
    return 1;
  }

  // Parse command-line arguments
  k = atoi(argv[1]);
  alpha = atof(argv[2]);

  // Validation
  if (k <= 0) {
    fprintf(stderr, "Error: Degree k=%d must be positive\n", k);
    return 1;
  }

  if (alpha <= 0) {
    fprintf(stderr, "Error: Alpha parameter (%.6f) must be positive\n", alpha);
    return 1;
  }

  printf("Random Regular Graph Generating Functions Analysis\n");
  printf("=================================================\n");
  printf("  Degree k: %d\n", k);
  printf("  Alpha parameter: %.6f\n", alpha);
  printf("  Time range: 0 to %.6f\n", TMAX);
  printf("  Allowed f values: ");

  for (int i = 0; i <= k; i++) {
    printf("%.3f ", (double)i / k);
  }
  printf("\n\n");

  // Loop over all allowed f values
  for (int i = 0; i <= k; i++) {
    double f = (double)i / k;
    int k1 = i;     // Number of weight-1 edges per node
    int k2 = k - i; // Number of weight-alpha edges per node

    printf("Processing f = %.3f (k1=%d, k2=%d)...\n", f, k1, k2);

    // Create output filename
    char filename[256];
    snprintf(filename, sizeof(filename), "GF_RRG_k%d_alpha%.2lf_n1%.d.dat", k,
             alpha, k1);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Error: Cannot open file %s\n", filename);
      continue;
    }

    // Write header with parameters
    fprintf(fp, "# Random Regular Graph Generating Functions\n");
    fprintf(fp, "# k=%d, alpha=%.6f, f=%.6f (k1=%d, k2=%d)\n", k, alpha, f, k1,
            k2);
    fprintf(fp, "# t p1 p2 p_bar m1 m2 m_bar S_corr S_uncorr\n");

    double t, p1_val, p2_val, p_bar_val;
    double m1_val, m2_val, m_bar_val;
    double s_corr_val, s_uncorr_val;

    for (t = 0.0; t <= TMAX; t += 0.0001) {
      // Calculate p values
      p1_val = p1(t);
      p2_val = p2(t);
      p_bar_val = f * p1_val + (1.0 - f) * p2_val; // Corrected formula

      // Solve the coupled system for m1 and m2
      solve_coupled_m1_m2(p1_val, p2_val, k1, k2, &m1_val, &m2_val);

      // Solve the independent equation for m_bar
      m_bar_val = solve_m_bar(p_bar_val);

      // Calculate S_corr and S_uncorr
      s_corr_val = S_corr(p1_val, m1_val, p2_val, m2_val, k1, k2);
      s_uncorr_val = S_uncorr(p_bar_val, m_bar_val);

      // Write results to file
      fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", t, p1_val,
              p2_val, p_bar_val, m1_val, m2_val, m_bar_val, s_corr_val,
              s_uncorr_val);
    }

    fclose(fp);
    printf("  Results saved to %s\n", filename);
  }

  printf("\nAll computations completed successfully!\n");
  printf("Generated %d output files for different f values.\n", k + 1);

  return 0;
}
