#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DELTA_P 1e-2

// Global parameters
int k;
double f;
int k1, k2;

// Function to evaluate the fixed point functions
void evaluate_functions(double p1, double p2, double m1, double m2, double *g1,
                        double *g2) {
  if (k1 > 0 && k2 > 0) {
    *g1 = 1.0 - pow(1.0 - p1 * m1, k1 - 1) * pow(1.0 - p2 * m2, k2);
    *g2 = 1.0 - pow(1.0 - p1 * m1, k1) * pow(1.0 - p2 * m2, k2 - 1);
  } else if (k1 == 0) {
    *g1 = 0.0;
    *g2 = 1.0 - pow(1.0 - p2 * m2, k - 1);
  } else if (k2 == 0) {
    *g1 = 1.0 - pow(1.0 - p1 * m1, k - 1);
    *g2 = 0.0;
  } else {
    *g1 = *g2 = 0.0;
  }
}

// Smart initial guess based on physical intuition
void get_smart_initial_guess(double p1, double p2, double *m1_init,
                             double *m2_init) {
  // Physical heuristic: start with small values, adjust for extreme cases
  *m1_init = 0.1;
  *m2_init = 0.1;

  // For extreme probability cases, use better guesses
  if (p1 < 0.01 && p2 > 0.5) {
    *m1_init = 0.01;
    *m2_init = fmin(0.8, p2);
  } else if (p2 < 0.01 && p1 > 0.5) {
    *m1_init = fmin(0.8, p1);
    *m2_init = 0.01;
  } else if (p1 > 0.8 && p2 > 0.8) {
    *m1_init = 0.7;
    *m2_init = 0.7;
  }
}

// Robust solver with multiple strategies
int solve_robust_m1_m2(double p1, double p2, double *m1_final,
                       double *m2_final) {
  const double tolerance = 1e-10;
  const int max_iterations = 5000;

  // Try multiple initial conditions and strategies
  double initial_guesses[][2] = {
      {0.1, 0.1},   // Standard
      {0.5, 0.5},   // Medium
      {0.01, 0.01}, // Small
      {0.05, 0.05}, // Small-medium
      {0.0, 0.0}    // Will be replaced with smart guess
  };

  // Get smart guess
  get_smart_initial_guess(p1, p2, &initial_guesses[4][0],
                          &initial_guesses[4][1]);

  int num_strategies = 5;

  for (int strategy = 0; strategy < num_strategies; strategy++) {
    double m1 = initial_guesses[strategy][0];
    double m2 = initial_guesses[strategy][1];

    // Method 1: Simple fixed point with adaptive damping
    double damping = 0.5; // Start with moderate damping
    int consecutive_improvements = 0;
    double prev_error = 1e10;

    for (int iter = 0; iter < max_iterations; iter++) {
      double g1, g2;
      evaluate_functions(p1, p2, m1, m2, &g1, &g2);

      // Adaptive damping based on convergence behavior
      double m1_new = (1.0 - damping) * m1 + damping * g1;
      double m2_new = (1.0 - damping) * m2 + damping * g2;

      // Ensure bounds
      m1_new = fmax(0.0, fmin(1.0, m1_new));
      m2_new = fmax(0.0, fmin(1.0, m2_new));

      double error = fabs(m1_new - m1) + fabs(m2_new - m2);

      // Check convergence
      if (error < tolerance) {
        *m1_final = m1_new;
        *m2_final = m2_new;
        return 1; // Success
      }

      // Adaptive damping adjustment
      if (error < prev_error) {
        consecutive_improvements++;
        if (consecutive_improvements > 5 && damping < 0.9) {
          damping = fmin(0.9, damping * 1.1); // Increase damping
        }
      } else {
        consecutive_improvements = 0;
        if (damping > 0.1) {
          damping = fmax(0.1, damping * 0.8); // Decrease damping
        }
      }

      prev_error = error;
      m1 = m1_new;
      m2 = m2_new;
    }

    // Method 2: Try different damping for this initial guess
    if (strategy < 3) { // Only for first few strategies
      m1 = initial_guesses[strategy][0];
      m2 = initial_guesses[strategy][1];
      damping = 0.9; // High damping for stability

      for (int iter = 0; iter < max_iterations / 2; iter++) {
        double g1, g2;
        evaluate_functions(p1, p2, m1, m2, &g1, &g2);

        double m1_new = (1.0 - damping) * m1 + damping * g1;
        double m2_new = (1.0 - damping) * m2 + damping * g2;

        m1_new = fmax(0.0, fmin(1.0, m1_new));
        m2_new = fmax(0.0, fmin(1.0, m2_new));

        double error = fabs(m1_new - m1) + fabs(m2_new - m2);

        if (error < tolerance) {
          *m1_final = m1_new;
          *m2_final = m2_new;
          return 1; // Success
        }

        m1 = m1_new;
        m2 = m2_new;
      }
    }
  }

  // If all strategies failed, return the best guess we have
  *m1_final = 0.0;
  *m2_final = 0.0;
  return 0; // Failed to converge
}

// Solve for m_bar (this is usually stable)
double solve_m_bar(double p_bar_val) {
  double m_old = 0.1, m_new;
  double tol = 1e-12;
  int max_iter = 1000;

  for (int iter = 0; iter < max_iter; iter++) {
    m_new = 1.0 - pow(1.0 - p_bar_val * m_old, k - 1);
    if (fabs(m_new - m_old) < tol)
      break;
    m_old = m_new;
  }

  return m_new;
}

// Calculate S_corr
double S_corr(double p1_val, double m1_val, double p2_val, double m2_val) {
  if (k1 == 0) {
    return 1.0 - pow(1.0 - p2_val * m2_val, k);
  } else if (k2 == 0) {
    return 1.0 - pow(1.0 - p1_val * m1_val, k);
  } else {
    return 1.0 -
           pow(1.0 - p1_val * m1_val, k1) * pow(1.0 - p2_val * m2_val, k2);
  }
}

// Calculate S_uncorr
double S_uncorr(double p_bar_val, double m_bar_val) {
  return 1.0 - pow(1.0 - p_bar_val * m_bar_val, k);
}

void print_usage(const char *program_name) {
  printf("Usage: %s k f\n", program_name);
  printf("\nRobust Parameter Space Solver\n");
  printf("Parameters:\n");
  printf("  k  - Degree of regular graph\n");
  printf("  f  - Fraction of weight-1 edges (will be rounded to i/k)\n");
  printf("\nExample: %s 4 0.5\n", program_name);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    print_usage(argv[0]);
    return 1;
  }

  k = atoi(argv[1]);
  f = atof(argv[2]);

  if (k <= 0 || f < 0.0 || f > 1.0) {
    fprintf(stderr, "Error: Invalid parameters\n");
    return 1;
  }

  // Round f to allowed value
  int f_index = (int)round(f * k);
  f = (double)f_index / k;
  k1 = f_index;
  k2 = k - k1;

  printf("Robust Parameter Space Analysis\n");
  printf("==============================\n");
  printf("k=%d, f=%.6f (k1=%d, k2=%d)\n", k, f, k1, k2);
  printf("Grid step: %.0e\n", DELTA_P);

  // Create output filename
  char filename[256];
  snprintf(filename, sizeof(filename), "robust_param_space_k%d_f%.3lf.dat", k,
           f);

  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Error: Cannot open file %s\n", filename);
    return 1;
  }

  // Write header
  fprintf(fp, "# Robust Parameter Space Analysis\n");
  fprintf(fp, "# k=%d, f=%.6f (k1=%d, k2=%d)\n", k, f, k1, k2);
  fprintf(fp, "# p1 p2 average_p m1 m2 mbar S_corr S_uncorr converged\n");

  int total_points = 0;
  int converged_points = 0;

  printf("Computing parameter space...\n");

  // Loop over parameter space
  for (double p1 = 0.0; p1 <= 1.0 + DELTA_P / 2; p1 += DELTA_P) {
    if (p1 > 1.0)
      p1 = 1.0;

    for (double p2 = 0.0; p2 <= 1.0 + DELTA_P / 2; p2 += DELTA_P) {
      if (p2 > 1.0)
        p2 = 1.0;

      // Solve the coupled system
      double m1, m2;
      int converged = solve_robust_m1_m2(p1, p2, &m1, &m2);

      // Calculate other quantities
      double p_bar_val = f * p1 + (1.0 - f) * p2;
      double average_p = (p1 + p2) / 2.0;
      double m_bar_val = solve_m_bar(p_bar_val);

      // Calculate S values
      double s_corr_val = S_corr(p1, m1, p2, m2);
      double s_uncorr_val = S_uncorr(p_bar_val, m_bar_val);

      // Write to file
      fprintf(fp, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d\n", p1, p2,
              average_p, m1, m2, m_bar_val, s_corr_val, s_uncorr_val,
              converged);

      total_points++;
      if (converged)
        converged_points++;
    }

    // Progress indicator
    if (((int)(p1 / DELTA_P)) % 100 == 0) {
      printf("Progress: p1=%.3f, Success rate: %.1f%%\n", p1,
             100.0 * converged_points / total_points);
    }
  }

  fclose(fp);

  printf("\nComputation completed!\n");
  printf("Total points: %d\n", total_points);
  printf("Converged: %d (%.1f%%)\n", converged_points,
         100.0 * converged_points / total_points);
  printf("Results saved to: %s\n", filename);

  return 0;
}
