#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE 1024
#define GRID_SIZE 101  // 0 to 1 with spacing 0.01
#define GRID_STEP 0.01

// Structure to hold all observables at a grid point
typedef struct {
    double p_bar;
    double size_of_GC;
    double err_on_size_of_GC;
    double avg_small_cluster;
    double err_on_avg_small_cluster;
} Observables;

// Bilinear interpolation for a single observable
double interpolate(double **data, double x, double y) {
    // Clamp to [0, 1] range
    if (x < 0) x = 0;
    if (x > 1) x = 1;
    if (y < 0) y = 0;
    if (y > 1) y = 1;
    
    // Find grid indices
    double x_idx = x / GRID_STEP;
    double y_idx = y / GRID_STEP;
    
    int i0 = (int)x_idx;
    int j0 = (int)y_idx;
    int i1 = (i0 + 1 < GRID_SIZE) ? i0 + 1 : i0;
    int j1 = (j0 + 1 < GRID_SIZE) ? j0 + 1 : j0;
    
    // Interpolation weights
    double fx = x_idx - i0;
    double fy = y_idx - j0;
    
    // Bilinear interpolation
    double v00 = data[i0][j0];
    double v10 = data[i1][j0];
    double v01 = data[i0][j1];
    double v11 = data[i1][j1];
    
    return (1-fx)*(1-fy)*v00 + fx*(1-fy)*v10 + 
           (1-fx)*fy*v01 + fx*fy*v11;
}

// Interpolate all observables at point (x, y)
Observables interpolate_all(double **p_bar, double **size_gc, double **err_size_gc,
                             double **avg_small, double **err_avg_small,
                             double x, double y) {
    Observables obs;
    obs.p_bar = interpolate(p_bar, x, y);
    obs.size_of_GC = interpolate(size_gc, x, y);
    obs.err_on_size_of_GC = interpolate(err_size_gc, x, y);
    obs.avg_small_cluster = interpolate(avg_small, x, y);
    obs.err_on_avg_small_cluster = interpolate(err_avg_small, x, y);
    return obs;
}

// Allocate 2D array
double** alloc_2d_array() {
    double **arr = (double **)malloc(GRID_SIZE * sizeof(double *));
    for (int i = 0; i < GRID_SIZE; i++) {
        arr[i] = (double *)calloc(GRID_SIZE, sizeof(double));
    }
    return arr;
}

// Free 2D array
void free_2d_array(double **arr) {
    for (int i = 0; i < GRID_SIZE; i++) {
        free(arr[i]);
    }
    free(arr);
}

// Generate output filename
void generate_output_filename(const char *input, double omega, double t_max, 
                              double dt, char *output) {
    // Extract basename from input path
    const char *basename = strrchr(input, '/');
    if (basename == NULL) {
        basename = input;
    } else {
        basename++; // Skip the '/'
    }
    
    // Create a copy to work with
    char temp_name[MAX_LINE];
    strcpy(temp_name, basename);
    
    // Replace "phase_diagram" with "parametric"
    char *pos = strstr(temp_name, "phase_diagram");
    if (pos != NULL) {
        char temp[MAX_LINE];
        strcpy(temp, pos + strlen("phase_diagram"));
        sprintf(pos, "parametric_omega%.2f_tmax%.2f_dt%.4f%s", 
                omega, t_max, dt, temp);
        strcpy(output, temp_name);
    } else {
        // If "phase_diagram" not found, just append
        sprintf(output, "parametric_omega%.2f_tmax%.2f_dt%.4f_%s", 
                omega, t_max, dt, basename);
    }
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_file> <Omega> <T_max> <delta_t>\n", argv[0]);
        fprintf(stderr, "Example: %s data.dat 2.5 6.28 0.01\n", argv[0]);
        return 1;
    }
    
    const char *input_file = argv[1];
    double omega = atof(argv[2]);
    double t_max = atof(argv[3]);
    double dt = atof(argv[4]);
    
    printf("Reading file: %s\n", input_file);
    printf("Parameters: Omega=%.4f, T_max=%.4f, dt=%.4f\n", omega, t_max, dt);
    
    // Allocate 2D arrays for all observables
    double **p_bar = alloc_2d_array();
    double **size_gc = alloc_2d_array();
    double **err_size_gc = alloc_2d_array();
    double **avg_small = alloc_2d_array();
    double **err_avg_small = alloc_2d_array();
    
    // Read data from file
    FILE *fin = fopen(input_file, "r");
    if (!fin) {
        perror("Error opening input file");
        return 1;
    }
    
    char line[MAX_LINE];
    int data_count = 0;
    
    while (fgets(line, sizeof(line), fin)) {
        // Skip comment lines
        if (line[0] == '#') continue;
        
        double p1, p2, pb, sgc, esgc, asc, easc;
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", 
                   &p1, &p2, &pb, &sgc, &esgc, &asc, &easc) == 7) {
            int i = (int)round(p1 / GRID_STEP);
            int j = (int)round(p2 / GRID_STEP);
            
            if (i >= 0 && i < GRID_SIZE && j >= 0 && j < GRID_SIZE) {
                p_bar[i][j] = pb;
                size_gc[i][j] = sgc;
                err_size_gc[i][j] = esgc;
                avg_small[i][j] = asc;
                err_avg_small[i][j] = easc;
                data_count++;
            }
        }
    }
    fclose(fin);
    
    printf("Read %d data points\n", data_count);
    
    // Generate output filename
    char output_file[MAX_LINE];
    generate_output_filename(input_file, omega, t_max, dt, output_file);
    
    // Open output file
    FILE *fout = fopen(output_file, "w");
    if (!fout) {
        perror("Error opening output file");
        return 1;
    }
    
    // Write header
    fprintf(fout, "# Parametric extraction along curve x(t)=1-|cos(t)|, y(t)=1-|cos(Omega*t)|\n");
    fprintf(fout, "# Omega=%.6f, T_max=%.6f, dt=%.6f\n", omega, t_max, dt);
    fprintf(fout, "# Source file: %s\n", input_file);
    fprintf(fout, "# t p1(t) p2(t) p_bar(t) size_of_GC(t) err_on_size_of_GC(t) ");
    fprintf(fout, "average_small_cluster_size(t) err_on_average_small_cluster(t)\n");
    
    // Sample the curve and compute observables
    int n_samples = 0;
    for (double t = 0; t <= t_max; t += dt) {
        // Compute p1(t) = x(t) and p2(t) = y(t)
        double p1_t = 1.0 - fabs(cos(t));
        double p2_t = 1.0 - fabs(cos(omega * t));
        
        // Interpolate all observables at (p1(t), p2(t))
        Observables obs = interpolate_all(p_bar, size_gc, err_size_gc,
                                          avg_small, err_avg_small,
                                          p1_t, p2_t);
        
        fprintf(fout, "%.6f %.6f %.6f %.6f %.6e %.6e %.6f %.6f\n",
                t, p1_t, p2_t, obs.p_bar, obs.size_of_GC, obs.err_on_size_of_GC,
                obs.avg_small_cluster, obs.err_on_avg_small_cluster);
        
        n_samples++;
    }
    
    fclose(fout);
    
    // Free memory
    free_2d_array(p_bar);
    free_2d_array(size_gc);
    free_2d_array(err_size_gc);
    free_2d_array(avg_small);
    free_2d_array(err_avg_small);
    
    printf("Extracted %d samples along parametric curve\n", n_samples);
    printf("Results saved to: %s\n", output_file);
    
    return 0;
}
