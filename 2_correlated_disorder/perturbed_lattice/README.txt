MODULAR CODE STRUCTURE - CORRELATED PERCOLATION SIMULATION
==========================================================

FILE STRUCTURE:
===============

Header Files (.h):
------------------
1. config.h          - Configuration structures, weight functions, and constants
2. union_find.h      - Union-Find data structure interface
3. random_utils.h    - Random number generation (Xoshiro256+, Gaussian)
4. lattice.h         - Lattice operations, edge management, and position perturbations
5. percolation.h     - Cluster finding and percolation statistics
6. output.h          - File output operations and filename generation

Implementation Files (.c):
--------------------------
1. main.c            - Main simulation driver with command-line argument parsing
2. union_find.c      - Union-Find implementation with path compression and union by rank
3. random_utils.c    - High-quality random number generators (Xoshiro256+, Box-Muller)
4. lattice.c         - Lattice creation, disorder generation, weight computation, and shuffling
5. percolation.c     - Cluster finding and susceptibility calculation (excluding largest cluster)
6. output.c          - Filename generation and results file writing

COMPILATION:
============

Using the provided Makefile:
    make              # Optimized build
    make debug        # Debug build with sanitizers
    make profile      # Profiling build
    make test         # Run quick tests
    make benchmark    # Run benchmarks
    make clean        # Clean build files

Manual compilation:
    gcc -O3 -march=native -ffast-math -o percolation main.c lattice.c output.c percolation.c random_utils.c union_find.c -lm

USAGE:
======

General syntax:
    ./percolation -L <size> -s <sigma> -w <weight_type> [options]

Required parameters:
    -L <int>    : Lattice size (L×L grid)
    -s <float>  : Position perturbation standard deviation (disorder strength)
    -w <type>   : Weight function type: 'power', 'exp', or 'bernoulli'

Weight function parameters:
    For 'power' (power law: ω(d) = d^α):
        -a <float> : Power law exponent α (default: 1.0)

    For 'exp' (exponential: ω(d) = Ω·exp(-d/λ)):
        -l <float> : Exponential decay length λ (required)
        -O <float> : Frequency prefactor Ω (required)

    For 'bernoulli' (threshold: ω(d) = θ(d-λ) + Ω·θ(λ-d)):
        -l <float> : Threshold distance λ (required)
        -O <float> : Frequency Ω for d < λ (required)

Optional parameters:
    -I <float>  : Initial simulation time T_INITIAL (default: 0.0)
    -F <float>  : Final simulation time T_FINAL (default: 10.0)
    -d <float>  : Time step dt (default: 0.01)
    -M <int>    : M1 - Number of activation realizations per time step (default: 100)
    -D <int>    : M2 - Number of disorder realizations (default: 100)
    -b          : Use open boundary conditions (default: periodic)
    -r <seed>   : Random seed for reproducibility (default: time-based)
    -h          : Show help message

EXAMPLES:
=========

Power law weight function:
    ./percolation -L 100 -s 0.1 -w power -a 2.0 -I 0.0 -F 10.0 -M 100 -D 100 -r 12345

Exponential weight function:
    ./percolation -L 100 -s 0.1 -w exp -l 1.5 -O 1.0 -I 0.0 -F 10.0 -M 100 -D 100 -r 12345

Bernoulli weight function:
    ./percolation -L 100 -s 0.1 -w bernoulli -l 1.0 -O 0.5 -I 0.0 -F 10.0 -M 100 -D 100 -r 12345

With open boundary conditions:
    ./percolation -L 100 -s 0.1 -w power -a 2.0 -b -r 12345

Quick test (small parameters):
    ./percolation -L 50 -s 0.1 -w power -a 2.0 -I 0.0 -F 5.0 -M 10 -D 10 -r 12345

Using default time range (0.0 to 10.0):
    ./percolation -L 100 -s 0.1 -w power -a 2.0 -r 12345


SIMULATION STRUCTURE:
====================

The simulation compares CORRELATED vs UNCORRELATED percolation:
- CORRELATED: Uses actual spatially-correlated edge weights from disorder realization
- UNCORRELATED: Uses same weights but shuffled randomly (destroys spatial correlations)

FOR each disorder realization (M2):
    Generate perturbed lattice positions (disorder strength σ)
    Compute edge weights based on perturbed distances
    Create CORRELATED version: original spatial arrangement
    Create UNCORRELATED version: shuffle weights randomly (once per disorder)

    FOR each time step from t=T_INITIAL to t=T_FINAL (with step dt):
        Compute activation probabilities p(t) = 1 - |cos(ω·t)| for all edges

        FOR each activation realization (M1):
            CORRELATED: Activate edges stochastically with correlated weights
            UNCORRELATED: Activate edges stochastically with shuffled weights

            Find clusters using Union-Find algorithm
            Compute statistics:
                - p(t): Fraction of active edges
                - P(t): Fraction of nodes in largest cluster
                - S(t): Susceptibility-related average cluster size
                        S = Σ(n_s·s²) / Σ(n_s·s), excluding largest cluster

            Accumulate for averaging over M1

        Average over M1 activation realizations

    Accumulate for averaging over M2

Average over M2 disorder realizations
Write results to file

PHYSICAL INTERPRETATION:
========================

Disorder Model:
    - Lattice sites have ideal positions (i,j)
    - Actual positions are perturbed: (i,j) → (i+δx, j+δy)
    - Perturbations δx, δy ~ N(0, σ²) are Gaussian random variables
    - Edge weights ω depend on actual perturbed distances

Weight Functions:
    - Power law: ω(d) = d^α
      Captures polynomial distance dependence

    - Exponential: ω(d) = Ω·exp(-d/λ)
      Models exponentially decaying interactions with prefactor Ω
      The prefactor Ω controls the overall frequency scale

    - Bernoulli: ω(d) = 1 if d≥λ, Ω if d<λ
      Binary threshold model with two distinct frequencies

Time Evolution:
    - Activation probability: p(t) = 1 - |cos(ω·t)|
    - Edges activate stochastically based on their frequency ω
    - System evolves through percolation transitions

Key Comparison:
    - CORRELATED: Preserves spatial correlations in disorder
    - UNCORRELATED: Same weight distribution but no spatial structure
    - Difference reveals the role of spatial correlations in percolation

OUTPUT FORMAT:
==============

Filename format:
    square_correlated_{weight}_L{L}_{BC}_sigma{σ}_[params]_Ti{T_INITIAL}_Tf{T_FINAL}_dt{dt}_M1{M1}_M2{M2}_seed{seed}.dat

    Where:
    - {weight}: power, exp, or bernoulli
    - {BC}: PBC (periodic) or OBC (open boundary conditions)
    - [params]: α for power, λ and Ω for exp, λ and Ω for bernoulli
    - Ti: Initial time
    - Tf: Final time

Output columns:
    1. t             - Time
    2. avg_prob      - Average activation probability ⟨1 - |cos(ω·t)|⟩

    Columns 3-5: CORRELATED (spatial correlations preserved)
    3. corr_p        - Fraction of active edges p(t)
    4. corr_P        - Fraction of nodes in largest cluster P(t)
    5. corr_S        - Average cluster size S(t) = Σ(n_s·s²)/Σ(n_s·s), excl. largest

    Columns 6-8: UNCORRELATED (shuffled weights, no spatial correlations)
    6. uncorr_p      - Fraction of active edges p(t)
    7. uncorr_P      - Fraction of nodes in largest cluster P(t)
    8. uncorr_S      - Average cluster size S(t), excl. largest

Averaging scheme:
    - Each observable is averaged over M1 activation realizations per time step
    - Then averaged over M2 disorder realizations
    - Total simulations per time point: 2 × M1 × M2
      (factor of 2 from correlated + uncorrelated versions)


TROUBLESHOOTING:
================

Common issues:

1. "Weight function type must be specified"
   → Use -w power, -w exp, or -w bernoulli

2. "Exponential weight requires lambda (-l) > 0"
   → Add -l parameter with positive value for exponential decay length

3. "Exponential weight requires Omega (-O) > 0"
   → Add -O parameter with positive value for frequency prefactor

4. "Bernoulli requires lambda (-l) > 0"
   → Add -l parameter with positive value for threshold distance

5. "Bernoulli requires Omega (-O) > 0"
   → Add -O parameter with positive value for low-distance frequency

6. "T_FINAL must be greater than T_INITIAL"
   → Ensure final time is larger than initial time

7. Compilation errors with AVX2
   → Remove -mavx2 flag from Makefile or update compiler

8. Very slow execution
   → Reduce L, M1, M2, or number of time steps ((T_FINAL-T_INITIAL)/dt)
   → Use make debug to check for issues

9. Memory allocation failures
   → Reduce lattice size L
   → System may have insufficient RAM for large L
