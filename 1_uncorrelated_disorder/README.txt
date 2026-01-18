MODULAR CODE STRUCTURE - UNCORRELATED PERCOLATION SIMULATION
============================================================

FILE STRUCTURE:
===============

Header Files (.h):
------------------
1. config.h          - Configuration structures and constants
2. union_find.h      - Union-Find data structure interface
3. random_utils.h    - Random number generation functions
4. lattice.h         - Lattice operations and edge management
5. percolation.h     - Cluster finding and statistics
6. output.h          - File output operations

Implementation Files (.c):
--------------------------
1. main.c            - Main simulation driver
2. union_find.c      - Union-Find implementation
3. random_utils.c    - Random number generators (Gaussian, Bernoulli)
4. lattice.c         - Lattice creation and weight generation
5. percolation.c     - Cluster finding and susceptibility calculation
6. output.c          - Filename generation and results writing

COMPILATION:
============

Using the provided Makefile:
    make

Or manually:
    gcc -Wall -O3 -std=c99 -o square_uncorrelated main.c union_find.c random_utils.c lattice.c percolation.c output.c -lm

USAGE:
======

Gaussian distribution:
    ./square_uncorrelated gaussian L mu sigma M D T_INITIAL T_FINAL dt seed

Bernoulli distribution:
    ./square_uncorrelated bernoulli L p Omega M D T_INITIAL T_FINAL dt seed

Parameters:
    L        : Lattice size (L×L)
    mu       : Mean frequency (Gaussian mode)
    sigma    : Frequency standard deviation (Gaussian mode)
    p        : Probability of frequency = 1 (Bernoulli mode)
    Omega    : High frequency value with probability 1-p (Bernoulli mode)
    M        : Number of disorder realizations
    D        : Number of percolation processes per time step
    T_INITIAL: Initial simulation time
    T_FINAL  : Final simulation time
    dt       : Time step
    seed     : Random seed for reproducibility

EXAMPLES:
=========

Gaussian:
    ./square_uncorrelated gaussian 100 1.0 0.5 10 100 0.0 10.0 0.1 12345

Bernoulli:
    ./square_uncorrelated bernoulli 100 0.5 2.0 10 100 0.0 10.0 0.1 12345


SIMULATION STRUCTURE:
====================

FOR each disorder realization (M):
    Generate random edge weights (Gaussian or Bernoulli)
    
    FOR each time step from T_INITIAL to T_FINAL:
        Calculate average activation probability
        
        FOR each percolation process (D):
            Activate edges stochastically
            Find clusters using Union-Find
            Compute statistics (active fraction, largest cluster, susceptibility)
            Accumulate for averaging

Compute final statistics (means and standard deviations)
Write results to file

OUTPUT FORMAT:
==============

Columns:
1. t             - Time
2. avg_prob      - Average activation probability (over M)
3. active_frac   - Mean fraction of active edges (over M×D)
4. active_std    - Standard deviation of active fraction
5. largest_frac  - Mean largest cluster fraction (over M×D)
6. largest_std   - Standard deviation of largest cluster fraction
7. suscept_avg   - Mean susceptibility: sum(s²·nₛ)/sum(s·nₛ) (over M×D)
8. suscept_std   - Standard deviation of susceptibility

