MODULAR CODE STRUCTURE - TWO-COLOUR BOND PERCOLATION WITH DISORDER PATTERN
==========================================================================

FILE STRUCTURE:
===============

Header Files (.h):
------------------
1. io.h            - Command-line parsing and file I/O operations
2. lattice.h       - Lattice structure and edge type assignment
3. percolation.h   - Percolation cluster analysis
4. utils.h         - Random number generation and statistics

Implementation Files (.c):
--------------------------
1. main.c          - Main simulation driver with disorder mode selection
2. io.c            - Input parsing, output file generation, header writing
3. lattice.c       - Square lattice creation and three disorder mode implementations
4. percolation.c   - Edge activation, cluster finding (Union-Find), statistics
5. utils.c         - RNG initialization, statistical calculations

Additional Tools:
-----------------
6. parametric_extract.c     - Standalone tool for parametric curve extraction
7. lattice_visualization.html - HTML visualization of lattice configurations

COMPILATION:
============

Main simulation:
    make                  # Build percolation_sim
    make clean            # Remove build artifacts
    make rebuild          # Clean and rebuild
    make run_example      # Run example simulation
    make test             # Run quick test

Manual compilation:
    gcc -O3 -Wall -Wextra -std=c99 -march=native -o percolation_sim \
        main.c lattice.c percolation.c utils.c io.c -lm

Parametric extraction tool:
    gcc -O3 -Wall -std=c99 -o parametric_extract parametric_extract.c -lm

USAGE:
======

Main simulation syntax:
    ./percolation_sim mode L M f p1_min p1_max p2_min p2_max delta_p seed

Parameters:
    mode    : Disorder mode - 'correlated', 'alternating', or 'uncorrelated'
    L       : Lattice size (L×L square lattice with periodic boundaries)
    M       : Number of Monte Carlo realizations
    f       : Fraction of type-1 edges (0-1, mapped to 0, 0.25, 0.5, 0.75, or 1.0)
    p1_min  : Minimum occupation probability for type-1 edges
    p1_max  : Maximum occupation probability for type-1 edges
    p2_min  : Minimum occupation probability for type-2 edges
    p2_max  : Maximum occupation probability for type-2 edges
    delta_p : Step size for probability sweep
    seed    : Random number generator seed (unsigned integer)

DISORDER MODES:
===============

The simulation supports three different ways to assign edge types:

1. CORRELATED MODE (mode='correlated'):
   Each node has exactly f×4 type-1 edges (hard constraint)

   Standard patterns for square lattice:
   - f = 0.00: All edges are type 2
   - f = 0.25: Checkerboard pattern (1 type-1, 3 type-2 per node)
   - f = 0.50: All horizontal edges type 1, all vertical type 2
   - f = 0.75: Three type-1 edges per node (horizontal + checkerboard vertical)
   - f = 1.00: All edges are type 1

   Properties:
   - Exact local constraints enforced at every node
   - Deterministic, spatially-structured patterns
   - No randomness in edge type assignment

2. ALTERNATING MODE (mode='alternating'):
   Edge types alternate along each lattice direction

   For f = 0.5:
   - Horizontal edges: type alternates along each row (random start per row)
   - Vertical edges: type alternates along each column (random start per column)
   - Creates striated patterns with random orientation

   Properties:
   - Perfect alternation along directions
   - Random starting types for each row/column
   - Works best for f = 0.5; other values give approximations
   - Warning: Odd L may cause PBC violations

3. UNCORRELATED MODE (mode='uncorrelated'):
   Each edge independently becomes type-1 with probability f

   Properties:
   - Complete spatial randomness
   - No local constraints
   - Edge types follow binomial distribution
   - Per-node type-1 count varies: mean = 4f, std ≈ 2√(f(1-f))

EXAMPLES:
=========

Correlated disorder with f=0.5 (horizontal type 1, vertical type 2):
    ./percolation_sim correlated 50 1000 0.5 0.0 1.0 0.0 1.0 0.05 12345

Alternating disorder with f=0.5 (striated pattern):
    ./percolation_sim alternating 100 500 0.5 0.0 1.0 0.0 1.0 0.05 42

Uncorrelated disorder with f=0.25 (random 25% type-1):
    ./percolation_sim uncorrelated 50 1000 0.25 0.0 1.0 0.0 1.0 0.1 999

Quick test (small parameters):
    ./percolation_sim alternating 10 10 0.5 0.5 0.5 0.5 0.5 0.1 42

Full phase diagram (101×101 grid):
    ./percolation_sim correlated 100 1000 0.5 0.0 1.0 0.0 1.0 0.01 12345

SIMULATION STRUCTURE:
====================

Initialization:
    Parse command-line arguments and validate parameters
    Map f to nearest standard value (0, 0.25, 0.5, 0.75, or 1.0)
    Initialize random number generator with provided seed
    Create square lattice with periodic boundary conditions

Edge Type Assignment (DONE ONCE):
    Assign edge types according to selected disorder mode
    - CORRELATED: Use deterministic patterns based on f
    - ALTERNATING: Use alternating patterns with random starts
    - UNCORRELATED: Randomly assign each edge with probability f
    Verify assignment (check constraints, count types, measure statistics)

Phase Diagram Calculation:
    FOR each p1 value from p1_min to p1_max (step delta_p):
        FOR each p2 value from p2_min to p2_max (step delta_p):
            Calculate p_bar = f × p1 + (1-f) × p2

            FOR each Monte Carlo realization (M):
                Activate edges stochastically:
                    - Type-1 edges active with probability p1
                    - Type-2 edges active with probability p2

                Find clusters using Union-Find algorithm

                Calculate cluster statistics:
                    - Giant component size (fraction of lattice)
                    - Small cluster susceptibility: χ = Σ(s²·n_s) / Σ(s·n_s)
                      where s = cluster size, n_s = number of clusters of size s
                      (excludes giant component)

                Store results for this realization

            Compute mean and standard error over M realizations
            Write data point to file

        Write blank line (for gnuplot pm3d format)

OUTPUT FORMAT:
==============

Filename format:
    phase_diagram_{mode}_L{L}_M{M}_f{f}.dat

    Examples:
    - phase_diagram_correlated_L100_M1000_f0.50.dat
    - phase_diagram_alternating_L50_M500_f0.25.dat
    - phase_diagram_uncorrelated_L100_M1000_f0.75.dat

File header:
    # Percolation simulation - {mode} disorder
    # L={L}, M={M}, f={f}, seed={seed}
    # p1 p2 p_bar size_of_GC err_on_size_of_GC average_small_cluster_size err_on_average_small_cluster

Output columns:
    1. p1                              - Type-1 edge occupation probability
    2. p2                              - Type-2 edge occupation probability
    3. p_bar                           - Effective probability: f×p1 + (1-f)×p2
    4. size_of_GC                      - Mean giant component size (fraction)
    5. err_on_size_of_GC               - Standard error on giant component size
    6. average_small_cluster_size      - Mean susceptibility χ = ⟨s²⟩/⟨s⟩
    7. err_on_average_small_cluster    - Standard error on susceptibility

Data format:
    - Space-separated columns
    - Blank lines separate rows (p1 values) for gnuplot pm3d
    - All statistics averaged over M Monte Carlo realizations



PARAMETRIC EXTRACTION TOOL:
============================

Purpose:
    Extract observables along parametric curves in the (p1, p2) phase diagram

    Specifically designed for curves of the form:
        p1(t) = 1 - |cos(t)|
        p2(t) = 1 - |cos(Ω·t)|

    where t ∈ [0, T_max] with step dt

Compilation:
    gcc -O3 -Wall -std=c99 -o parametric_extract parametric_extract.c -lm

Usage:
    ./parametric_extract <input_file> <Omega> <T_max> <delta_t>

Parameters:
    input_file : Phase diagram data file (output from main simulation)
    Omega      : Frequency ratio parameter (controls p2 oscillation)
    T_max      : Maximum time parameter
    delta_t    : Time step for sampling curve

Example:
    ./parametric_extract phase_diagram_alternating_L100_M1000_f0.50.dat 2.5 6.28 0.01

Method:
    - Reads phase diagram data into 101×101 grid
    - Samples parametric curve at discrete time points
    - Uses bilinear interpolation to estimate observables at each (p1(t), p2(t))
    - Outputs time series data

Output filename:
    parametric_omega{Ω}_tmax{T_max}_dt{dt}_{original_filename}

Output columns:
    1. t                               - Time parameter
    2. p1(t)                           - Type-1 probability at time t
    3. p2(t)                           - Type-2 probability at time t
    4. p_bar(t)                        - Effective probability
    5. size_of_GC(t)                   - Giant component size
    6. err_on_size_of_GC(t)            - Error on giant component
    7. average_small_cluster_size(t)   - Susceptibility
    8. err_on_average_small_cluster(t) - Error on susceptibility


HTML LATTICE VISUALIZATION:
===========================

HTML visualization:
    - Open lattice_visualization.html in web browser
    - Visualize lattice edge type assignments
    - Interactive display of different disorder modes
    - Useful for understanding spatial patterns

TROUBLESHOOTING:
================

Common issues:

1. "Invalid mode" error
   → Use 'correlated', 'alternating', or 'uncorrelated'

2. "L, M, and delta_p must be positive" error
   → Check that all parameters are positive values

3. Very long execution time
   → Reduce grid resolution (increase delta_p)
   → Reduce lattice size L
   → Reduce number of realizations M
   → Use smaller p ranges

4. Memory issues
   → Reduce lattice size L (memory scales as L²)
   → Check system RAM availability

5. Compilation errors
   → Ensure -lm flag for math library
   → Use C99 standard: -std=c99
   → Check that all source files are present

6. Unexpected f values in output
   → Input f is automatically mapped to nearest standard value
   → Check console output for "Mapped f value"

7. PBC violations warning (alternating mode, odd L)
   → Use even L for perfect alternation
   → Or accept slight imperfections with odd L

8. Parametric extraction: "Error opening input file"
   → Verify phase diagram file exists and path is correct
   → Check file permissions

FILE DEPENDENCIES:
==================

main.c requires:
    - io.h, lattice.h, percolation.h, utils.h

lattice.c requires:
    - lattice.h, utils.h

percolation.c requires:
    - percolation.h, lattice.h, utils.h

io.c requires:
    - io.h

utils.c requires:
    - utils.h

parametric_extract.c:
    - Standalone program (no dependencies on other modules)
