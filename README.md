# Dynamical_Entanglement_Percolation
This repositories contains the codes used to perform the numerical simulations in the paper "Dynamical entanglement percolation iwith spatially correlated disorder"

Each subfolder contains a README file, and each program is "self-explanatory", in the sense that the user is asked for the correct parameters in input when executing the programs.

The folder "1_uncorrelated_disorder" considers uncorrelated random frequencies with either gaussian or Bernoulli distribution.

The folder "2_correlated_disorder" contains the codes for the perturbed square lattice in which frequencies depend on the eucliden distances, and we assume a functional dependence of omega on d of three possible types: (i) exponential, (ii) algebraic, (iii) threshold-like. Model (iii) is the one used to produce the results in the main.

The folder "3_two_colour_model" contains the codes to simulate the two-colour bond percolation model, running Montecarlo simulations for the square lattice and the numerical solution of the recursive equations for the mean-field solution.
