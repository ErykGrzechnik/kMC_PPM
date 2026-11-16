# kMC_PPM
Repository of the kinetic Monte Carlo algorithm used to fit reaction rate constants to the post-polymerization modification reactions of linear polymers. 

Monte Carlo SSA Kinetic Parameter Optimization

1. Overview

This repository contains a computational algorithm designed to extract macroscopic kinetic rate constants (k_a, k_{b0}, k_b, k_c, k_d, k_e) for polymer modification reactions assumed to be irreversible, bimolecular and second order. 
The system utilizes a 1D lattice-based Stochastic Simulation Algorithm (SSA, or Gillespie algorithm) to model the step-by-step stochastic progression of functional group conversions along discrete polymer chains.
This script employs a Massive Swarm Particle Swarm Optimization (PSO) approach. The optimizer utilizes an analytically approximated boundary matrix based on initial experimental reaction rate to narrow the optimization range.

2. System Requirements

To execute this script, the following software environment is required:

- MATLAB: Version 9.11.0.1837725 (R2021b) Update 2 (or newer).

- MATLAB Add-Ons:

  * Global Optimization Toolbox: Required for the `particleswarm` function.

  * MATLAB Coder / Custom C++ Compiler setup: Required to interface with the external C++ engine.

  * C++ Compiler: A compatible C++ compiler must be installed and configured within MATLAB (e.g., MinGW-w64 or MSVC for Windows, GCC for Linux).

3. Input Data Format

The script requires experimental conversion data formatted as an Excel spreadsheet. A reference file, `data_template.xlsx`, is provided.

Column 1: Time (the script assumes the data is in seconds. If needed, change the time_to_seconds_multiplier variable to convert input data from hours to seconds etc. ).

Column 2: Fractional or percentage conversion (0 to 1, or 0 to 100).

Column 3: Error bars given as percent points (optional).

Cell E1 (Row 1, Column 5): The molar equivalent of the modifying reagent used in the experiment (e.g., `2.15`). If left blank, the script defaults to 2.0.

4. To test the script put all of the provided files in the working directory. Using the 'data_template' as input with min_log_count = 20000, max_log_count = 25000, DP = 100, N_molecules = 1000, and variants_to_test = {'5p'}
the optimization should finish in about 5 minutes with final parameters close to the following:
k_a	2.70E-05
k_{b0}	2.70E-05
k_b	2.35E-05
k_c	1.70E-05
k_d	1.07E-05
k_e	1.11E-06
MSE	0.4684

Those reaction rate constants are expressed in (M^-1 * s^-1)

5. C++ Engine Compilation

The core computational workload is directed to C++ engine (`kMC_PPM_engine.cpp`). A pre-compiled `.mexw64` binary (for Windows) is provided. 
However, if you encounter architecture compatibility errors, or if you are running on a different operating system, you must compile the source file manually.

To compile the engine and enable CPU multithreading (OpenMP) for the massive swarm evaluations, execute the appropriate command in the MATLAB Command Window while in your working directory:

For Windows (MSVC Compiler):
mex COMPFLAGS="$COMPFLAGS /openmp" kMC_PPM_engine.cpp

(This generates `kMC_PPM_engine.mexw64`)

For Linux (GCC Compiler):
mex CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" kMC_PPM_engine.cpp

(This generates `kMC_PPM_engine.mexa64`)

6. C++ Engine Architecture & Logic

The provided `kMC_PPM_engine.cpp` file is designed for high computational throughput required by PSO operations.

* Initialization & Randomization: Each particle evaluates a discrete matrix of reaction rate constants. 
To avoid pseudo-random resonance, the `<chrono>` library dynamically seeds a Mersenne Twister (`std::mt19937`) unique to every simulated parameter set.

* O(1) Swap-and-Pop Memory: To avoid O(N) lattice scanning delays, the engine maintains six dynamically updating memory vectors (`lists[0]` through `lists[5]`),
 mapping directly to unreacted monomer states (unhindered a, or sterically hindered b0, b, c, d, e). When a reaction occurs,
 the specific monomer is "popped" from its current state list and the last element of the list is swapped into its place. 
This allows the Gillespie loop to randomly select and execute reactions in absolute O(1) constant time.

* Local Neighborhood Updates: Upon a successful reaction, the engine calculates structural state shifts for only the immediate neighboring repeat units.

* Parallel Execution: The highest-level swarm evaluation loop is decorated with `#pragma omp parallel for`,
 forcing the C++ engine to distribute particle evaluations across all available CPU cores simultaneously,
 resulting in near-linear scaling with processor thread count.
