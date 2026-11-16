# Monte Carlo SSA Kinetic Parameter Optimization

## 1. Overview
This repository contains a computational algorithm designed to extract macroscopic kinetic rate constants (k_a, k_b, k_c, k_d, k_e) for polymer modification reactions assumed to be irreversible, bimolecular, and second order. 

The system utilizes a 1D lattice-based Stochastic Simulation Algorithm (SSA, or Gillespie algorithm) to model the step-by-step stochastic progression of functional group conversions along discrete polymer chains. This script employs a Particle Swarm Optimization (PSO) approach. The optimizer utilizes an analytically approximated boundary matrix based on the initial experimental reaction rate to dynamically bound the optimization range.

## 2. Physical Model Assumptions & Constraints
The physical architecture of the polymer simulation assumes:
- Reactive Sites: Every repeating unit (r.u.) is exactly one reactive site.
- Steric Hindrance Penalty: The local reactivity of a repeat unit is dictated strictly by the reaction states of its immediate neighbors within a window of ±3 repeating units. 
- State Mapping: Unhindered units are evaluated as rate `a`. Neighbor conversions at distances of 3, 2, 1, and adjacent on both sides map to sterically hindered rates `b`, `c`, `d`, and `e` respectively. 

By altering the `variants_to_test` variable in the script, users can enforce predefined equality constraints to simplify the optimization:
- '1p'      : a = b = c = d = e
- '2p'      : a = b = c = d > e
- '3pABC'   : a = b = c > d > e
- '3pBCD'   : a > b = c = d > e
- '3pAB_CD' : a = b > c = d > e
- '4pAB'    : a = b > c > d > e
- '4pBC'    : a > b = c > d > e
- '4pCD'    : a > b > c = d > e
- '5p'      : a > b > c > d > e (Fully Unconstrained Baseline)

## 3. System Requirements
To execute this script, the following software environment is required:
- MATLAB: Version 9.11.0.1837725 (R2021b) Update 2 (or newer).
- Optimization Toolbox (v9.2 was used).
- Global Optimization Toolbox (v4.6 was used): Required for the `particleswarm` function.
- C++ Compiler (MATLAB Support for MinGW-w64 C/C++ Compiler v21.2.0 was used): A compatible C++ compiler must be installed and configured within MATLAB (`mex -setup c++`) to interface with the external engine. MinGW-w64 or MSVC are recommended for Windows; GCC is recommended for Linux.

## 4. Input Data Format
The script requires experimental conversion data formatted as an Excel spreadsheet. A reference file, `data_template.xlsx`, is provided.
- Column 1: Time. The script assumes the data is in seconds. If needed, change the `time_to_seconds_multiplier` variable in the `.m` file to convert input data from hours/minutes to seconds.
- Column 2: Fractional or percentage conversion (0 to 1, or 0 to 100).
- Column 3: Error bars given as percent points (optional).
- Cell E1 (Row 1, Column 5): The molar equivalent of the modifying reagent used in the experiment (e.g., `2.15`). If left blank, the script defaults to 2.0.

Note: You must explicitly update the variables `DP`, `N_molecules`, and `initial_concentration_M` directly inside the `.m` script to match your specific experiment or modify it to read those values directly from input file.

## 5. Non-Reproducible Stochastic Randomization
By design, the C++ engine utilizes the `<chrono>` library alongside a bitwise XOR thread-index hash (`timeSeed ^ (p * 0x9E3779B9ULL)`) to explicitly decouple the Mersenne Twister seeds for every simulated particle. Because adjacent particles receive decorrelated seeds, execution paths are entirely non-reproducible.

To test the script, ensure all provided files are in the working directory and execute the script using `data_template.xlsx`. Using `DP = 100`, `N_molecules = 1000`, `max_log_count = 25000`, `main_swarm_size = 500`, and `variants_to_test = {'5p'}`. Depending on CPU core availability, the MSE should stochastically settle below 0.5, generating optimal constants near the following values (expressed in M^-1 s^-1):
 k_a : 2.7E-05
 k_b : 2.4E-05
 k_c : 1.6E-05
 k_d : 1.1E-05
 k_e : 1.1E-06

Upon completion, the engine exports an Excel statistical report and a `.mat` file (`Swarm_History...`). The `.mat` file contains the complete matrix of every discrete kinetic configuration evaluated, which is used for downstream post-analysis and parameter landscape mapping.

## 6. C++ Engine Compilation
The core computational workload is handled by the code shown in `kMC_PPM_engine.cpp` file. A pre-compiled `.mexw64` binary (for Windows) is provided. If you encounter architecture compatibility errors or are running on a different OS, you must compile the source file manually. 

Note that the `UseParallel` flag is disabled in the MATLAB script because the C++ engine internally manages multithreaded evaluation. To enable this OpenMP CPU multithreading for the massive swarm, execute the appropriate command in the MATLAB Command Window while in your working directory:

For Windows (MSVC Compiler):
```matlab
mex COMPFLAGS="$COMPFLAGS /openmp" kMC_PPM_engine.cpp

(This generates kMC_PPM_engine.mexw64)

For Linux (GCC Compiler):
```matlab
mex CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" kMC_PPM_engine.cpp

(This generates kMC_PPM_engine.mexa64)

## 7. C++ Engine Architecture & Logic
The provided kMC_PPM_engine.cpp is designed for computational throughput required by PSO operations.

O(1) Swap-and-Pop Memory: To avoid O(N) lattice scanning delays, the engine maintains five dynamically updating memory vectors (lists[0] through lists[4]). When a reaction occurs, the specific monomer is "popped" from its current state list and the last element of the list is swapped into its place. This allows the Gillespie loop to randomly select and execute reactions in absolute constant time.

Local Propensity Updates: Upon a successful reaction, the engine calculates structural state shifts solely for the immediate neighboring repeat units (within a radius of ±3 units).

Parallel Execution: The highest-level swarm evaluation loop is decorated with #pragma omp parallel for, forcing the C++ engine to distribute independent PSO particle evaluations concurrently across all available CPU cores.