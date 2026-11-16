%% Monte Carlo SSA Kinetic Parameter Optimization
% Developed by Eryk Grzechnik, Rafał Kopiasz and Dominik Jańczewski at
% Warsaw University of Technology.
% Description: This script performs global parameter optimization for macroscopic
% kinetic rate constants (k_a through k_e) of a polymer modification reaction
% assumed to be irreversible, bimolecular and second order.
% It utilizes a 1D lattice Stochastic Simulation Algorithm (SSA) executed by
% a compiled C++ MEX engine. The optimization relies on a single continuous
% Particle Swarm Optimization (PSO) to explore a bounded parameter space
% derived from an analytical approximation of the initial rate constant (k_a).

clc; clear; close all;

%% 1. USER SETTINGS & PHYSICAL CONSTANTS
% File Input/Output Configuration
target_file_pattern = {'data_template'};  % Substring to identify target experimental datasets
omit_file_pattern = {'files_to_exclude'}; % Substring to exclude specific files from the batch
reverse_order = 0;                        % Set to 1 to process files in reverse alphabetical order

% PSO Logging Criteria
max_log_count     = 25 * 1000;            % Hard stop: maximum allowed logs per mode to prevent infinite runs

% Physical System & Simulation Parameters
DP = 100;                                 % Degree of polymerization (lattice size per molecule)
N_molecules = 1000;                       % Number of polymer chains simulated per Monte Carlo evaluation
N_Avogadro = 6.022e23;
initial_concentration_M = 0.116;          % Initial concentration of the polymer in molarity (M)
time_to_seconds_multiplier = 1;           % Multiplier to convert raw experimental time units to seconds

% Model Complexity Selection
% Available variants
% '1p'      : a = b = c = d = e
% '2p'      : a = b = c = d > e
% '3pABC'   : a = b = c > d > e
% '3pBCD'   : a > b = c = d > e
% '3pAB_CD' : a = b > c = d > e
% '4pAB'    : a = b > c > d > e
% '4pBC'    : a > b = c > d > e
% '4pCD'    : a > b > c = d > e
% '5p'      : a > b > c > d > e  (Fully Unconstrained Baseline)

variants_to_test = {'5p'};

% Analytical Approximation Settings
prefit_cutoff_conv = 25;                  % Maximum conversion (%) used to approximate the analytical initial k_a
approx_up = 1e2;                          % Upper bound multiplier for the PSO search space (relative to k_a)
approx_down = 1e-5;                       % Lower bound multiplier for the PSO search space (relative to k_a)

% Particle Swarm Optimization (PSO) Settings
main_swarm_size = 500;                    % Number of particles evaluated concurrently per iteration
main_stall_limit = 100;                   % Iterations required without MSE improvement to trigger a stall termination
main_iter_limit = 1000;                   % Maximum theoretical PSO iterations
print_threshold = 10;                     % MSE threshold below which the engine prints particle updates to the console

fprintf('Starting Optimization at: %s\n', datetime('now','Format','HH:mm:ss'));

%% 2. DATASET IDENTIFICATION & QUEUEING
% Scans the active directory for Excel files matching the designated physical constraints.
all_excel = dir('*.xlsx');
file_names = string({all_excel.name});

valid_mask = contains(file_names, target_file_pattern) & ...
    ~contains(file_names, 'stat_test') & ...
    ~contains(file_names, omit_file_pattern);

selected_files = all_excel(valid_mask);

if reverse_order
    selected_files = selected_files(end:-1:1);
end

num_files = length(selected_files);
if num_files == 0
    error('No valid experimental .xlsx datasets found in the current folder.');
end
fprintf('Found %d valid datasets. Initializing Pipeline...\n', num_files);

%% 3. BATCH PROCESSING LOOP
tic;
for f = 1:num_files

    current_file = selected_files(f).name;
    [~, base_name, ~] = fileparts(current_file);

    fprintf('PROCESSING FILE %d/%d: %s\n', f, num_files, base_name);

    % Initialize Export Matrix
    vars_to_save = {'k_a', 'k_b', 'k_c', 'k_d', 'k_e', 'MSE', 'Threshold_Source'};
    stat_matrix = cell(length(vars_to_save) + 1, length(variants_to_test) + 1);
    stat_matrix{1, 1} = 'Parameter';
    stat_matrix(2:end, 1) = vars_to_save';
    stat_matrix(1, 2:end) = variants_to_test;

    % Preprocess Experimental Data
    exp_data = readmatrix(current_file);
    if size(exp_data, 2) >= 5 && ~isnan(exp_data(1,5))
        eq_reagent = exp_data(1,5);
    else
        eq_reagent = 2.0;
    end

    N_reagent = round(eq_reagent * DP * N_molecules);
    Vol_L     = (DP * N_molecules) / (N_Avogadro * initial_concentration_M);

    time_raw = exp_data(:,1);
    conv_raw = exp_data(:,2);
    nan_mask = isnan(time_raw) | isnan(conv_raw);

    time_exp_seconds = time_raw(~nan_mask) * time_to_seconds_multiplier;
    experimental_conversion = conv_raw(~nan_mask);

    if max(experimental_conversion) <= 1.01
        experimental_conversion = experimental_conversion * 100;
    end

    % Analytical k_a Approximation
    % Identifies the early-stage conversion points to mathematically approximate the unhindered rate constant
    valid_initial_idx = find(time_exp_seconds > 0 & experimental_conversion > 0);
    if isempty(valid_initial_idx)
        error('Invalid dataset %s: No data points found with t > 0 and conversion > 0.', base_name);
    end

    candidates = valid_initial_idx(experimental_conversion(valid_initial_idx) <= prefit_cutoff_conv);
    if ~isempty(candidates)
        [~, max_rel_idx] = max(experimental_conversion(candidates));
        idx_approx = candidates(max_rel_idx);
    else
        idx_approx = valid_initial_idx(1);
    end

    time_approx = time_exp_seconds(idx_approx);
    fractional_conversion = experimental_conversion(idx_approx) / 100;

    Ca_0 = initial_concentration_M;
    Cb_0 = initial_concentration_M * eq_reagent;
    Ca_t = Ca_0 * (1 - fractional_conversion);

    Cb_t_raw = Cb_0 - (Ca_0 * fractional_conversion);
    if Cb_t_raw < 1e-12
        warning('Reagent concentration dropped below zero during approximation. Applying 1e-12 minimum constraint.');
        Cb_t = 1e-12;
    else
        Cb_t = Cb_t_raw;
    end

    if abs(Cb_0 - Ca_0) > 1e-6
        ka_approx = (1 / (time_approx * (Cb_0 - Ca_0))) * log((Ca_0 * Cb_t) / (Cb_0 * Ca_t));
    else
        ka_approx = (1 / time_approx) * (1 / Ca_t - 1 / Ca_0);
    end

    % MODEL EVALUATION LOOP
    for m = 1:length(variants_to_test)
        current_mode = variants_to_test{m};
        fprintf('\n INITIATING %s MODEL EVALUATION \n', upper(current_mode));

        % Define uniform wide log-bounds for all 5 parameters based on the analytical approximation
        wide_lb = ka_approx * approx_down;
        wide_ub = ka_approx * approx_up;
        curr_lb = log10(repmat(wide_lb, 1, 5));
        curr_ub = log10(repmat(wide_ub, 1, 5));

        fprintf('  PSO Bounds | Lower: %.2e | Upper: %.2e\n', wide_lb, wide_ub);

        % Reset the persistent swarm logging function for the new evaluation
        dummy_output_fcn = @(optimValues, state) swarmLogger(optimValues, state, print_threshold, current_mode, max_log_count, main_swarm_size);
        dummy_output_fcn([], 'clear');

        % Configure the cost function to interface with the C++ SSA engine
        objective_func = @(log10_K_search) evaluate_mse(log10_K_search, current_mode, DP, ...
            N_molecules, N_reagent, time_exp_seconds, experimental_conversion, N_Avogadro, Vol_L);

        options = optimoptions('particleswarm', 'SwarmSize', main_swarm_size, ...
            'MaxIterations', main_iter_limit, 'StallIterLimit', main_stall_limit, ...
            'FunctionTolerance', 0, 'Display', 'none', ...
            'UseParallel', false, 'UseVectorized', true, 'OutputFcn', dummy_output_fcn);

        % Execute Particle Swarm Optimization
        seq_timer = tic;
        [log10_K_search, best_mse, ~, seq_output] = particleswarm(objective_func, 5, curr_lb, curr_ub, options);
        seq_elapsed = toc(seq_timer);

        % Process and Constrain Final Results
        best_K_macro = expand_and_constrain_matrix(log10_K_search, current_mode);

        if any(best_K_macro <= wide_lb * 1.05)
            warning('One or more optimal parameters hit the lower search boundary. Consider extending approx_down.');
        end

        actual_iters = max(1, seq_output.iterations);
        time_per_iter = seq_elapsed / actual_iters;

        fprintf(' Run Finished: %s | Duration: %.1fs | Speed: %.2fs/iter\n', datestr(now, 'HH:MM:SS'), seq_elapsed, time_per_iter);
        fprintf(' FINAL %s MSE: %.4f\n', upper(current_mode), best_mse);

        % Map optimal parameters to the statistical export matrix
        for col = 1:5
            stat_matrix{col + 1, m + 1} = best_K_macro(col);
        end
        stat_matrix{size(stat_matrix,1) - 1, m + 1} = best_mse;
        stat_matrix{size(stat_matrix,1), m + 1} = 'Pure PSO';

        % Extract the accumulated log and save the swarm history space to .mat
        current_log = dummy_output_fcn([], 'get');
        log_filename = sprintf('Swarm_History_%s_%s_%dmols_%s.mat', base_name, current_mode, N_molecules, datestr(now,'yyyymmdd_HHMMSS'));
        swarm_log = current_log;
        save(log_filename, 'swarm_log');
        fprintf(' Swarm history log saved as: %s (%d discrete parameter combinations)\n', log_filename, size(current_log, 1));

    end % End Modes Loop

    % Export final statistics to Excel
    base_filename = sprintf('%s_stat_test', base_name);
    output_file = sprintf('%s_1_%d_mols.xlsx', base_filename, N_molecules);
    file_counter = 1;
    while isfile(output_file)
        file_counter = file_counter + 1;
        output_file = sprintf('%s_%d_%d_mols.xlsx', base_filename, file_counter, N_molecules);
    end
    writecell(stat_matrix, output_file);
    fprintf('Saved optimal parameter matrix to %s.\n', output_file);

end % End Files Loop

total_time = toc;
fprintf('BATCH PROCESSING COMPLETE in %.2f minutes.\n', total_time/60);


%% HELPER FUNCTIONS

function mse_array = evaluate_mse(log10_K_search, mode, DP, Nm, Nr, t_eval, exp_data, Na, V)
% Evaluates the Mean Squared Error (MSE) of the kinetic rate constants by passing them to the C++ Monte Carlo Engine
HARD_PENALTY = 1e8;
PENALTY_MULTIPLIER = 1e14;

K_macro = expand_and_constrain_matrix(log10_K_search, mode);
k_micro_matrix = K_macro / (Na * V);

% CALLING THE EXECUTABLE C++ MEX FUNCTION
mse_array = kMC_PPM_engine(k_micro_matrix, DP, Nm, Nr, t_eval, exp_data);

% Apply strict monotonicity penalty to enforce hierarchy (k_a >= k_b >= k_c >= k_d >= k_e)
penalty = max(0, K_macro(:,2) - K_macro(:,1)) + ...
    max(0, K_macro(:,3) - K_macro(:,2)) + ...
    max(0, K_macro(:,4) - K_macro(:,3)) + ...
    max(0, K_macro(:,5) - K_macro(:,4));

invalid_mask = penalty > 0;
mse_array(invalid_mask) = HARD_PENALTY + (penalty(invalid_mask) * PENALTY_MULTIPLIER);
end

function stop = swarmLogger(optimValues, state, print_thresh, mode, max_log_count, swarm_size)
% Custom output function to aggressively index and log valid particle evaluations for post-optimization statistical analysis
PENALTY_THRESHOLD = 1e7;

persistent currentBestMSE;
persistent swarm_log;
persistent log_idx;

if strcmp(state, 'clear')
    currentBestMSE = Inf;
    swarm_log = NaN(max_log_count + (swarm_size * 2), 6); % 5 parameters + 1 MSE
    log_idx = 1;
    stop = false;
    return;
elseif strcmp(state, 'get')
    stop = swarm_log(1:log_idx-1, :); % Trim NaNs before returning
    return;
end

if isempty(currentBestMSE)
    currentBestMSE = Inf;
    swarm_log = NaN(max_log_count + (swarm_size * 2), 6);
    log_idx = 1;
end

stop = false;
current_mse = optimValues.bestfval;
current_iter = optimValues.iteration;

% Isolate and process valid evaluations (filters out heavily penalized non-monotonic vectors)
valid_mask = optimValues.swarmfvals < PENALTY_THRESHOLD;
if any(valid_mask)
    valid_swarm = optimValues.swarm(valid_mask, :);
    valid_mse = optimValues.swarmfvals(valid_mask);
    num_valid = size(valid_swarm, 1);

    linear_K = expand_and_constrain_matrix(valid_swarm, mode);

    % Dynamic overflow protection
    if (log_idx + num_valid - 1) > size(swarm_log, 1)
        swarm_log = [swarm_log; NaN(swarm_size * 2, 6)];
    end

    swarm_log(log_idx : log_idx + num_valid - 1, 1:5) = linear_K;
    swarm_log(log_idx : log_idx + num_valid - 1, 6) = valid_mse;
    log_idx = log_idx + num_valid;
end

if current_mse < currentBestMSE
    if current_mse <= print_thresh
        best_lin = expand_and_constrain_matrix(optimValues.bestx, mode);
        fprintf('   -> [Iter %d] MSE: %.4f | K: [%.2e, %.2e, %.2e, %.2e, %.2e]\n', ...
            current_iter, current_mse, best_lin(1), best_lin(2), best_lin(3), best_lin(4), best_lin(5));
    end
    currentBestMSE = current_mse;
end

% Terminate the optimization prematurely if sufficient statistical data is acquired
if (log_idx - 1) >= max_log_count
    stop = true;
end
end

function K_macro_matrix = expand_and_constrain_matrix(log10_K_search_matrix, mode)
% Processes an Nx5 matrix of log10 parameters, mapping them to linear space and applying hierarchical equality constraints
num_particles = size(log10_K_search_matrix, 1);
K_macro_matrix = 10.^log10_K_search_matrix;

for p = 1:num_particles
    K_macro_matrix(p, :) = apply_constraints(K_macro_matrix(p, :), mode);
end
end

function K_vec = apply_constraints(K_vec, mode)
% Applies mathematical constraints to force hierarchical parameter structures.
% Vector index mapping: 1=ka, 2=kb, 3=kc, 4=kd, 5=ke

switch mode

    % === BASELINE (5 PARAMETERS) ===
    case '5p'
        % Fully unconstrained aside from monotonicity penalty

        % === 4 PARAMETER MODELS (Baseline + 1 constraint) ===
    case '4pAB' % a=b
        K_vec(2) = K_vec(1);
    case '4pBC' % b=c
        K_vec(3) = K_vec(2);
    case '4pCD' % c=d
        K_vec(4) = K_vec(3);

        % === 3 PARAMETER MODELS (Baseline + 2 constraints) ===
    case '3pABC' % a=b=c
        K_vec(2:3) = K_vec(1);
    case '3pBCD' % b=c=d
        K_vec(3:4) = K_vec(2);
    case '3pAB_CD' % a=b AND c=d
        K_vec(2) = K_vec(1);
        K_vec(4) = K_vec(3);

        % === 2 PARAMETER MODELS ===
    case '2p' % a=b=c=d
        K_vec(2:4) = K_vec(1);

        % === EQUIVALENCE CONSTRAINT ===
    case '1p'
        K_vec(2:5) = K_vec(1);
end
end