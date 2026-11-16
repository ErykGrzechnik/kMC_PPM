#include "mex.h"
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // INPUT VALIDATION & ROBUSTNESS
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("kMC_PPM:inputs", "Exactly 6 inputs required (k_matrix, DP, Nm, Nr_initial, t_eval, exp_data).");
    }
    if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 5) {
        mexErrMsgIdAndTxt("kMC_PPM:inputs", "k_matrix must be a double array with exactly 5 columns.");
    }
    
    size_t t_length = mxGetNumberOfElements(prhs[4]);
    if (t_length != mxGetNumberOfElements(prhs[5])) {
        mexErrMsgIdAndTxt("kMC_PPM:inputs", "t_eval and exp_data arrays must have identical lengths.");
    }

    double* k_matrix = mxGetPr(prhs[0]);
    size_t num_particles = mxGetM(prhs[0]); 
    
    int DP = (int)mxGetScalar(prhs[1]);
    int Nm = (int)mxGetScalar(prhs[2]);
    double Nr_initial = mxGetScalar(prhs[3]);
    
    double* t_eval = mxGetPr(prhs[4]);
    double* exp_data = mxGetPr(prhs[5]);
    
    plhs[0] = mxCreateDoubleMatrix(num_particles, 1, mxREAL);
    double* mse_out = mxGetPr(plhs[0]);
    
    int N_total_sites = Nm * DP;

    #pragma omp parallel for
    for (int p = 0; p < num_particles; ++p) {
        try {
            double k[5];
            for (int i = 0; i < 5; ++i) {
                k[i] = k_matrix[p + i * num_particles];
            }
            
            // Strictly Stochastic Randomization using Bitwise XOR Decorrelation
            uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            std::mt19937 rng(timeSeed ^ (p * 0x9E3779B9ULL)); 
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            
            // O(1) SWAP-AND-POP MEMORY VECTORS
            std::vector<int> state(N_total_sites, 0); // 0=a, 1=b, 2=c, 3=d, 4=e, 5=reacted
            std::vector<int> pos(N_total_sites, 0);   // Tracks index position within the state list
            std::vector<int> lists[5];
            
            for (int i = 0; i < 5; ++i) {
                lists[i].reserve(N_total_sites);
            }
            
            // Initialize all sites as unhindered (State 0 / a)
            for (int i = 0; i < N_total_sites; ++i) {
                pos[i] = i;
                lists[0].push_back(i);
            }

            // Swap-and-Pop Lambda: Moves a monomer between state lists in O(1) time
            auto change_state = [&](int curr, int new_s) {
                int old_s = state[curr];
                if (old_s == new_s || old_s == 5) return; 
                
                // Pop from old list
                int p_idx = pos[curr];
                int last_monomer = lists[old_s].back();
                lists[old_s][p_idx] = last_monomer;
                pos[last_monomer] = p_idx;
                lists[old_s].pop_back();
                
                // Push to new list
                state[curr] = new_s;
                if (new_s != 5) {
                    pos[curr] = lists[new_s].size();
                    lists[new_s].push_back(curr);
                }
            };

            // Local Propensity Evaluator Lambda (Restricted to +/- 3 neighbors for 5 parameters)
            auto calc_state = [&](int curr) {
                if (state[curr] == 5) return 5;
                int j = curr % DP;
                
                bool L1 = (j > 0)     && (state[curr - 1] == 5);
                bool R1 = (j < DP - 1) && (state[curr + 1] == 5);
                bool L2 = (j > 1)     && (state[curr - 2] == 5);
                bool R2 = (j < DP - 2) && (state[curr + 2] == 5);
                bool L3 = (j > 2)     && (state[curr - 3] == 5);
                bool R3 = (j < DP - 3) && (state[curr + 3] == 5);
                
                if (L1 && R1) return 4;
                if (L1 || R1) return 3;
                if (L2 || R2) return 2;
                if (L3 || R3) return 1;
                return 0;
            };
            
            int N_modif = 0;
            double t = 0.0;
            double N_reag = Nr_initial;
            
            std::vector<double> t_events;
            std::vector<double> conv_events;
            t_events.reserve(N_total_sites + 1);
            conv_events.reserve(N_total_sites + 1);
            t_events.push_back(0.0);
            conv_events.push_back(0.0);
            
            // GILLESPIE SSA LOOP
            while (N_reag > 0 && N_modif < N_total_sites) {
                
                double a_rates[5];
                double a0 = 0.0;
                
                for (int s = 0; s < 5; ++s) {
                    a_rates[s] = k[s] * lists[s].size() * N_reag;
                    a0 += a_rates[s];
                }
                
                if (a0 <= 0.0) break;
                
                double u1 = unif(rng);
                if(u1 < 1e-15) u1 = 1e-15; 
                t -= std::log(u1) / a0;
                
                double r = unif(rng) * a0;
                int pick = -1;
                double running_sum = 0.0;
                
                // Roulette wheel channel selection
                for (int s = 0; s < 5; ++s) {
                    running_sum += a_rates[s];
                    if (r < running_sum && lists[s].size() > 0) {
                        std::uniform_int_distribution<int> dist(0, lists[s].size() - 1);
                        pick = lists[s][dist(rng)];
                        break;
                    }
                }
                
                // Execute Reaction & Local Updates
                if (pick != -1) {
                    change_state(pick, 5);
                    
                    // Update ONLY the +/- 3 local neighbors (reduced from 4 for 5-parameter limit)
                    int m = pick / DP;
                    int j = pick % DP;
                    for (int d = -3; d <= 3; ++d) {
                        if (d == 0) continue;
                        int nj = j + d;
                        if (nj >= 0 && nj < DP) {
                            int ncurr = m * DP + nj;
                            if (state[ncurr] != 5) {
                                int new_s = calc_state(ncurr);
                                change_state(ncurr, new_s);
                            }
                        }
                    }
                    
                    N_modif++;
                    N_reag--;
                    t_events.push_back(t);
                    conv_events.push_back(((double)N_modif / (double)N_total_sites) * 100.0);
                }
            }
            
            std::vector<double> sim_conv(t_length, 0.0);
            int current_event_idx = 0;
            for (size_t eval_idx = 0; eval_idx < t_length; ++eval_idx) {
                double current_t = t_eval[eval_idx];
                while (current_event_idx + 1 < t_events.size() && t_events[current_event_idx + 1] <= current_t) {
                    current_event_idx++;
                }
                sim_conv[eval_idx] = conv_events[current_event_idx];
            }
            
            double particle_mse = 0.0;
            for (size_t t_i = 0; t_i < t_length; ++t_i) {
                double diff = sim_conv[t_i] - exp_data[t_i];
                particle_mse += (diff * diff);
            }
            
            mse_out[p] = particle_mse / t_length;
        }
        catch(...) {
            // Failsafe catch block to prevent OpenMP memory allocation errors from hard-crashing MATLAB
            mse_out[p] = 1e8; 
        }
    }
}