% analyse results from parameter sweep

clear all;
close all;

% parameter ranges to use
D_d = [0.001,0.005,linspace(0.01,0.1,10)]; 
h = linspace(1,4,7);
gamma = linspace(1,9,17);

% define time span
T = 2000;
dt = 0.1;
Nt = T/dt;     % number of time points
t = linspace(0,T,Nt+1); % all time points

% parameter sweep using given parameter ranges with sequential solver
[mean_solns, peak_offsets, mean_period_length, total_combinations, is_parallel] = ...
          Hes1_tissue_model_sequential_solve(D_d, h, gamma);

% check for which parameters, oscillations have the required mean period
% length of 2-3 hours (i.e. 120-180 minutes)
possible_parameter_indices = find(mean_period_length(3,:)>120 & mean_period_length(3,:)<180);

possible_parameter_values = total_combinations(possible_parameter_indices,:)
