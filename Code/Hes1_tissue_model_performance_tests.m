% test parameter values

original_values = [0.003; 1; 8]; % from masters thesis

% parameter sweeps with different numbers of parameters

% 8 total iterations (combinations of all parameter values)
min_values = [linspace(0.001, 0.1, 2); linspace(1,9,2); linspace(1,9,2)];

% 189 total iterations
mid_D_d = [0.001, 0.01, 0.1];
mid_h = linspace(1,4,7); 
mid_gamma = linspace(5,9,9);

% 1428 total iterations
max_D_d = [0.001,0.005,linspace(0.01,0.1,10)]; 
max_h = linspace(1,4,7);
max_gamma = linspace(1,9,17);

%% Test 1 Sequential 8 Iterations

Hes1_tissue_model_sequential_solve(min_values(1,:), min_values(2,:), min_values(3,:));

%% Test 2 Sequential 189 Iterations

Hes1_tissue_model_sequential_solve(mid_D_d, mid_h, mid_gamma);

%% Test 3 Sequential 756 Iterations

Hes1_tissue_model_sequential_solve(max_D_d, mid_h, mid_gamma);

%% Test 4 Sequential 1428 Iterations

Hes1_tissue_model_sequential_solve(max_D_d, max_h, max_gamma);

%% Test 5 Parallel 8 Iterations

Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:));

%% Test 6 Parallel 189 Iterations

Hes1_tissue_model_parallel_solve(mid_D_d, mid_h, mid_gamma);

%% Test 7 Parallel 756 Iterations

Hes1_tissue_model_parallel_solve(max_D_d, mid_h, mid_gamma);

%% Test 8 Parallel 1428 Iterations

Hes1_tissue_model_parallel_solve(max_D_d, max_h, max_gamma);
