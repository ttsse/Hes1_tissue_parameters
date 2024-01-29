% test parameter values

original_values = [0.003; 1; 8]; % from masters thesis

% minimum parameter sweep with two values each
min_values = [linspace(0.001, 0.1, 2); linspace(1,9,2); linspace(1,9,2)];

% zero parameter values
zero_D_d = [0; 1; 8];
zero_h = [0.003; 0; 8];
zero_gamma = [0.003; 1; 0];
all_zero = [0; 0; 0];

% negative parameter values
neg_D_d = [-0.003; 1; 8];
neg_h = [0.003; -1; 8];
neg_gamma = [0.003; 1; -8];
all_neg = [-0.003; -1; -8];

% parameter values empty
no_D_d = [1; 8];
no_h = [0.003; 8];
no_gamma = [0.003; 1];

% tolerance for comparison of results
tol = 1e-6;

%% Test 1 Require Positive Parameter Inputs

% sequential solver

% check that right errors are output for zero values
assertError(@() Hes1_tissue_model_sequential_solve(zero_D_d(1,:), zero_D_d(2,:), zero_D_d(3,:)), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(zero_h(1,:), zero_h(2,:), zero_h(3,:)), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(zero_gamma(1,:), zero_gamma(2,:), zero_gamma(3,:)), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(all_zero(1,:), all_zero(2,:), all_zero(3,:)), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

% check that right errors are output for negative values
assertError(@() Hes1_tissue_model_sequential_solve(neg_D_d(1,:), neg_D_d(2,:), neg_D_d(3,:)), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(neg_h(1,:), neg_h(2,:), neg_h(3,:)), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(neg_gamma(1,:), neg_gamma(2,:), neg_gamma(3,:)), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_sequential_solve(all_neg(1,:), all_neg(2,:), all_neg(3,:)), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

% parallel solver with auto graining

% check that right errors are output for zero values
assertError(@() Hes1_tissue_model_parallel_solve(zero_D_d(1,:), zero_D_d(2,:), zero_D_d(3,:),0), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(zero_h(1,:), zero_h(2,:), zero_h(3,:),0), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(zero_gamma(1,:), zero_gamma(2,:), zero_gamma(3,:),0), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(all_zero(1,:), all_zero(2,:), all_zero(3,:),0), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

% check that right errors are output for negative values
assertError(@() Hes1_tissue_model_parallel_solve(neg_D_d(1,:), neg_D_d(2,:), neg_D_d(3,:),0), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(neg_h(1,:), neg_h(2,:), neg_h(3,:),0), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(neg_gamma(1,:), neg_gamma(2,:), neg_gamma(3,:),0), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(all_neg(1,:), all_neg(2,:), all_neg(3,:),0), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

% parallel solver with specified graining

% check that right errors are output for zero values
assertError(@() Hes1_tissue_model_parallel_solve(zero_D_d(1,:), zero_D_d(2,:), zero_D_d(3,:),1), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(zero_h(1,:), zero_h(2,:), zero_h(3,:),1), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(zero_gamma(1,:), zero_gamma(2,:), zero_gamma(3,:),1), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(all_zero(1,:), all_zero(2,:), all_zero(3,:),1), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

% check that right errors are output for negative values
assertError(@() Hes1_tissue_model_parallel_solve(neg_D_d(1,:), neg_D_d(2,:), neg_D_d(3,:),1), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(neg_h(1,:), neg_h(2,:), neg_h(3,:),1), ...
    'MathBiology:negativeParameters','h needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(neg_gamma(1,:), neg_gamma(2,:), neg_gamma(3,:),1), ...
    'MathBiology:negativeParameters','gamma needs to be positive.');

assertError(@() Hes1_tissue_model_parallel_solve(all_neg(1,:), all_neg(2,:), all_neg(3,:),1), ...
    'MathBiology:negativeParameters','D_d needs to be positive.');

%% Test 2 Positive Number of Iterations

% sequential solver

% original values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_sequential_solve(original_values(1,:), original_values(2,:), original_values(3,:));
assert(~isempty(total_combinations))

% minimum parameter sweep values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_sequential_solve(min_values(1,:), min_values(2,:), min_values(3,:));
assert(~isempty(total_combinations))

% check that right errors are output for missing values
assertError(@() Hes1_tissue_model_parallel_solve([], no_D_d(1,:), no_D_d(2,:)), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_h(1,:), [], no_h(2,:)), ...
    'MATLAB:notEnoughInputs', 'h cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_gamma(1,:), no_gamma(2,:), []), ...
    'MATLAB:notEnoughInputs', 'gamma cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve([], [], []), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');


% parallel solver with auto graining

% original values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 0);
assert(~isempty(total_combinations))

% minimum parameter sweep values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 0);
assert(~isempty(total_combinations))

% check that right errors are output for missing values
assertError(@() Hes1_tissue_model_parallel_solve([], no_D_d(1,:), no_D_d(2,:), 0), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_h(1,:), [], no_h(2,:), 0), ...
    'MATLAB:notEnoughInputs', 'h cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_gamma(1,:), no_gamma(2,:), [], 0), ...
    'MATLAB:notEnoughInputs', 'gamma cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve([], [], [], 0), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');

% parallel solver with specified graining

% original values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 1);
assert(~isempty(total_combinations))

% minimum parameter sweep values
[~, ~, ~, total_combinations, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 1);
assert(~isempty(total_combinations))

% check that right errors are output for missing values
assertError(@() Hes1_tissue_model_parallel_solve([], no_D_d(1,:), no_D_d(2,:), 1), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_h(1,:), [], no_h(2,:), 1), ...
    'MATLAB:notEnoughInputs', 'h cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve(no_gamma(1,:), no_gamma(2,:), [], 1), ...
    'MATLAB:notEnoughInputs', 'gamma cannot be empty.');

assertError(@() Hes1_tissue_model_parallel_solve([], [], [], 1), ...
    'MATLAB:notEnoughInputs', 'D_d cannot be empty.');

%% Test 3 Non-Negative Mean Results

% sequential solver

% original values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_sequential_solve(original_values(1,:), original_values(2,:), original_values(3,:));
assert(all(mean_solns(:)))

% minimum parameter sweep values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_sequential_solve(min_values(1,:), min_values(2,:), min_values(3,:));
assert(all(mean_solns(:)))

% parallel solver with auto graining

% original values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 0);
assert(all(mean_solns(:)))

% minimum parameter sweep values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 0);
assert(all(mean_solns(:)))

% parallel solver with specified graining

% original values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 1);
assert(all(mean_solns(:)))

% minimum parameter sweep values
[mean_solns, ~, ~, ~, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 1);
assert(all(mean_solns(:)))

%% Test 4 Check Parallel Usage

% check that sequential code runs sequentially (i.e. not using parallel
% workers)

% original values
[~, ~, ~, ~, is_parallel_seq] = Hes1_tissue_model_sequential_solve(original_values(1,:), original_values(2,:), original_values(3,:));
assert(all(is_parallel_seq(:) == 0))

% minimum parameter sweep values
[~, ~, ~, ~, is_parallel_seq] = Hes1_tissue_model_sequential_solve(min_values(1,:), min_values(2,:), min_values(3,:));
assert(all(is_parallel_seq(:) == 0))

% check that parallel code (auto-grained) runs in parallel (i.e. using parallel workers)

% original values
[~, ~, ~, ~, is_parallel_par] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 0);
assert(all(is_parallel_par(:) == 1))

% minimum parameter sweep values
[~, ~, ~, ~, is_parallel_par] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 0);
assert(all(is_parallel_par(:) == 1))

% check that parallel code (specified graining) runs in parallel (i.e. using parallel workers)

% original values
[~, ~, ~, ~, is_parallel_par] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 1);
assert(all(is_parallel_par(:) == 1))

% minimum parameter sweep values
[~, ~, ~, ~, is_parallel_par] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 1);
assert(all(is_parallel_par(:) == 1))

%% Test 5 Sequential-Parallel Comparison

% original values

% run sequential and parallel runs (both with auto and specified graining)
[mean_solns_seq, peak_offsets_seq, mean_period_length_seq, total_combinations_seq, ~] = Hes1_tissue_model_sequential_solve(original_values(1,:), original_values(2,:), original_values(3,:));
[mean_solns_par, peak_offsets_par, mean_period_length_par, total_combinations_par, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 0);
[mean_solns_par_fix, peak_offsets_par_fix, mean_period_length_par_fix, total_combinations_par_fix, ~] = Hes1_tissue_model_parallel_solve(original_values(1,:), original_values(2,:), original_values(3,:), 1);

% difference between sequential and parallel runs
abs_mean_solns = abs(mean_solns_seq-mean_solns_par);
abs_mean_solns_fix = abs(mean_solns_seq-mean_solns_par_fix);
abs_peak_offsets = abs(peak_offsets_seq-peak_offsets_par);
abs_peak_offsets_fix = abs(peak_offsets_seq-peak_offsets_par_fix);
abs_mean_period_length = abs(mean_period_length_seq-mean_period_length_par);
abs_mean_period_length_fix = abs(mean_period_length_seq-mean_period_length_par_fix);
abs_total_combinations = abs(total_combinations_seq-total_combinations_par);
abs_total_combinations_fix = abs(total_combinations_seq-total_combinations_par_fix);

% check that results for sequential and parallel runs are equal (up to a
% tolerance)
assert(all(abs_mean_solns(:) < tol))
assert(all(abs_mean_solns_fix(:) < tol))
assert(all(abs_peak_offsets(:) < tol))
assert(all(abs_peak_offsets_fix(:) < tol))
assert(all(abs_mean_period_length(:) < tol))
assert(all(abs_mean_period_length_fix(:) < tol))
assert(all(abs_total_combinations(:) < tol))
assert(all(abs_total_combinations_fix(:) < tol))

% minimum parameter sweep values

% run sequential and parallel runs
[mean_solns_seq, peak_offsets_seq, mean_period_length_seq, total_combinations_seq, ~] = Hes1_tissue_model_sequential_solve(min_values(1,:), min_values(2,:), min_values(3,:));
[mean_solns_par, peak_offsets_par, mean_period_length_par, total_combinations_par, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 0);
[mean_solns_par_fix, peak_offsets_par_fix, mean_period_length_par_fix, total_combinations_par_fix, ~] = Hes1_tissue_model_parallel_solve(min_values(1,:), min_values(2,:), min_values(3,:), 1);

% difference between sequential and parallel runs
abs_mean_solns = abs(mean_solns_seq-mean_solns_par);
abs_mean_solns_fix = abs(mean_solns_seq-mean_solns_par_fix);
abs_peak_offsets = abs(peak_offsets_seq-peak_offsets_par);
abs_peak_offsets_fix = abs(peak_offsets_seq-peak_offsets_par_fix);
abs_mean_period_length = abs(mean_period_length_seq-mean_period_length_par);
abs_mean_period_length_fix = abs(mean_period_length_seq-mean_period_length_par_fix);
abs_total_combinations = abs(total_combinations_seq-total_combinations_par);
abs_total_combinations_fix = abs(total_combinations_seq-total_combinations_par_fix);

% check that results for sequential and parallel runs are equal (up to a
% tolerance)
assert(all(abs_mean_solns(:) < tol))
assert(all(abs_mean_solns_fix(:) < tol))
assert(all(abs_peak_offsets(:) < tol))
assert(all(abs_peak_offsets_fix(:) < tol))
assert(all(abs_mean_period_length(:) < tol))
assert(all(abs_mean_period_length_fix(:) < tol))
assert(all(abs_total_combinations(:) < tol))
assert(all(abs_total_combinations_fix(:) < tol))
