% run tests

% unit tests
result_unit = runtests('Hes1_tissue_model_unit_tests')

% performance tests (save results after every part)
result_perf_seq = runperf('Hes1_tissue_model_performance_tests_seq');
path1 = './Results';
matrixname = ['result_perf_seq.mat'];
save(fullfile(path1,matrixname), 'result_perf_seq');

result_perf_par = runperf('Hes1_tissue_model_performance_tests_par');
matrixname = ['result_perf_par.mat'];
save(fullfile(path1,matrixname), 'result_perf_par');

result_perf_par_fix = runperf('Hes1_tissue_model_performance_tests_par_fix');
matrixname = ['result_perf_par_fix.mat'];
save(fullfile(path1,matrixname), 'result_perf_par_fix');

% performance summary
performance_summary_table = sampleSummary([result_perf_seq,result_perf_par,result_perf_par_fix])
matrixname = ['perf_summary.mat'];
save(fullfile(path1,matrixname), 'performance_summary_table');

% plot means results of performance test
iterations = [8; 189; 756; 1428];
seq_pos_std = performance_summary_table.Mean(1:4)+performance_summary_table.StandardDeviation(1:4);
seq_neg_std = performance_summary_table.Mean(1:4)-performance_summary_table.StandardDeviation(1:4);
par_pos_std = performance_summary_table.Mean(5:8)+performance_summary_table.StandardDeviation(5:8);
par_neg_std = performance_summary_table.Mean(5:8)-performance_summary_table.StandardDeviation(5:8);
par_fix_pos_std = performance_summary_table.Mean(9:12)+performance_summary_table.StandardDeviation(9:12);
par_fix_neg_std = performance_summary_table.Mean(9:12)-performance_summary_table.StandardDeviation(9:12);
figs = figure(1); clf,
x1 = [iterations; flipud(iterations)];
inBetween = [seq_neg_std; flipud(seq_pos_std)];
fill(x1, inBetween, [1,0.702,0.655], 'HandleVisibility','off');
hold on
plot(iterations, performance_summary_table.Mean(1:4), 'r-', 'LineWidth', 2, 'DisplayName', 'sequential')
hold on
plot(iterations, seq_pos_std, 'r--', 'LineWidth', 1.5, 'HandleVisibility','off')
hold on
plot(iterations, seq_neg_std, 'r--', 'LineWidth', 1.5, 'HandleVisibility','off')
hold on
x2 = [iterations; flipud(iterations)];
inBetween = [par_neg_std; flipud(par_pos_std)];
fill(x2, inBetween, [0,0.7,0.9], 'HandleVisibility','off');
hold on
plot(iterations, par_pos_std, 'Color', [0.196, 0.332, 0.482], 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility','off')
hold on
plot(iterations, par_neg_std, 'Color', [0.196, 0.332, 0.482], 'LineStyle', '--',  'LineWidth', 2, 'HandleVisibility','off')
hold on
plot(iterations, performance_summary_table.Mean(5:8), 'Color', [0.196, 0.332, 0.482], 'LineStyle', '-',  'LineWidth', 2, 'DisplayName', 'auto parallel')
hold on
x3 = [iterations; flipud(iterations)];
inBetween = [par_fix_neg_std; flipud(par_fix_pos_std)];
fill(x3, inBetween, [0,0.4,0.8], 'HandleVisibility','off');
hold on
plot(iterations, par_fix_pos_std, 'Color', [0, 0.2, 0.4], 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility','off')
hold on
plot(iterations, par_fix_neg_std, 'Color', [0, 0.2, 0.4], 'LineStyle', '--',  'LineWidth', 2, 'HandleVisibility','off')
hold on
plot(iterations, performance_summary_table.Mean(9:12), 'Color', [0, 0.2, 0.4], 'LineStyle', '-',  'LineWidth', 2, 'DisplayName', 'specified parallel')
ax = gca;
ax.FontSize = 13;
xlabel('Number of iterations', 'FontSize', 16)
ylabel('Mean performance (sec)', 'FontSize', 16)
legend('Location','northwest')

filename1 = ['performance_sequential_parallel_fix.png'];
saveas(figs,fullfile(path1,filename1))