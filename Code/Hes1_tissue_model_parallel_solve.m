function [mean_solns, peak_offsets, mean_period_length, total_combinations, is_parallel] = ...
          Hes1_tissue_model_parallel_solve(D_d, h, gamma,fixed)

tic
close all;
format long;

% check that input is in the right format
% first need to be sure there are enough inputs
if nargin < 3
    % use your suitable error id here
    error('MyError:InadequateInput','Not enough inputs');      
end

% then check if the inputs have empty entries
%  before checking if they are nonsense
if isempty(D_d)
    error('MATLAB:notEnoughInputs', 'D_d cannot be empty.');
elseif ~isnumeric(D_d) % also need to make sure they have numeric entries
    % use your suitable error id here
    error('MyError:nonNumericInputs', 'D_d entries must be numeric.');
elseif D_d<=0
    error('MathBiology:negativeParameters', 'D_d needs to be positive.');
end

if isempty(h)
    error('MATLAB:notEnoughInputs', 'h cannot be empty.');
elseif ~isnumeric(h)
    % use your suitable error id here
    error('MyError:nonNumericInputs', 'h must be numeric.');
elseif h<=0
    error('MathBiology:negativeParameters', 'h needs to be positive.');
end 

if isempty(gamma)
    error('MATLAB:notEnoughInputs', 'gamma cannot be empty');
elseif ~isnumeric(gamma)
    % use your suitable error id here
    error('MyError:nonNumericInputs', 'gamma must be numeric.');
elseif gamma<=0
    error('MathBiology:negativeParameters', 'gamma needs to be positive.');
end 

% define parameters of the system
alpha_d = 0.05;
alpha_m = 0.2;
alpha_p = 0.05;
alpha_n = 0.1;
mu_d = log(2)/50;
mu_m = log(2)/24.1;
mu_p = log(2)/22.3;
mu_n = log(2)/22;

% define parameters for numerical solution of x and t
a = 0;
b = 5.12; % full size of fetus (in mm) at E10.5
T = 2000;
dx = 0.1;
dt = 0.1;

Nx = (b-a)/dx; % number of space points
Nt = T/dt;     % number of time points

x = linspace(a,b,Nx+1); % all space points
t = linspace(0,T,Nt+1); % all time points

% make vector with all combinations possible of the three parameter arrays
C = {gamma,h,D_d};
D = C;
[D{:}] = ndgrid(C{:});
total_combinations = cell2mat(cellfun(@(m)m(:),D,'uni',0));

mean_solns = zeros(size(total_combinations,1),length(t),4);
peak_offsets = zeros(3,size(total_combinations,1));
mean_period_length = zeros(4,size(total_combinations,1));
is_parallel = zeros(size(total_combinations,1),1);

% Start a new parallel pool one if one does not exist
pool = gcp();
if isempty(pool)
    pool = parpool('Processes');
end

% set parfor Options
if fixed == 1
    opts = parforOptions(pool', 'RangePartitionMethod', 'fixed', ...
    'SubrangeSize', 200);
else
    opts = parforOptions(pool', 'RangePartitionMethod', 'auto');
end

% solve PDE system for all given parameter values in parallel
parfor(i = 1:size(total_combinations,1),opts)
    
    s = 0; % symmetry constant for pdepe
    parameters = [alpha_d, alpha_m, alpha_p, alpha_n, ...
        mu_d, mu_m, mu_p, mu_n, total_combinations(i,3), ...
        total_combinations(i,2), total_combinations(i,1)];

    % solve PDE for given parameter values
    soln = pdepe(s,@(t,x,u,DuDx) tissue_pde(t,x,u,DuDx,parameters),@homtissueics,@tissuebcs,x,t);
    d = soln(:,:,1); % Dll1 solution
    m = soln(:,:,2); % Hes1 mRNA solution
    p = soln(:,:,3); % Hes1 protein solution
    n = soln(:,:,4); % Ngn2 solution

    % mean solution for each molecule (this should not differ significantly
    % from individual solutions at each space point since the initial
    % conditions are homogeneous)
    mean_solns(i,:,:) = [mean(d,2), mean(m,2), mean(p,2), mean(n,2)];

    % find local maxima for each molecule
    [~,d_index] = findpeaks(d(:,end));
    [~,m_index] = findpeaks(m(:,end));
    [~,p_index] = findpeaks(p(:,end));
    [~,n_index] = findpeaks(n(:,end));

    % check first value for all molecules to check for a maximum
    if d(1,end) > d(2,end)
        d_index = [1; d_index];
    end

    if m(1,end) > m(2,end)
        m_index = [1; m_index];
    end

    if p(1,end) > p(2,end)
        p_index = [1; p_index];
    end

    if n(1,end) > n(2,end)
        n_index = [1; n_index];
    end

    % find local minima for each molecule
    [~,d_index_min] = findpeaks(-1 * d(:,end));
    [~,m_index_min] = findpeaks(-1 * m(:,end));
    [~,p_index_min] = findpeaks(-1 * p(:,end));
    [~,n_index_min] = findpeaks(-1 * n(:,end));

    % check first value for all molecules to check for a minimum
    if d(1,end) < d(2,end)
        d_index_min = [1; d_index_min];
    end

    if m(1,end) < m(2,end)
        m_index_min = [1; m_index_min];
    end

    if p(1,end) < p(2,end)
        p_index_min = [1; p_index_min];
    end

    if n(1,end) < n(2,end)
        n_index_min = [1; n_index_min];
    end

    % check the amplitude of the last 5 oscillations/maximum-minimum pairs
    % found to determine if the oscillations are stable or not
    if length(d_index) > 5 && length(d_index_min) > 5
    % last 5 maxima of Dll1
    d_last_5_max = d(d_index(end-5:end),end);

    % last 5 minima of Dll1
    d_last_5_min = d(d_index_min(end-5:end),end);

    % mean peak prominence for the last 5 peaks
    mean_peak_prom = mean(d_last_5_max - d_last_5_min);

    else
        mean_peak_prom = 0;
    end

    % mean level of Dll1 from time t=360 minutes to T=2000 minutes
    mean_d_after_start = mean(d(3601:end,end));

    % compare same number of peaks for all variables
    min_peaks = min([length(d_index), length(p_index), ...
        length(m_index), length(n_index)]);

    % check that oscillations are stable
    % if the mean difference between the last 5 minima and maxima of Dll1
    % is less than 5% the mean value of Dll1 after the initial start up
    % period (calculated after t=360 minutes), then the oscillations are
    % considered unstable (or if there are too few maxima in the beginning
    % - threshold of maxima set to 5)
    if (mean_peak_prom > 0.05*mean_d_after_start) && (min_peaks > 5)

        % truncate no. of indices to the shortest vector length
        d_index = d_index(1:min_peaks);
        p_index = p_index(1:min_peaks);
        m_index = m_index(1:min_peaks);
        n_index = n_index(1:min_peaks);

        d_period = zeros(length(d_index)-1,1);
        m_period = zeros(length(d_index)-1,1);
        p_period = zeros(length(d_index)-1,1);
        n_period = zeros(length(d_index)-1,1);

        % find period of oscillation for all variables (in minutes)
        for w=2:length(d_period)
            d_period(w-1) = d_index(w)- d_index(w-1);
            m_period(w-1) = m_index(w)- m_index(w-1);
            p_period(w-1) = p_index(w)- p_index(w-1);
            n_period(w-1) = n_index(w)- n_index(w-1);
        end

        mean_period_length(:,i) = t(round(mean([d_period, m_period, p_period, n_period])));

        % find offset of local maxima between Hes1 mRNA and protein
        peak_offset_m_p = sum(abs(m_index-p_index))/min_peaks;
        % find offset of local maxima between Hes1 protein and Dll1
        peak_offset_d_p = sum(abs(d_index-p_index))/min_peaks;
        % find offset of local maxima between Hes1 protein and Ngn2
        peak_offset_n_p = sum(abs(n_index-p_index))/min_peaks;

        % offsets between oscillations of different molecules (in minutes)
        peak_offsets(:,i) = t(round([peak_offset_m_p, peak_offset_d_p, peak_offset_n_p]));

        % plot time behaviour of system averaged over all space points
        figs1 = figure('visible','off');
        plot(t, mean(d,2), '-g', 'Linewidth', 4, 'Displayname', 'Dll1')
        hold on
        plot(t, mean(m,2), '-r', 'Linewidth', 4, 'Displayname', 'Hes1 mRNA')
        hold on
        plot(t, mean(p,2), '-b', 'Linewidth', 4, 'Displayname', 'Hes1 protein')
        hold on
        plot(t, mean(n,2), '-k', 'Linewidth', 4, 'Displayname', 'Ngn2')
        hold off
        legend('Fontsize', 13)
        xlabel('time (min)', 'Fontsize', 16)
        ylabel('expression', 'Fontsize', 16)
        xticks(0:120:T);
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'Fontsize',13)

        % plot behaviour for all space points
        figs2 = figure('visible','off');
        steps = 180;
        subplot(2,2,1);
        fig1 = pcolor(t,x,transpose(p));
        set(fig1, 'EdgeColor', 'none');
        title('A', 'Fontsize', 16)
        xlabel('time (minutes)', 'Fontsize', 16)
        ylabel('x (\mum)', 'Fontsize', 16)
        xticks(0:steps:T);
        set(gca,'layer','top')

        subplot(2,2,2);
        fig2 = pcolor(t,x,transpose(m));
        set(fig2, 'EdgeColor', 'none');
        title('B', 'Fontsize', 16)
        xlabel('time (minutes)', 'Fontsize', 16)
        ylabel('x (\mum)', 'Fontsize', 16)
        xticks(0:steps:T);
        set(gca,'layer','top')

        subplot(2,2,3);
        fig3 = pcolor(t,x,transpose(d));
        set(fig3, 'EdgeColor', 'none');
        title('C', 'Fontsize', 16)
        xlabel('time (minutes)', 'Fontsize', 16)
        ylabel('x (\mum)', 'Fontsize', 16)
        xticks(0:steps:T);
        set(gca,'layer','top')

        subplot(2,2,4);
        fig4 = pcolor(t,x,transpose(n));
        set(fig4, 'EdgeColor', 'none');
        title('D', 'Fontsize', 16)
        xlabel('time (minutes)', 'Fontsize', 16)
        ylabel('x (\mum)', 'Fontsize', 16)
        xticks(0:steps:T);
        set(gca,'layer','top')

        hp4 = get(subplot(2,2,4),'Position');
        colorbar('Position', [hp4(1)+hp4(3)+0.025  hp4(2)  0.025  hp4(2)+hp4(3)*2.1])

        % save figures
        path1 = './Figures';
        filename1 = ['parallel_average_index_' num2str(i) '_D_d_ ' num2str(total_combinations(i,3)) '_h_' num2str(total_combinations(i,2)) '_gamma_' num2str(total_combinations(i,1)) '.png'];
        filename2 = ['parallel_index_' num2str(i) '_D_d_ ' num2str(total_combinations(i,3)) '_h_' num2str(total_combinations(i,2)) '_gamma_' num2str(total_combinations(i,1)) '.png'];
        saveas(figs1,fullfile(path1,filename1))
        saveas(figs2,fullfile(path1,filename2))
        close(figs1)
        close(figs2)

        % otherwise mean period length is 0 as there are no stable
        % oscillations and, thus, no average period length. The same applies to
        % the peak offset
    else
        mean_period_length(:,i) = [0, 0, 0, 0];
        peak_offsets(:,i) = [0,0,0];
    end

    % check that code actually runs in parallel
    is_parallel(i) = parallel.internal.pool.isPoolThreadWorker||~isempty(getCurrentJob);
end

% save results
path2 = './Results';
if fixed == 1
    matrixname = ['parallel_results_' num2str(length(total_combinations)) '_specified.mat'];
else
    matrixname = ['parallel_results_' num2str(length(total_combinations)) '_auto.mat'];
end
save(fullfile(path2,matrixname), 'total_combinations', 'mean_period_length', 'mean_solns', 'peak_offsets');
tEnd = toc