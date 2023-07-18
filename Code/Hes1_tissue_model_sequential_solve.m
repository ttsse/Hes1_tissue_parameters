function [mean_solns, peak_offsets, mean_period_length] = Hes1_tissue_model_sequential_solve()

clear all;
close all;
format long;

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

% parameters used in parameter sweep
D_d = linspace(0.001,0.1,2);
h = linspace(1,9,2);
gamma = linspace(1,9,2);
total_combinations = length(D_d)*length(h)*length(gamma);

mean_solns = zeros(total_combinations,length(t),4);
peak_offsets = zeros(3,total_combinations);
mean_period_length = zeros(4,total_combinations);

index = 1;

for i = 1:length(D_d)

    for j = 1:length(h)

        for k = 1:length(gamma)
            s = 0;
            parameters = [alpha_d, alpha_m, alpha_p, alpha_n, ...
                mu_d, mu_m, mu_p, mu_n, D_d(i), h(j), gamma(k)];

            soln = pdepe(s,@(t,x,u,DuDx) tissue_pde(t,x,u,DuDx,parameters),@homtissueics,@tissuebcs,x,t);
            d = soln(:,:,1); % Dll1 solution
            m = soln(:,:,2); % Hes1 mRNA solution
            p = soln(:,:,3); % Hes1 protein solution
            n = soln(:,:,4); % Ngn2 solution

            mean_solns(index,:,:) = [mean(d,2), mean(m,2), mean(p,2), mean(n,2)];


            % plot time behaviour of system averaged over all space points
            figs1 = figure('visible','off');
            plot(t, mean_solns(i,:,1), '-g', 'Linewidth', 4, 'Displayname', 'Dll1')
            hold on
            plot(t, mean_solns(i,:,2), '-r', 'Linewidth', 4, 'Displayname', 'Hes1 mRNA')
            hold on
            plot(t, mean_solns(i,:,3), '-b', 'Linewidth', 4, 'Displayname', 'Hes1 protein')
            hold on
            plot(t, mean_solns(i,:,4), '-k', 'Linewidth', 4, 'Displayname', 'Ngn2')
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
            path = './Figures';
            filename1 = ['sequential_average_index_' num2str(index) '_D_d_ ' num2str(D_d(i)) '_h_' num2str(h(j)) '_gamma_' num2str(gamma(k)) '.png'];
            filename2 = ['sequential_index_' num2str(index) '_D_d_ ' num2str(D_d(i)) '_h_' num2str(h(j)) '_gamma_' num2str(gamma(k)) '.png'];
            saveas(figs1,fullfile(path,filename1))
            saveas(figs2,fullfile(path,filename2))
            close(figs1)
            close(figs2)

            % find offset of peaks of oscillations
            d_index = [];
            m_index = [];
            p_index = [];
            n_index = [];

            % special case for N: check if the local maximum is at the first time step
            if n(1,end)>n(2,end)
                n_index(end+1) = t(1);
            end

            for v=2:(length(t)-1)

                % find local maxima of Dll1
                if d(v,end)>d(v-1,end) && d(v,end)>d(v+1,end)
                    d_index = [d_index,t(v)];
                end
                % find local maxima of Hes1 mRNA
                if m(v,end)>m(v-1,end) && m(v,end)>m(v+1,end)
                    m_index = [m_index,t(v)];
                end
                % find local maxima of Hes1 protein
                if p(v,end)>p(v-1,end) && p(v,end)>p(v+1,end)
                    p_index = [p_index,t(v)];
                end
                % find local maxima of Ngn2
                if n(v,end)>n(v-1,end) && n(v,end)>n(v+1,end)
                    n_index = [n_index,t(v)];
                end
            end

            % compare same number of peaks for all variables
            min_peaks = min([length(d_index), length(m_index), ...
                length(p_index), length(n_index)]);
            d_index = d_index(1:(end-(length(d_index)-min_peaks)));
            m_index = m_index(1:(end-(length(m_index)-min_peaks)));
            p_index = p_index(1:(end-(length(p_index)-min_peaks)));
            n_index = n_index(1:(end-(length(n_index)-min_peaks)));

            d_period = zeros(length(d_index)-1,1);
            m_period = zeros(length(d_index)-1,1);
            p_period = zeros(length(d_index)-1,1);
            n_period = zeros(length(d_index)-1,1);

            % find period of oscillation for all variables
            for w=2:length(d_period)
                d_period(w-1) = d_index(w)- d_index(w-1);
                m_period(w-1) = m_index(w)- m_index(w-1);
                p_period(w-1) = p_index(w)- p_index(w-1);
                n_period(w-1) = n_index(w)- n_index(w-1);
            end

            mean_period_length(:,index) = [mean(d_period), mean(m_period), mean(p_period), mean(n_period)];

            % find mean offset of maxima between Hes1 mRNA and protein
            peak_offset_m_p = mean(abs(m_index-p_index));
            % find mean offset of maxima between Hes1 protein and Dll1
            peak_offset_d_p = mean(abs(d_index-p_index));
            % find mean offset of maxima between Hes1 protein and Ngn2
            peak_offset_n_p = mean(abs(n_index-p_index));

            peak_offsets(:,index) = [peak_offset_m_p, peak_offset_d_p, peak_offset_n_p];
            index = index + 1;

        end
    end
end