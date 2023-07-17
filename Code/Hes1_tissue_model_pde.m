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
D_d = 0.003;
h = 1;
gamma = 8;
parameters = [alpha_d, alpha_m, alpha_p, alpha_n, ...
    mu_d, mu_m, mu_p, mu_n, D_d, h, gamma];

% define parameters for numerical solution of x and t
s = 0;
a = 0;
b = 5.12; % full size of fetus (in mm) at E10.5
T = 2000;
dx = 0.1;
dt = 0.1;

Nx = (b-a)/dx; % number of space points
Nt = T/dt;     % number of time points

x = linspace(a,b,Nx+1);
t = linspace(0,T,Nt+1);

soln = pdepe(s,@(t,x,u,DuDx) tissue(t,x,u,DuDx,parameters),@homtissueics,@tissuebcs,x,t);
d = soln(:,:,1); % Dll1 solution
p = soln(:,:,2); % Hes1 protein solution
m = soln(:,:,3); % Hes1 mRNA solution
n = soln(:,:,4); % Ngn2 solution

% find average behaviour (averaged over all space points)
d_average = zeros(length(d),1);
p_average = zeros(length(d),1);
m_average = zeros(length(d),1);
n_average = zeros(length(d),1);

for i=1:length(d)
    d_average(i) = sum(d(i,:))/length(x);
    p_average(i) = sum(p(i,:))/length(x);
    m_average(i) = sum(m(i,:))/length(x);
    n_average(i) = sum(n(i,:))/length(x);
end

% plot time behaviour of system averaged over all space points
figure()
plot(t(1:length(d)),d_average, '-g', 'Linewidth', 4, 'Displayname', 'Dll1')
hold on
plot(t(1:length(d)),p_average, '-b', 'Linewidth', 4, 'Displayname', 'Hes1 protein')
hold on
plot(t(1:length(d)),m_average, '-r', 'Linewidth', 4, 'Displayname', 'Hes1 mRNA')
hold on
plot(t(1:length(d)),n_average, '-k', 'Linewidth', 4, 'Displayname', 'Ngn2')
hold off
legend('Fontsize', 13)
xlabel('time (min)', 'Fontsize', 16)
ylabel('expression', 'Fontsize', 16)
xticks(0:120:T);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'Fontsize',13)

% plot behaviour for all space points
figure()
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

    % find offset of peaks of oscillations
    d_index = [];
    p_index = [];
    m_index = [];
    n_index = [];

    if n(1,end)>n(2,end)
        n_index(end+1) = t(1);
    end

    for i=2:(length(t)-1)

        % find local maxima of Dll1
        if d(i,end)>d(i-1,end) && d(i,end)>d(i+1,end)
            d_index = [d_index,t(i)];
        end
        % find local maxima of Hes1 protein
        if p(i,end)>p(i-1,end) && p(i,end)>p(i+1,end)
            p_index = [p_index,t(i)];
        end
        % find local maxima of Hes1 mRNA
        if m(i,end)>m(i-1,end) && m(i,end)>m(i+1,end)
            m_index = [m_index,t(i)];
        end
        % find local maxima of Ngn2
        if n(i,end)>n(i-1,end) && n(i,end)>n(i+1,end)
            n_index = [n_index,t(i)];
        end
    end

    % compare same number of peaks for all variables
    min_peaks = min([length(d_index), length(p_index), ...
        length(m_index), length(n_index)]);
    d_index = d_index(1:(end-(length(d_index)-min_peaks)));
    p_index = p_index(1:(end-(length(p_index)-min_peaks)));
    m_index = m_index(1:(end-(length(m_index)-min_peaks)));
    n_index = n_index(1:(end-(length(n_index)-min_peaks)));

    d_period = zeros(length(d_index)-1,1);
    p_period = zeros(length(d_index)-1,1);
    m_period = zeros(length(d_index)-1,1);
    n_period = zeros(length(d_index)-1,1);

    % find period of oscillation for all variables
    for i=2:length(d_period)
        d_period(i-1) = d_index(i)- d_index(i-1);
        m_period(i-1) = m_index(i)- m_index(i-1);
        p_period(i-1) = p_index(i)- p_index(i-1);
        n_period(i-1) = n_index(i)- n_index(i-1);
    end

    avg_d_period = sum(d_period)/length(d_period)
    avg_m_period = sum(m_period)/length(d_period)
    avg_p_period = sum(p_period)/length(d_period)
    avg_n_period = sum(n_period)/length(d_period)


    % find offset of local maxima between Hes1 mRNA and protein
    peak_offset_m_p = sum(abs(m_index-p_index))/length(m_index)
    % find offset of local maxima between Hes1 protein and Dll1
    peak_offset_d_p = sum(abs(d_index-p_index))/length(d_index)
    % find offset of local maxima between Hes1 protein and Ngn2
    peak_offset_n_p = sum(abs(n_index-p_index))/length(n_index)

% implementation of the tissue PDE model
function [c,f,s] = tissue(x,t,u,DuDx,parameters)
c = [1; 1; 1; 1];
f = [parameters(9); 0; 0; 0] .* DuDx; % diffusion
% System equations correlating to molecules u = [D, P, M, N]
% parameters = [alpha_d, alpha_m, alpha_p, alpha_n, ...
%    mu_d, mu_m, mu_p, mu_n, D_d, h, gamma];
D = parameters(1) * u(4) - parameters(5) * u(1);
P = parameters(3) * u(3) - parameters(7) * u(2);
M = parameters(2) * u(1) /(1 + u(2)^parameters(10)) - parameters(6) * u(3);
N = parameters(4) / (1 + u(2)^parameters(11)) - parameters(8) * u(4);
s = [D; P; M; N];
end

% homogeneous initial conditions
function u0 = homtissueics(x)
u0 = [1; 2; 2; 1];
end

% boundary conditions (zero flux boundary condition)
function [pl,ql,pr,qr] = tissuebcs(xl,ul,xr,ur,t)
pl = [0; 0; 0; 0];
ql = [1; 1; 1; 1];
pr = [0; 0; 0; 0];
qr = [1; 1; 1; 1];
end