% implementation of the tissue PDE model
function [c,f,s] = tissue_pde(x,t,u,DuDx,parameters)
% System equations correlating to molecules u = [D, M, P, N]
% parameters = [alpha_d, alpha_m, alpha_p, alpha_n, ...
%    mu_d, mu_m, mu_p, mu_n, D_d, h, gamma];
c = [1; 1; 1; 1];
f = [parameters(9); 0; 0; 0] .* DuDx; % diffusion
D = parameters(1) * u(4) - parameters(5) * u(1);
M = parameters(2) * u(1) /(1 + u(3)^parameters(10)) - parameters(6) * u(2);
P = parameters(3) * u(2) - parameters(7) * u(3);
N = parameters(4) / (1 + u(3)^parameters(11)) - parameters(8) * u(4);
s = [D; M; P; N];
end