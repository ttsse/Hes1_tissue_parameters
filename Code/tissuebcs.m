% boundary conditions (zero flux boundary condition)
function [pl,ql,pr,qr] = tissuebcs(xl,ul,xr,ur,t)
pl = [0; 0; 0; 0];
ql = [1; 1; 1; 1];
pr = [0; 0; 0; 0];
qr = [1; 1; 1; 1];
end